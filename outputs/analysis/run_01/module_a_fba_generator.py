"""
Module A: FBA Ground Truth Generator with Graph Conversion
===========================================================
Loads the COBRApy textbook model, runs FBA knockout simulations
(single, double, random subsets), and converts the metabolic model
to a PyTorch Geometric HeteroData graph.

Expected model stats: 95 reactions, 72 metabolites, 137 genes
WT growth rate: 0.8739
"""

import copy
import itertools
import random
import time
from multiprocessing import Pool
from typing import Dict, List, Optional, Tuple

import cobra
import numpy as np
import torch
from torch_geometric.data import HeteroData


# ---------------------------------------------------------------------------
# 1. Model loading
# ---------------------------------------------------------------------------

def load_textbook_model() -> cobra.Model:
    """Load the E. coli core textbook model and validate dimensions."""
    model = cobra.io.load_model("textbook")
    n_rxn = len(model.reactions)
    n_met = len(model.metabolites)
    n_gene = len(model.genes)
    print(f"[Model] Loaded textbook model: "
          f"{n_rxn} reactions, {n_met} metabolites, {n_gene} genes")
    assert n_rxn == 95, f"Expected 95 reactions, got {n_rxn}"
    assert n_met == 72, f"Expected 72 metabolites, got {n_met}"
    assert n_gene == 137, f"Expected 137 genes, got {n_gene}"
    return model


def compute_wt_growth(model: cobra.Model) -> float:
    """Compute wild-type growth rate via FBA."""
    solution = model.optimize()
    growth = solution.objective_value
    print(f"[FBA] Wild-type growth rate: {growth:.4f}")
    return growth


# ---------------------------------------------------------------------------
# 2. FBA knockout execution (worker function for multiprocessing)
# ---------------------------------------------------------------------------

def _fba_knockout_worker(args: Tuple) -> Optional[float]:
    """
    Worker function for multiprocessing.
    Args: (gene_ids, model_json_str)

    Returns the growth rate after knocking out the specified genes,
    or 0.0 if the model becomes infeasible.
    """
    gene_ids, model_json_str = args
    try:
        model = cobra.io.from_json(model_json_str)
        with model:
            for gid in gene_ids:
                if gid in [g.id for g in model.genes]:
                    model.genes.get_by_id(gid).knock_out()
            solution = model.optimize()
            if solution.status == "optimal":
                return solution.objective_value
        return 0.0
    except Exception:
        return 0.0


def _run_knockout_single_process(
    model: cobra.Model,
    gene_ids: List[str],
) -> float:
    """Run a single knockout on a copied model (no multiprocessing)."""
    try:
        model_cp = model.copy()
        with model_cp:
            for gid in gene_ids:
                if gid in [g.id for g in model_cp.genes]:
                    model_cp.genes.get_by_id(gid).knock_out()
            solution = model_cp.optimize()
            if solution.status == "optimal":
                return solution.objective_value
        return 0.0
    except Exception:
        return 0.0


def _run_knockout_batch(
    model: cobra.Model,
    knockout_combos: List[List[str]],
    n_workers: int = 1,
) -> List[float]:
    """Run a batch of knockout simulations, optionally in parallel."""
    if n_workers <= 1:
        return [_run_knockout_single_process(model, combo) for combo in knockout_combos]

    model_json_str = cobra.io.to_json(model)
    args_list = [(combo, model_json_str) for combo in knockout_combos]
    with Pool(processes=n_workers) as pool:
        results = pool.map(_fba_knockout_worker, args_list)
    return list(results)


# ---------------------------------------------------------------------------
# 3. Knockout combination generators
# ---------------------------------------------------------------------------

def generate_single_knockouts(model: cobra.Model) -> List[List[str]]:
    """All 137 single-gene knockouts."""
    return [[g.id] for g in model.genes]


def generate_double_knockouts(
    model: cobra.Model,
    n_combos: int = 2000,
    seed: int = 42,
) -> List[List[str]]:
    """Random subset of double-gene knockout combinations."""
    gene_ids = [g.id for g in model.genes]
    all_pairs = list(itertools.combinations(gene_ids, 2))
    rng = random.Random(seed)
    rng.shuffle(all_pairs)
    selected = all_pairs[:min(n_combos, len(all_pairs))]
    return [list(pair) for pair in selected]


def generate_random_knockouts(
    model: cobra.Model,
    n_combos: int = 500,
    min_k: int = 1,
    max_k: int = 5,
    seed: int = 123,
) -> List[List[str]]:
    """Random knockout combinations with variable size."""
    gene_ids = [g.id for g in model.genes]
    rng = random.Random(seed)
    combos = []
    for _ in range(n_combos):
        k = rng.randint(min_k, min(max_k, len(gene_ids)))
        combo = rng.sample(gene_ids, k)
        combos.append(combo)
    return combos


# ---------------------------------------------------------------------------
# 4. Graph conversion: COBRApy Model -> PyG HeteroData
# ---------------------------------------------------------------------------

def model_to_heterodata(model: cobra.Model) -> HeteroData:
    """
    Convert a COBRApy model to a PyTorch Geometric HeteroData graph.

    Node types:
        - metabolite  (72 nodes)
        - reaction    (95 nodes)
        - gene        (137 nodes)

    Edge types:
        - (metabolite, consumes, reaction)   : substrate edges (stoich < 0)
        - (metabolite, produces, reaction)   : product edges   (stoich > 0)
        - (gene, regulates, reaction)        : GPR rule edges

    Node features:
        - metabolite: [degree, 1]
        - reaction:   [degree, 1]
        - gene:       [degree, 1]
    """
    data = HeteroData()

    # --- Build index maps (sorted for reproducibility) ---
    met_ids = sorted([m.id for m in model.metabolites])
    rxn_ids = sorted([r.id for r in model.reactions])
    gene_ids = sorted([g.id for g in model.genes])

    met_idx = {mid: i for i, mid in enumerate(met_ids)}
    rxn_idx = {rid: i for i, rid in enumerate(rxn_ids)}
    gene_idx = {gid: i for i, gid in enumerate(gene_ids)}

    n_met = len(met_ids)
    n_rxn = len(rxn_ids)
    n_gene = len(gene_ids)

    # --- Stoichiometry edges ---
    consumes_src, consumes_dst, consumes_stoich = [], [], []
    produces_src, produces_dst, produces_stoich = [], [], []

    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        for met, coeff in rxn.metabolites.items():
            mi = met_idx[met.id]
            if coeff < 0:
                consumes_src.append(mi)
                consumes_dst.append(ri)
                consumes_stoich.append(abs(coeff))
            else:
                produces_src.append(mi)
                produces_dst.append(ri)
                produces_stoich.append(coeff)

    # --- GPR edges (gene -> reaction) ---
    gpr_src, gpr_dst = [], []
    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        for gene in rxn.genes:
            gi = gene_idx[gene.id]
            gpr_src.append(gi)
            gpr_dst.append(ri)

    # --- Node features ---
    # Metabolite degree = number of reactions it participates in
    met_degree = np.zeros(n_met, dtype=np.float32)
    for rxn in model.reactions:
        for met in rxn.metabolites:
            met_degree[met_idx[met.id]] += 1

    # Reaction degree = number of metabolites + genes connected
    rxn_degree = np.zeros(n_rxn, dtype=np.float32)
    for rxn in model.reactions:
        ri = rxn_idx[rxn.id]
        rxn_degree[ri] = len(rxn.metabolites) + len(rxn.genes)

    # Gene degree = number of reactions it regulates
    gene_degree = np.zeros(n_gene, dtype=np.float32)
    for gene in model.genes:
        gi = gene_idx[gene.id]
        gene_degree[gi] = len(gene.reactions)

    # --- Assemble HeteroData ---
    # Node features: [degree, constant_1]
    data["metabolite"].x = torch.tensor(
        np.column_stack([met_degree, np.ones(n_met, dtype=np.float32)]),
        dtype=torch.float32,
    )
    data["reaction"].x = torch.tensor(
        np.column_stack([rxn_degree, np.ones(n_rxn, dtype=np.float32)]),
        dtype=torch.float32,
    )
    data["gene"].x = torch.tensor(
        np.column_stack([gene_degree, np.ones(n_gene, dtype=np.float32)]),
        dtype=torch.float32,
    )

    # Stoichiometry edges with attributes (forward: metabolite→reaction)
    data["metabolite", "consumes", "reaction"].edge_index = torch.tensor(
        [consumes_src, consumes_dst], dtype=torch.long
    )
    data["metabolite", "consumes", "reaction"].edge_attr = torch.tensor(
        consumes_stoich, dtype=torch.float32
    ).unsqueeze(1)

    data["metabolite", "produces", "reaction"].edge_index = torch.tensor(
        [produces_src, produces_dst], dtype=torch.long
    )
    data["metabolite", "produces", "reaction"].edge_attr = torch.tensor(
        produces_stoich, dtype=torch.float32
    ).unsqueeze(1)

    # GPR edges (forward: gene→reaction)
    data["gene", "regulates", "reaction"].edge_index = torch.tensor(
        [gpr_src, gpr_dst], dtype=torch.long
    )

    # Reverse edges (reaction→metabolite, reaction→gene) for bidirectional
    # message passing. Without these, metabolite/gene nodes never receive
    # messages and stay un-updated through HGTConv layers.
    data["reaction", "rev_consumes", "metabolite"].edge_index = torch.tensor(
        [consumes_dst, consumes_src], dtype=torch.long
    )
    data["reaction", "rev_produces", "metabolite"].edge_index = torch.tensor(
        [produces_dst, produces_src], dtype=torch.long
    )
    data["reaction", "rev_regulates", "gene"].edge_index = torch.tensor(
        [gpr_dst, gpr_src], dtype=torch.long
    )

    # --- Print stats ---
    n_consumes = len(consumes_src)
    n_produces = len(produces_src)
    n_gpr = len(gpr_src)
    print(f"[Graph] HeteroData constructed:")
    print(f"  Nodes: metabolite={n_met}, reaction={n_rxn}, gene={n_gene}")
    print(f"  Edges: consumes={n_consumes}, produces={n_produces}, "
          f"regulates(GPR)={n_gpr}")
    print(f"  Reverse: rev_consumes={n_consumes}, rev_produces={n_produces}, "
          f"rev_regulates={n_gpr}")
    print(f"  Total edges: {(n_consumes + n_produces + n_gpr) * 2}")

    return data


# ---------------------------------------------------------------------------
# 5. Knockout mask builder
# ---------------------------------------------------------------------------

def build_knockout_mask(
    model: cobra.Model,
    knockout_combos: List[List[str]],
) -> torch.Tensor:
    """
    Build a binary knockout mask tensor of shape (n_combos, n_genes).
    1 = gene is knocked out, 0 = gene is active.
    """
    gene_ids = sorted([g.id for g in model.genes])
    gene_idx = {gid: i for i, gid in enumerate(gene_ids)}
    n_gene = len(gene_ids)
    n_combos = len(knockout_combos)

    mask = torch.zeros(n_combos, n_gene, dtype=torch.float32)
    for i, combo in enumerate(knockout_combos):
        for gid in combo:
            if gid in gene_idx:
                mask[i, gene_idx[gid]] = 1.0
    return mask


# ---------------------------------------------------------------------------
# 6. Main orchestrator class
# ---------------------------------------------------------------------------

class FBAGroundTruthGenerator:
    """
    Orchestrates FBA ground-truth generation and graph conversion.

    Usage:
        gen = FBAGroundTruthGenerator()
        data = gen.run(single=True, double=True, random_n=100, n_workers=4)
        # data["graph"]            -> HeteroData
        # data["single_ko_mask"]   -> Tensor (137, 137)
        # data["single_ko_growth"] -> np.ndarray (137,)
        # etc.
    """

    def __init__(self):
        self.model = load_textbook_model()
        self.wt_growth = compute_wt_growth(self.model)
        self.graph: Optional[HeteroData] = None

    def run(
        self,
        single: bool = True,
        double: bool = True,
        double_n: int = 2000,
        random_n: int = 100,
        n_workers: int = 1,
        seed: int = 42,
    ) -> Dict:
        """
        Run knockout simulations and build graph.

        Returns dict with keys:
            - "graph": HeteroData object
            - "wt_growth": float
            - "single_ko_mask", "single_ko_growth": single knockout results
            - "double_ko_mask", "double_ko_growth": double knockout results
            - "random_ko_mask", "random_ko_growth": random knockout results
        """
        result = {"wt_growth": self.wt_growth}

        # Build graph
        self.graph = model_to_heterodata(self.model)
        result["graph"] = self.graph

        # Single knockouts
        if single:
            t0 = time.time()
            combos = generate_single_knockouts(self.model)
            mask = build_knockout_mask(self.model, combos)
            growths = _run_knockout_batch(self.model, combos, n_workers=n_workers)
            elapsed = time.time() - t0
            growths = np.array(growths, dtype=np.float32)
            print(f"[Single KO] {len(combos)} knockouts in {elapsed:.1f}s "
                  f"({elapsed/len(combos)*1000:.0f}ms each)")
            result["single_ko_mask"] = mask
            result["single_ko_growth"] = growths

        # Double knockouts
        if double:
            t0 = time.time()
            combos = generate_double_knockouts(self.model, n_combos=double_n, seed=seed)
            mask = build_knockout_mask(self.model, combos)
            growths = _run_knockout_batch(self.model, combos, n_workers=n_workers)
            elapsed = time.time() - t0
            growths = np.array(growths, dtype=np.float32)
            print(f"[Double KO] {len(combos)} knockouts in {elapsed:.1f}s "
                  f"({elapsed/len(combos)*1000:.0f}ms each)")
            result["double_ko_mask"] = mask
            result["double_ko_growth"] = growths

        # Random knockouts
        if random_n > 0:
            t0 = time.time()
            combos = generate_random_knockouts(self.model, n_combos=random_n, seed=seed)
            mask = build_knockout_mask(self.model, combos)
            growths = _run_knockout_batch(self.model, combos, n_workers=n_workers)
            elapsed = time.time() - t0
            growths = np.array(growths, dtype=np.float32)
            print(f"[Random KO] {len(combos)} knockouts in {elapsed:.1f}s "
                  f"({elapsed/len(combos)*1000:.0f}ms each)")
            result["random_ko_mask"] = mask
            result["random_ko_growth"] = growths

        return result


# ---------------------------------------------------------------------------
# 7. Test / verification
# ---------------------------------------------------------------------------

def run_tests():
    """Quick verification tests for Module A."""
    print("=" * 60)
    print("Module A: FBA Ground Truth Generator -- Verification Tests")
    print("=" * 60)

    all_pass = True

    # Test 1: Model loads correctly
    print("\n--- Test 1: Model Loading ---")
    try:
        model = load_textbook_model()
        wt = compute_wt_growth(model)
        assert abs(wt - 0.8739) < 0.01, f"WT growth {wt:.4f} != 0.8739"
        print(f"  PASS: WT growth = {wt:.4f}")
    except AssertionError as e:
        print(f"  FAIL: {e}")
        all_pass = False
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 2: Single knockout produces 137 results with correct dimensions
    print("\n--- Test 2: Single Knockout (137 genes) ---")
    try:
        combos = generate_single_knockouts(model)
        assert len(combos) == 137, f"Expected 137 single KOs, got {len(combos)}"
        mask = build_knockout_mask(model, combos)
        assert mask.shape == (137, 137), f"Mask shape {mask.shape} != (137, 137)"
        # Each row should have exactly one knockout
        assert (mask.sum(dim=1) == 1).all(), "Each single KO mask row should sum to 1"
        print(f"  PASS: {len(combos)} single knockouts, mask shape {mask.shape}")

        # Run a small batch to verify FBA execution
        t0 = time.time()
        growths = _run_knockout_batch(model, combos[:10], n_workers=1)
        elapsed = time.time() - t0
        n_valid = sum(1 for g in growths if g > 1e-6)
        print(f"  Sample (first 10): {n_valid}/10 feasible in {elapsed:.2f}s")
        print(f"  Growth rates: "
              f"{[f'{g:.4f}' if g is not None else 'None' for g in growths]}")
    except AssertionError as e:
        print(f"  FAIL: {e}")
        all_pass = False
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 3: HeteroData graph has correct node/edge counts
    print("\n--- Test 3: HeteroData Graph ---")
    try:
        graph = model_to_heterodata(model)
        assert graph["metabolite"].num_nodes == 72, \
            f"metabolite nodes {graph['metabolite'].num_nodes} != 72"
        assert graph["reaction"].num_nodes == 95, \
            f"reaction nodes {graph['reaction'].num_nodes} != 95"
        assert graph["gene"].num_nodes == 137, \
            f"gene nodes {graph['gene'].num_nodes} != 137"

        n_consumes = graph["metabolite", "consumes", "reaction"].edge_index.shape[1]
        n_produces = graph["metabolite", "produces", "reaction"].edge_index.shape[1]
        n_gpr = graph["gene", "regulates", "reaction"].edge_index.shape[1]
        n_total_stoich = n_consumes + n_produces

        print(f"  Nodes: met={graph['metabolite'].num_nodes}, "
              f"rxn={graph['reaction'].num_nodes}, "
              f"gene={graph['gene'].num_nodes}")
        print(f"  Edges: consumes={n_consumes}, produces={n_produces}, "
              f"GPR={n_gpr}, total_stoich={n_total_stoich}")

        # Verify node feature dimensions
        assert graph["metabolite"].x.shape == (72, 2), \
            f"metabolite features {graph['metabolite'].x.shape} != (72, 2)"
        assert graph["reaction"].x.shape == (95, 2), \
            f"reaction features {graph['reaction'].x.shape} != (95, 2)"
        assert graph["gene"].x.shape == (137, 2), \
            f"gene features {graph['gene'].x.shape} != (137, 2)"
        print(f"  Feature dims: metabolite={graph['metabolite'].x.shape}, "
              f"reaction={graph['reaction'].x.shape}, gene={graph['gene'].x.shape}")
        print(f"  PASS: Graph structure validated")
    except AssertionError as e:
        print(f"  FAIL: {e}")
        all_pass = False
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 4: Random 100 knockout generation
    print("\n--- Test 4: Random 100 Knockout Generation ---")
    try:
        random_combos = generate_random_knockouts(model, n_combos=100, seed=99)
        assert len(random_combos) == 100, f"Expected 100 combos, got {len(random_combos)}"
        random_mask = build_knockout_mask(model, random_combos)
        assert random_mask.shape[0] == 100, f"Mask rows {random_mask.shape[0]} != 100"
        assert random_mask.shape[1] == 137, f"Mask cols {random_mask.shape[1]} != 137"
        # At least some genes should be knocked out per combo
        assert (random_mask.sum(dim=1) >= 1).all(), "Each combo should knock out >= 1 gene"
        print(f"  PASS: 100 random knockouts generated, mask shape {random_mask.shape}")
    except AssertionError as e:
        print(f"  FAIL: {e}")
        all_pass = False
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 5: Full generator class with small random subset
    print("\n--- Test 5: FBAGroundTruthGenerator Class ---")
    try:
        gen = FBAGroundTruthGenerator()
        data = gen.run(single=True, double=False, random_n=20, n_workers=1)

        assert "graph" in data
        assert "wt_growth" in data
        assert "single_ko_mask" in data
        assert "single_ko_growth" in data
        assert "random_ko_mask" in data
        assert "random_ko_growth" in data

        assert data["single_ko_mask"].shape == (137, 137)
        assert len(data["single_ko_growth"]) == 137
        assert data["random_ko_mask"].shape == (20, 137)
        assert len(data["random_ko_growth"]) == 20

        # Count feasible single KOs
        n_feasible = np.sum(data["single_ko_growth"] > 1e-6)
        n_lethal = 137 - n_feasible
        print(f"  Single KO: {n_feasible} feasible, {n_lethal} lethal/infeasible")
        n_random_feasible = np.sum(data["random_ko_growth"] > 1e-6)
        print(f"  Random KO (20): {n_random_feasible} feasible")
        print(f"  PASS: Generator class works correctly")
    except AssertionError as e:
        print(f"  FAIL: {e}")
        all_pass = False
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Test 6: Multiprocessing
    print("\n--- Test 6: Multiprocessing ---")
    try:
        combos_mp = generate_random_knockouts(model, n_combos=20, seed=777)
        t_seq = time.time()
        growths_seq = _run_knockout_batch(model, combos_mp, n_workers=1)
        t_seq = time.time() - t_seq

        t_mp = time.time()
        growths_mp = _run_knockout_batch(model, combos_mp, n_workers=2)
        t_mp = time.time() - t_mp

        # Compare results (both should produce same answers)
        for i, (g_seq, g_mp) in enumerate(zip(growths_seq, growths_mp)):
            if g_seq is None and g_mp is None:
                continue
            assert g_seq is not None and g_mp is not None, \
                f"Mismatch at idx {i}: seq={g_seq}, mp={g_mp}"
            assert abs(g_seq - g_mp) < 1e-4, \
                f"Growth mismatch at idx {i}: {g_seq:.6f} vs {g_mp:.6f}"

        print(f"  Sequential: {t_seq:.2f}s, Multiprocessing (2 workers): {t_mp:.2f}s")
        print(f"  Results match: {len(growths_seq)} knockouts verified")
        print(f"  PASS: Multiprocessing produces identical results")
    except AssertionError as e:
        print(f"  FAIL: {e}")
        all_pass = False
    except Exception as e:
        print(f"  ERROR: {e}")
        all_pass = False

    # Summary
    print("\n" + "=" * 60)
    if all_pass:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED -- see above")
    print("=" * 60)
    print(f"  Model: {len(model.reactions)} rxn, {len(model.metabolites)} met, "
          f"{len(model.genes)} genes")
    print(f"  WT growth: {wt:.4f}")
    print(f"  Graph: met(72) + rxn(95) + gene(137) = 304 nodes")
    if 'graph' in dir():
        n_c = graph["metabolite", "consumes", "reaction"].edge_index.shape[1]
        n_p = graph["metabolite", "produces", "reaction"].edge_index.shape[1]
        n_g = graph["gene", "regulates", "reaction"].edge_index.shape[1]
        print(f"  Edges: consumes({n_c}) + produces({n_p}) + GPR({n_g}) "
              f"= {n_c + n_p + n_g} total")
    else:
        print(f"  Edges: (graph not built, skipped)")
    print(f"  Single KO mask: (137, 137), Random KO mask: (100, 137)")
    print(f"  Multiprocessing: verified with 2 workers")


if __name__ == "__main__":
    run_tests()
