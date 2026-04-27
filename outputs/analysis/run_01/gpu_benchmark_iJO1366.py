"""GPU benchmark with iJO1366-scale graph"""
import torch
from torch_geometric.nn import HGTConv
from torch_geometric.data import HeteroData
import time
import random

random.seed(42)
torch.manual_seed(42)

# iJO1366 scale: 1,367 genes, 2,583 reactions, 1,774 metabolites
data = HeteroData()
n_met, n_rxn, n_gene = 1774, 2583, 1367
data['metabolite'].x = torch.randn(n_met, 64)
data['reaction'].x = torch.randn(n_rxn, 64)
data['gene'].x = torch.randn(n_gene, 64)

# Generate edges proportional to real model
n_met_rxn = 6000
n_gene_rxn = 4000
edges_met_rxn = [[random.randint(0, n_met-1), random.randint(0, n_rxn-1)] for _ in range(n_met_rxn)]
edge_idx = torch.tensor(edges_met_rxn, dtype=torch.long).t().contiguous()
data['metabolite', 'consumes', 'reaction'].edge_index = edge_idx
data['reaction', 'rev_consumes', 'metabolite'].edge_index = torch.stack([edge_idx[1], edge_idx[0]])

edges_gene_rxn = [[random.randint(0, n_gene-1), random.randint(0, n_rxn-1)] for _ in range(n_gene_rxn)]
edge_idx2 = torch.tensor(edges_gene_rxn, dtype=torch.long).t().contiguous()
data['gene', 'regulates', 'reaction'].edge_index = edge_idx2
data['reaction', 'rev_regulates', 'gene'].edge_index = torch.stack([edge_idx2[1], edge_idx2[0]])

node_types = ['metabolite', 'reaction', 'gene']
metadata = (
    node_types,
    [
        ('metabolite', 'consumes', 'reaction'),
        ('reaction', 'rev_consumes', 'metabolite'),
        ('gene', 'regulates', 'reaction'),
        ('reaction', 'rev_regulates', 'gene'),
    ]
)

# GPU benchmark
conv = HGTConv(64, 64, metadata, heads=4).cuda()
data_gpu = data.to('cuda')
h_dict = {nt: data_gpu[nt].x for nt in node_types}
edge_dict = {k: data_gpu[k].edge_index for k in data_gpu.edge_types}

for _ in range(5):
    out = conv(h_dict, edge_dict)
torch.cuda.synchronize()

torch.cuda.reset_peak_memory_stats()
t0 = time.time()
for _ in range(50):
    out = conv(h_dict, edge_dict)
torch.cuda.synchronize()
t1 = time.time()
gpu_ms = (t1 - t0) / 50 * 1000
gpu_mem = torch.cuda.max_memory_allocated() / 1024 / 1024
print(f"iJO1366 HGTConv forward GPU: {gpu_ms:.2f} ms/call")
print(f"iJO1366 GPU peak memory: {gpu_mem:.1f} MB")

# CPU benchmark
conv_cpu = HGTConv(64, 64, metadata, heads=4)
h_cpu = {nt: data[nt].x for nt in node_types}
edge_cpu = {k: data[k].edge_index for k in data.edge_types}

for _ in range(5):
    out = conv_cpu(h_cpu, edge_cpu)

t0 = time.time()
for _ in range(50):
    out = conv_cpu(h_cpu, edge_cpu)
t1 = time.time()
cpu_ms = (t1 - t0) / 50 * 1000
print(f"iJO1366 HGTConv forward CPU: {cpu_ms:.2f} ms/call")
print(f"GPU speedup: {cpu_ms / gpu_ms:.1f}x")

# Full training estimate
fwd_per_epoch = 5000 * 2  # 5000 samples, fwd+bwd
total_gpu = fwd_per_epoch * 100 * gpu_ms / 1000
total_cpu = fwd_per_epoch * 100 * cpu_ms / 1000
print(f"\nFull training (100 epochs, 5000 samples):")
print(f"  GPU: {total_gpu:.1f} s ({total_gpu/60:.1f} min)")
print(f"  CPU: {total_cpu:.1f} s ({total_cpu/60:.1f} min)")
