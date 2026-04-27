"""GPU benchmark for HGTConv on ai_env"""
import torch
from torch_geometric.nn import HGTConv
from torch_geometric.data import HeteroData
import time
import random

random.seed(42)
torch.manual_seed(42)

data = HeteroData()
n_met, n_rxn, n_gene = 72, 95, 137
data['metabolite'].x = torch.randn(n_met, 32)
data['reaction'].x = torch.randn(n_rxn, 32)
data['gene'].x = torch.randn(n_gene, 32)

edges_met_rxn = []
for _ in range(200):
    edges_met_rxn.append([random.randint(0, n_met-1), random.randint(0, n_rxn-1)])
edge_idx = torch.tensor(edges_met_rxn, dtype=torch.long).t().contiguous()
data['metabolite', 'consumes', 'reaction'].edge_index = edge_idx
data['reaction', 'rev_consumes', 'metabolite'].edge_index = torch.stack([edge_idx[1], edge_idx[0]])

edges_gene_rxn = []
for _ in range(160):
    edges_gene_rxn.append([random.randint(0, n_gene-1), random.randint(0, n_rxn-1)])
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

conv = HGTConv(32, 32, metadata, heads=4).cuda()
data = data.to('cuda')
h_dict = {nt: data[nt].x for nt in node_types}
edge_dict = {k: data[k].edge_index for k in data.edge_types}

for _ in range(10):
    out = conv(h_dict, edge_dict)
torch.cuda.synchronize()

torch.cuda.reset_peak_memory_stats()
t0 = time.time()
for _ in range(200):
    out = conv(h_dict, edge_dict)
torch.cuda.synchronize()
t1 = time.time()
gpu_ms = (t1 - t0) / 200 * 1000
gpu_mem = torch.cuda.max_memory_allocated() / 1024 / 1024
print(f"HGTConv forward GPU (textbook): {gpu_ms:.2f} ms/call")
print(f"GPU peak memory: {gpu_mem:.1f} MB")

conv_cpu = HGTConv(32, 32, metadata, heads=4)
h_cpu = {nt: data[nt].x.cpu() for nt in node_types}
edge_cpu = {k: data[k].edge_index.cpu() for k in data.edge_types}
for _ in range(10):
    out = conv_cpu(h_cpu, edge_cpu)
t0 = time.time()
for _ in range(200):
    out = conv_cpu(h_cpu, edge_cpu)
t1 = time.time()
cpu_ms = (t1 - t0) / 200 * 1000
print(f"HGTConv forward CPU (textbook): {cpu_ms:.2f} ms/call")
print(f"GPU speedup: {cpu_ms / gpu_ms:.1f}x")

n_met_l, n_rxn_l, n_gene_l = 1774, 2583, 1367
scale = (n_gene_l / n_gene) * (n_rxn_l / n_rxn) * (n_met_l / n_met)
print(f"\n--- iJO1366 estimate ---")
print(f"Node scale factor: {scale:.1f}x")
print(f"Estimated GPU forward: {gpu_ms * scale:.0f} ms")
print(f"Estimated CPU forward: {cpu_ms * scale:.0f} ms")

fwd_per_epoch = 837 * 2
total_gpu = fwd_per_epoch * 100 * gpu_ms / 1000
total_cpu = fwd_per_epoch * 100 * cpu_ms / 1000
print(f"\nFull training (100 epochs, 837 samples):")
print(f"  GPU: {total_gpu:.1f} s ({total_gpu/60:.1f} min)")
print(f"  CPU: {total_cpu:.1f} s ({total_cpu/60:.1f} min)")
