#%%
import scanpy as sc
from CSPclust import process,getpathways,compute_clusters

sc.settings.verbosity = 0             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

#%%
adata = sc.read(r"data\liver\liver.h5ad")

# %%
adata = process(adata)
# %%
gene_pathway_matrix, adatacheck, all_pathways, all_genes = getpathways(adata)
# %%
adata = compute_clusters(gene_pathway_matrix, adatacheck,  all_pathways, all_genes)
# %%
