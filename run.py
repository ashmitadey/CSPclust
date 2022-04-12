#%%
import scanpy as sc
from CSPclust import process,getpathways,compute_clusters
import warnings
warnings.filterwarnings("ignore")

sc.settings.verbosity = 0             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

#%%
adata = sc.read(r"data\pbmc\pbmc3k.h5ad")

# %%
adata = process(adata)
# %%
gene_pathway_matrix, adatacheck, all_pathways, genesyms, found_pathways = getpathways(adata)
adata = compute_clusters(gene_pathway_matrix, adatacheck, all_pathways, genesyms, found_pathways)
# %%
