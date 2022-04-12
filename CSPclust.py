import numpy as np
import pandas as pd
import scanpy as sc
import requests, json

sc.settings.verbosity = 0             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


def process(adata):
    adata.var_names_make_unique()
    sc.pl.highest_expr_genes(adata, n_top=20, save='.pdf')

    # Filtering

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        jitter=0.4, multi_panel=True, save='.pdf')
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)   
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.pl.highly_variable_genes(adata, save='.pdf')
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    return adata

def getpathways(adata):
    adatacheck = adata.to_df()
    genesyms = list(adatacheck.columns)
    len(genesyms)
    column_names = ["Gene", "Pathway_Name", "adjusted_p_value"]
    pathwaydf = pd.DataFrame(columns = column_names)
    pathwaydf['Gene'] = genesyms
    pathwaydf.head()
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(genesyms)
    description = 'My Gene List'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')
    data = json.loads(response.text)    
    print(data)
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    
    user_list_id = data['userListId'] # Enter this ID from the output of previous cell

    gene_set_library = 'KEGG_2016'
    response = requests.get(
    ENRICHR_URL + query_string % (user_list_id, gene_set_library))

    if not response.ok:
        raise Exception('Error fetching enrichment results')

    data = json.loads(response.text)
    for result in data['KEGG_2016']:
        pathway_name = result[1].split("Homo sapiens")[0]
        adjusted_p_value = result[6]
        genelist = result[5]
        for gene in genelist:
            if (len(pathwaydf['Gene'][pathwaydf['Gene'] == gene].index.tolist()) >0):
                rowNum = pathwaydf['Gene'][pathwaydf['Gene'] == gene].index.tolist()[0]
                pathwaydf.loc[pathwaydf.index[rowNum], 'Pathway_Name'] = pathway_name
                pathwaydf.loc[pathwaydf.index[rowNum], 'adjusted_p_value'] = adjusted_p_value
    found_pathways = pathwaydf[~pathwaydf['Pathway_Name'].isnull()]
    not_found_pathways = pathwaydf[pathwaydf['Pathway_Name'].isnull()]        
    len(found_pathways['Pathway_Name'].unique())
    all_pathways = found_pathways['Pathway_Name'].unique()
    all_genes = found_pathways['Gene'].unique()
    print(len(all_pathways))
    print(len(all_genes))
    final_matrix = pd.DataFrame(index=genesyms, columns=all_pathways)
    for index,row in pathwaydf.iterrows():
        rowIndex = row['Gene']
        colIndex = row['Pathway_Name']
        cellValue = row['adjusted_p_value']
        newRow = []
        for pathway in all_pathways:
            if colIndex == pathway:
                newRow.append(cellValue)
            else:
                newRow.append(0.0)
        #print("inserted...." + rowIndex)
        #     print(newRow)
        final_matrix.loc[rowIndex] = newRow
    
    final_matrix
    gene_pathway_matrix = final_matrix.transpose()
    return gene_pathway_matrix, adatacheck, all_pathways, genesyms, found_pathways

def compute_clusters(gene_pathway_matrix, adatacheck, all_pathways, genesyms, found_pathways):
    c = np.corrcoef(adatacheck.to_numpy().astype('float64'), gene_pathway_matrix.to_numpy().astype('float64'))
    
    pd.DataFrame(c).to_csv('./corrmatrix.csv')
    
    cell_pathway_matrix = pd.DataFrame(index=list(adatacheck.index.values), columns=all_pathways, dtype=np.float64)
    for i in range(0,cell_pathway_matrix.shape[0]):
        cell_pathway_matrix.iloc[i] = c[i][cell_pathway_matrix.shape[0]:]
    
    adatanew = sc.AnnData(X=cell_pathway_matrix,
                        obs=cell_pathway_matrix.iloc[:,0:0],
                        var=cell_pathway_matrix.iloc[0:0,:].transpose())
    sc.tl.pca(adatanew, svd_solver='arpack')
    sc.pl.pca(adatanew, color='Metabolic pathways ')
    sc.pp.neighbors(adatanew, n_neighbors=5, n_pcs=45)
    sc.tl.umap(adatanew)
    sc.tl.leiden(adatanew)
    sc.set_figure_params(figsize=[4,4])
    sc.pl.umap(adatanew, color=['leiden'], save='.pdf')
    len(adatanew.obs['leiden'].cat.categories)
    sc.tl.rank_genes_groups(adatanew, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adatanew, n_genes=10, sharey=False, save='.pdf')
    sc.tl.rank_genes_groups(adatanew, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adatanew, n_genes=10, sharey=False, save='.pdf')
    adatanew.obs['leiden'].to_csv('./clusters_cells.csv')
    maxIndices = adatacheck.idxmax()
    matnew = pd.DataFrame(columns=['Gene','Cell','Cluster','Expression_Value'])

    gene_cluster = []

    for gene in genesyms:
        newRow = {'Gene': gene, 'Cell': maxIndices[gene], 'Cluster': adatanew.obs['leiden'][maxIndices[gene]], 'Expression_Value': adatacheck.loc[maxIndices[gene]][gene]}
        matnew = matnew.append(newRow, ignore_index = True)
    maxIndices2 = cell_pathway_matrix.idxmax()
    matnew2 = pd.DataFrame(columns=['Pathway','Cell','Cluster','Genes','Correlation_Value'])

    # gene_cluster = []

    for pathway in list(cell_pathway_matrix.columns):
        newRow = {'Pathway': pathway, 'Cell': maxIndices2[pathway], 'Cluster': adatanew.obs['leiden'][maxIndices2[pathway]], 'Genes': found_pathways[found_pathways['Pathway_Name'] == pathway]['Gene'].values.tolist(), 'Correlation_Value': cell_pathway_matrix.loc[maxIndices2[pathway]][pathway]}
        matnew2 = matnew2.append(newRow, ignore_index = True)
    matnew2 = matnew2.sort_values(by = ['Cluster', 'Correlation_Value'], ascending = [True, False])
    final_df2 = pd.DataFrame()
    # topn = 5

    for i in range(0,len(adatanew.obs['leiden'].cat.categories)):
        final_df2[f'cluster_{i}_Pathways'] = pd.Series(matnew2[matnew2['Cluster'] == str(i)]['Pathway'].tolist())
        final_df2[f'cluster_{i}_Genes'] = pd.Series(matnew2[matnew2['Cluster'] == str(i)]['Genes'].tolist())

    final_df2.to_csv('./pathway_gene_ranks.csv')
    sc.set_figure_params(figsize=[10,6])
    for i in range(0,5):
        pathwaylisttemp = pd.Series(matnew2[matnew2['Cluster'] == str(i)]['Pathway'][:8].tolist())
    #for ipathway in pathwaylisttemp:
        sc.pl.violin(adatanew, pathwaylisttemp, groupby='leiden')
    # sc.pl.violin(adatanew, ['HIF-1 signaling pathway  '], groupby='leiden', save='.pdf')
    return adatanew