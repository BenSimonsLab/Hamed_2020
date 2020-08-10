import numpy as np
import pandas as pd
import scanpy as sc

from pyscenic.rss import regulon_specificity_scores


def scenic_enrichment(fileID, leiden_res):
    file_path = "data/brain-development_"+fileID+".h5ad"
    adata = sc.read_h5ad(file_path)
    auc_mtx = pd.read_csv('data/scenic/auc_mtx_'+fileID+'.csv.gz', index_col=0)
    
    rss_cellType = regulon_specificity_scores(auc_mtx, adata.obs[leiden_res][auc_mtx.index.values].values)
    
    rss_cellType.to_csv("data/scenic/rss_"+fileID+"_"+leiden_res+".csv.gz")
    
    return



scenic_enrichment('normal_subset_ependymal_cleaned', 'leiden_0_17')
