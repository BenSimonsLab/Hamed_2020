#!/usr/bin/env python3

#####################
## Scanpy Pipeline ##
#####################

import numpy as np
import pandas as pd
import scanpy as sc


def run_adata_processing(adata, fileID):
    # expects cleaned and scaled adata as input
    
    # setup directories
    import os
    if not os.path.exists("output/"+fileID+"/"):
        os.makedirs("output/"+fileID+"/")
        
    sc.settings.figdir = 'output/'+fileID+'/'
    
    # HVGs
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var['highly_variable']].copy()
    # regress counts
    sc.pp.regress_out(adata, ['n_counts'])
    sc.pp.scale(adata, max_value=10)
    # PCA
    sc.tl.pca(adata)
    # neighbourhood graph
    sc.pp.neighbors(adata, n_neighbors=15)
    
    # UMAP
    sc.tl.umap(adata)
    
    # Leiden clustering
    list_leiden_res = [0.05, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0]
    
    for leiden_res in list_leiden_res:
        sc.tl.leiden(adata, resolution=leiden_res, key_added='leiden_'+str(leiden_res).replace(".", "_"))
        adata.obs['annot_leiden_'+str(leiden_res).replace(".", "_")] = adata.obs['leiden_'+str(leiden_res).replace(".", "_")]
        
    # save results
    adata.write('data/brain-regeneration_'+fileID+'.h5ad')
    return


def get_merged_adata(sampleIDs):
    list_adata = []
    
    for sampleID in sampleIDs:
        # print(sampleID)
        adata_temp = sc.read_10x_h5('data/samples/' + sampleID + '/outs/filtered_feature_bc_matrix.h5')
        adata_temp.var_names_make_unique()
        list_adata.append(adata_temp)
        
    adata_merged = list_adata[0].concatenate(list_adata[1:],
                                             batch_key = 'sampleID',
                                             batch_categories = sampleIDs)
    
    adata_merged.obs['cellID'] = [cellID.replace('-1-', '-') for cellID in adata_merged.obs.index.tolist()]
    adata_merged.obs.index = pd.Index(adata_merged.obs['cellID'])
    
    return(adata_merged)


def run_pipe(sampleIDs, fileID):
    adata = get_merged_adata(sampleIDs)
    
    # import metadata    
    metadata_exp = pd.read_csv("data/metadata_experiment.tsv", sep='\t')
    
    meta_merged = pd.merge(adata.obs[['cellID', 'sampleID']], metadata_exp, on='sampleID', how='left')
    scrublet_scores = pd.read_csv("data/scrublet-scores/scrublet_scores.csv")
    meta_merged = pd.merge(meta_merged, scrublet_scores, on='cellID', how='left')
    meta_merged = meta_merged.set_index('cellID')
    
    adata.obs = meta_merged
    
    # filter cells based on QC
    metadata_qc = pd.read_csv('data/metadata_QC_pass_normal.csv.gz')
    
    filter_QC = [x in metadata_qc.cellID.values for x in adata.obs_names]
    adata = adata[filter_QC, :].copy()
    
    
    # remove sex specific genes (Zeisel 2018)
    sex_genes = ['Xist', 'Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d']
    filter_sex_genes = [x not in sex_genes for x in adata.var_names]
    adata = adata[:, filter_sex_genes]
    
    # not really needed, but generates n_genes slot
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    mito_genes = adata.var_names.str.startswith('mt-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    
    # Normalise and log1p        
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    
    run_adata_processing(adata, fileID)
    return


# merge samples and run pipeline

run_pipe(sampleIDs=['sampleID_1', 'sampleID_2'], fileID='normal_E12_5')
run_pipe(sampleIDs=['sampleID_3', 'sampleID_4'], fileID='normal_E14_5')
run_pipe(sampleIDs=['sampleID_5', 'sampleID_6'], fileID='normal_E16_5')

run_pipe(sampleIDs=['sampleID_1', 'sampleID_2',
                    'sampleID_3', 'sampleID_4',
                    'sampleID_5', 'sampleID_6',
                    'sampleID_7', 'sampleID_8',
                    'sampleID_9', 'sampleID_10',
                    'sampleID_11', 'sampleID_12',
                    'sampleID_13', 'sampleID_14',
                    'sampleID_15', 'sampleID_16',
                    'sampleID_17', 'sampleID_18',
                    'sampleID_19', 'sampleID_20',
                    'sampleID_21', 'sampleID_22'], fileID='normal')


def run_pipe_subset(main_fileID, subsetID, subset_leidenID, subset_clusterIDs):
        
        fileID = main_fileID+'_'+subsetID
        file_path_main = "data/brain-regeneration_"+main_fileID+".h5ad"
        
        adata = sc.read_h5ad(file_path_main)
        
        filter_sub = [x in subset_clusterIDs for x in adata.obs[subset_leidenID]]
        
        adata_tmp = adata[filter_sub,:]
        adata_sub = sc.AnnData(adata_tmp.raw.X, obs=adata_tmp.obs, var=adata_tmp.raw.var)
        
        adata_sub.raw = adata_sub
        
        run_adata_processing(adata_sub, fileID)
        
        return


# split up by cell types and remove debris

run_pipe_subset("normal_E12_5", "subset_RGCs", "leiden_0_1", ['0'])
run_pipe_subset("normal_E14_5", "subset_RGCs", "leiden_0_1", ['2'])
run_pipe_subset("normal_E16_5", "subset_RGCs", "leiden_0_1", ['3'])

run_pipe_subset("normal", "subset_neuronal_glial", "leiden_0_2", ['1', '2', '3', '4', '5', '6', '8', '9', '10', '11', '15', '16', '18'])

run_pipe_subset("normal", "subset_ependymal", "leiden_0_7", ['18'])
run_pipe_subset("normal_subset_ependymal", "cleaned", "leiden_0_6", ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])



def run_pipe_merged_subset(subsetIDs, fileID):
        
        list_adata = []
        
        for subsetID in subsetIDs:
            # print(sampleID)
            file_path_main = "data/brain-regeneration_"+subsetID+".h5ad"
            adata_temp = sc.read_h5ad(file_path_main)
            list_adata.append(adata_temp)
        
        adata_merged = list_adata[0].concatenate(list_adata[1:],
                                                 batch_key = 'subsetID',
                                                 batch_categories = subsetIDs,
                                                 index_unique=None)
                                                 
        adata = sc.AnnData(adata_merged.raw.X, obs=adata_merged.obs, var=adata_merged.raw.var)
        
        adata.raw = adata
        
        run_adata_processing(adata, fileID)
        
        return


# merge any subsets

run_pipe_merged_subset(["normal_E12_5_subset_RGCs",
                        "normal_E14_5_subset_RGCs",
                        "normal_E16_5_subset_RGCs"],
                            "normal_prenatal_RGCs")

