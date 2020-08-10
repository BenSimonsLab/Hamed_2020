####################
## Add annotation ##
####################

import numpy as np
import pandas as pd
import scanpy as sc


def add_annotation(fileID, leiden_res, annot_dict, annot_key="annot_"):
    # add annotation based on cluster numbers
    file_path = "data/brain-regeneration_"+fileID+".h5ad"
    adata = sc.read_h5ad(file_path)
    
    adata.obs[annot_key+leiden_res] = adata.obs[leiden_res].values.rename_categories(annot_dict)
    
    # save results
    adata.write('data/brain-regeneration_'+fileID+'.h5ad')
    
    return


# Examples of annotations

# Ependymal cells
annot_ependymal = {'0': 'Adult ependymal cells',
                   '1': 'Juvenile ependymal cells',
                   '2': 'Juvenile ependymal cells [cycling]'}
add_annotation("normal_subset_ependymal_cleaned", "leiden_0_15", annot_ependymal)


# OPCs
annot_OPCs = {'0': 'Juvenile OPCs',
              '1': 'Juvenile OPCs [cycling]',
              '2': 'Adult OPCs'}
add_annotation("normal_subset_OPCs_oligos_cleaned_subset_OPCs", "leiden_0_1", annot_OPCs)



# Neuroblasts
annot_neuroblasts = {'0': 'GE NBs [1]',
                     '1': 'Hippocampal NBs [1]',
                     '2': 'GE NBs [2]',
                     '3': 'Hippocampal NBs [2]',
                     '4': 'Early EmDienNBs',
                     '5': 'Early EmCorNBs',
                     '6': 'Early EmSthNBs',
                     '7': 'Early EmNBs'}
add_annotation("normal_subset_neuroblasts", "leiden_0_15", annot_neuroblasts)


annot_complete = {'0': 'Adult microglia [1]',
                  '1': 'GE NBs [1]',
                  '2': 'Hippocampal NBs [1]',
                  '3': 'Juvenile RG & TAPs',
                  '4': 'Oligodendrocytes [1]',
                  '5': 'Quiescent NSCs [1]',
                  '6': 'Juvenile microglia',
                  '7': 'Embryonic RG',
                  '8': 'Early EmNBs',
                  '9': 'OPCs',
                  '10': 'Gliogenic precursors (APCs)',
                  '11': 'Gliogenic precursors (aNSCs)',
                  '12': 'GE NBs [2]',
                  '13': 'Quiescent NSCs [2]',
                  '14': 'Hippocampal NBs [2]',
                  '15': 'ImStNeurons',
                  '16': 'Adult microglia [2]',
                  '17': 'Ependymal cells',
                  '18': 'Macrophages',
                  '19': 'VLMC',
                  '20': 'T-cells',
                  '21': 'GABAergic INs',
                  '22': 'Oligodendrocytes [2]',
                  '23': 'Erythrocytes',
                  '24': 'Choroid plexus epithelia',
                  '25': 'ImPreMDs',
                  '26': 'Endothelial cells',
                  '27': 'Myeloid-DSCs'}
add_annotation("normal", "leiden_0_52", annot_complete)

