# Imports
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import numpy as np
import pickle
import multiprocessing
import itertools
from functools import reduce
from sklearn.preprocessing import LabelBinarizer

# Helper Functions
def my_mapper(mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent'):
    meta = sm.samap.adata.obs
    all_mz_cluster = meta[mz_col].unique()
    all_mm_cluster = meta[mm_col].unique()
    all_mz_cluster.sort()
    all_mm_cluster.sort()
    lb = LabelBinarizer(sparse_output=True)
    knn_mz  = lb.fit_transform(sm.samap.adata.obs[mz_col]).T.dot(knn)
    knn_sum = pd.DataFrame(lb.fit_transform(sm.samap.adata.obs[mm_col]).T.dot(knn_mz.T).todense())
    mz_counts = meta[mz_col].value_counts().sort_index()
    mm_counts = meta[mm_col].value_counts().sort_index()
    mz_mm_counts = pd.DataFrame(np.array(mz_counts) * np.array(mm_counts).reshape(len(mm_counts), 1))
    knn_mean = knn_sum / mz_mm_counts
    knn_mean.index = all_mm_cluster
    knn_mean.columns = all_mz_cluster
    knn_mean = knn_mean.drop("unassigned", axis=0)
    knn_mean = knn_mean.drop("unassigned", axis=1)
    return(knn_mean)

def perm_mapper(seed_num, mz_col = 'mz_struct_b2_vdc', mm_col = 'mm_ABA_parent'):
    np.random.seed(seed_num)
    sm.samap.adata.obs['mz_perm'] = sm.samap.adata.obs[mz_col]
    sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_perm'] != 'unassigned', 'mz_perm'] = np.random.permutation(sm.samap.adata.obs.loc[sm.samap.adata.obs['mz_perm'] != 'unassigned', 'mz_perm'])
    sm.samap.adata.obs['mm_perm'] = sm.samap.adata.obs[mm_col]
    sm.samap.adata.obs.loc[sm.samap.adata.obs['mm_perm'] != 'unassigned', 'mm_perm'] = np.random.permutation(sm.samap.adata.obs.loc[sm.samap.adata.obs['mm_perm'] != 'unassigned', 'mm_perm'])
    perm_map = my_mapper('mz_perm', 'mm_perm')
    return(perm_map)

def my_p(x):
    return( 1 - (np.count_nonzero(perm_dist<x) / perm_dist.size))

# Main
# Read SAMap results
with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/bb_bird_sm_ran_3.pkl', 'rb') as out_file:
     sm = pickle.load(out_file)

# Find real mapping values
knn = sm.samap.adata.obsp['knn']
nonzero_mask = np.array(knn[knn.nonzero()] > 0)[0]
rows = knn.nonzero()[0][nonzero_mask]
cols = knn.nonzero()[1][nonzero_mask]
knn[rows, cols] = 1

# Add in bird metadata
sm.samap.adata.obs['tg_region_cluster'] = sm.samap.adata.obs['tg_region'] + "_" + sm.samap.adata.obs['tg_cluster_orig2']
sm.samap.adata.obs.loc[sm.samap.adata.obs['tg_region_cluster'] == "unassigned_unassigned", 'tg_region_cluster'] = "unassigned"

real = my_mapper(mz_col = 'mz_good_names', mm_col = 'tg_cluster_orig2')
real.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/bb_bird_cluster_mapping_mine4.csv")

# Permutations
perm_nums = 1000
with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    perm_list = mp_pool.starmap(perm_mapper, zip(list(range(1, perm_nums+1)), ['mz_good_names'] * perm_nums, ['tg_cluster_orig2'] * perm_nums))

# Find permutation p-values
def p_col(x):
    all_col_values = list()
    for i in range(0, len(perm_list)):
        all_col_values.extend(perm_list[i].iloc[:,x].tolist())
    ge_mat = np.greater(np.array([all_col_values]), np.array([real.iloc[:,x]]).T)
    p = ge_mat.sum(axis=1) / (len(perm_list)*real.shape[0])
    return(p)

with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    p_list = mp_pool.map(p_col, range(0, real.shape[1]))

perm_p = pd.DataFrame(p_list, columns = real.index, index = real.columns).T
perm_p.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/bb_bird_cluster_mapping_mine_p4.csv")

