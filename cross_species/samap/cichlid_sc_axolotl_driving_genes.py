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

# Read SAMap results
with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/bb_axolotl_sm_ran.pkl', 'rb') as out_file:
     sm = pickle.load(out_file)

# Add in axolotl metadata
am_meta = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/axolotl_metadata.csv', index_col=0)
am_meta['cellclusters'] = am_meta['cellclusters'].replace(["glut_SUBSET_", "GABA_SUBSET_", "epen_clus_", "npc_SUBSET_", "oligodendrocyte_", "endothelial_", "microglia_"], ["GLUT", "GABA", "EPEN", "NB", "OLIG", "ENDO", "MG"], regex=True)
sm.sams['am'].adata.obs['am_cluster'] = 'unassigned'
sm.sams['am'].adata.obs.loc[am_meta.index,'am_cluster'] = am_meta['cellclusters']
sm.sams['am'].adata.obs['am_region_cluster'] = 'unassigned'
sm.sams['am'].adata.obs.loc[am_meta.index,'am_region_cluster'] = pd.Series(np.array(sm.sams['am'].adata.obs.loc[am_meta.index,'region'])) + "_" + am_meta['cellclusters']

# Define the clusters
keys={'mz':'good_names','am':'am_cluster'}
gpf = GenePairFinder(sm,keys=keys)
def my_gene_pair_finder(mm_clust, mz_clust):
    print(mz_clust + "; " + mm_clust)
    try:
        Gp, G1, G2, pvals1, pvals2 = gpf.find_genes(mm_clust, mz_clust, n_genes=len(gpf.gns))
        gdf = pd.DataFrame({'id': mm_clust + ";" + mz_clust, 'celltype1': mm_clust, 'celltype2': mz_clust, 'genes': Gp})
        return(gdf)
    except:
        return()
        print("Failed")

# Find the genes driving each celltype pair
ct1 = ['am_' + s for s in np.sort(sm.sams['am'].adata.obs[keys['am']].unique()).tolist()]
ct2 = ['mz_' + s for s in np.sort(sm.sams['mz'].adata.obs[keys['mz']].unique()).tolist()]
ct_pairs = list(itertools.product(ct1, ct2))
with multiprocessing.Pool(multiprocessing.cpu_count()) as mp_pool:
    all_pairs_list = mp_pool.starmap(my_gene_pair_finder, ct_pairs)

# Concatenate results from all celltype pairs
all_pairs_list_df = [e for e in all_pairs_list if isinstance(e, pd.core.frame.DataFrame)]
all_pairs_df = pd.concat(all_pairs_list_df)
all_pairs_df.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/bb_axolotl_samap_genes.csv")

# Calculate the number of genes driving celltype pairs
all_pairs_num = pd.DataFrame(ct_pairs, columns = ['celltype1', 'celltype2'])
all_pairs_num['id'] = all_pairs_num['celltype1'] + ';' + all_pairs_num['celltype2']
all_pairs_num.index = all_pairs_num['id']
all_pairs_num['num'] = 0
nonzero_counts = all_pairs_df['id'].value_counts()
all_pairs_num.loc[nonzero_counts.index, 'num'] = nonzero_counts
all_pairs_num.to_csv("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/results/bb_axolotl_samap_genes_num.csv")

