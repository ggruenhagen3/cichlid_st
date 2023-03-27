# conda activate SAMap
# Imports
from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import numpy as np

# BLAST Results
map_dir = '/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/bin/samap_directory/mouse_mz/maps/'

# Protein -> Gene table
mm_prot = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/mouse_prot_table3_saunders.csv', sep=None)
mz_prot = pd.read_csv('/storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/mz_prot_table3.csv', sep=None)
mm_prot_tup = list(zip(mm_prot["protein_id"], mm_prot["gene"]))
mz_prot_tup = list(zip(mz_prot["protein_id"], mz_prot["gene"]))
protein_to_gene_names = {'mz': mz_prot_tup, 'mm':mm_prot_tup}

# Objects containing gene x counts
fn1 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/data/bb.h5ad'
fn2 = '/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/mouse_w_pc_down_norm.h5ad'
filenames = {'mz':fn1,'mm':fn2}
sam1=SAM()
sam1.load_data(fn1)
sam2=SAM()
sam2.load_data(fn2)
sams = {'mz':fn1,'mm':fn2}
sm = SAMAP(sams, names=protein_to_gene_names, f_maps = map_dir)
sm.run()
samap = sm.samap

# Save samap
import pickle
with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/bb_saunders_samap_3.pkl', 'wb') as out_file:
    pickle.dump(samap, out_file)


with open('/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/bb_saunders_sm_ran_3.pkl', 'wb') as out_file:
    pickle.dump(sm, out_file)

# Export as h5ad so it can be loaded into Seurat
samap.adata.write_h5ad("/storage/home/hcoda1/6/ggruenhagen3/scratch/bcs/data/bb_saunders_3.h5ad")
