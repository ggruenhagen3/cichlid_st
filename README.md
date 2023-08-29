# cichlid_st: Scripts for Cichlid Spatial Transcriptomics Analysis
## Overview
The scripts in this github repository are for the analysis of spatial transcriptomics and single nucleus RNA sequencing data from the cichlid forebrain. These datasets can be found here [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217619) (spatial transcriptomics dataset is the GSE217615 SubSeries and the snRNA-seq dataset is the GSE211477 SubSeries). Additionally, this repository contains scripts for the comparison of these data to forebrain datasets from other vertebrate species. The following are links to the datasets used for each species: [axolotl](https://zenodo.org/record/6390083), [bird](https://cloud.biohpc.swmed.edu/index.php/s/nLicEtkmjGGmRF8), [turtle](https://public.brain.mpg.de/Laurent/ReptilePallium2018/), [mouse scRNA-seq (Zeisel et al)](http://mousebrain.org/adolescent/downloads.html), [mouse scRNA-seq (Saunders et al)](http://dropviz.org/), and [mouse spatial transcriptomics (Ortiz et al)](https://www.molecularatlas.org/download-data). Finally, these scripts contain hard-coded file paths and have not been designed nor tested on other systems.

This repository contains 3 directories and the code in each will be explained in detail below.
1. [clustering](https://github.com/ggruenhagen3/cichlid_st/tree/main/clustering): scripts for clustering spatial transcriptomics data
2. [cross_species](https://github.com/ggruenhagen3/cichlid_st/tree/main/cross_species): scripts for cross-species comparisons of transcriptomics data
3. [misc](https://github.com/ggruenhagen3/cichlid_st/tree/main/misc): miscellaneous scripts

## clustering
This folder contains code to find the near-optimal clustering parameters of our cichlid spatial transcriptomics data utilizing the ChooseR approach. First, `hb_params.bash` is used create a SLURM job script for every combination of the parameters. Thse SLURM jobs call `hb_chooser.R`, which will perform bootstraps of the data and determine the silhoutte width of the combination of parameters. Finally, `hb_chooser_collect.R` collects the outputs from `hb_chooser.r` and is used to evaluate which combination of parameters is near-optimal.

## cross_species
This folder contains scripts to compare forebrain transcriptomics data across multiple vertebrate species. The 

## misc
This folder contains many scripts for varied purposes, mainly visualization of results.
- `cell2location.py` uses cell2location to predict the location and abundance of celltypes from cichlid snRNA-seq data in spatial transcriptomics data
- `cell2location_plots.R` plots the results from cell2location, including the composition of each cichlid brain region by celltype
- `downsample_zeisel.R` downsamples the mouse scRNA-seq dataset from Zeisel et al by randomly removing 10% of UMIs, allowing for a comparison to the original dataset
- `plot_cichlid_zeisel_umap.R` plots cichlid nuclei and mouse cells together in UMAP space after integration performed by SAMap
- `samap_celltype_mapping_plot.R` plots the similarity score and significance of celltypes between species after integration by SAMap
