# Read Input
args = commandArgs(trailingOnly=TRUE)
batch_num   = as.numeric(args[1])
min_dist    = as.numeric(args[2])
n_neighbors = as.numeric(args[3])
num_cores   = 5
if (length(args) > 3) { num_cores = as.numeric(args[4]) }
print(paste0("Using parameters: batch_num=", batch_num, ", min_dist=", min_dist, ", n_neighbors=", n_neighbors, ", num_cores=", num_cores))

# All hard-coded paths in this block
# git_dir_path: path to git directory
# obj_path: unclustered Seurat Object
# results_path: path to store results in 
git_dir_path = "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/chooser/my/cluster.stability/"
obj_path = "/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/st_c2b2_010622.rds"
results_path = paste0("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/chooser/chooser_spatial/examples/")

# Load Libraries
source(paste0(git_dir_path, "clustering/chooseR_code/pipeline.R"))
library(Seurat)
library(ggplot2)
library(parallel)
`%>%` <- magrittr::`%>%`

# Load Object
obj = readRDS(obj_path)

# Define the number of PCs to use, and which assay and reduction to use.
npcs <- 2
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
assay <- "SCT"
reduction <- "umap"
results_path = paste0(results_path, "batch_data_", batch_num, "/")

# Run pipeline
pipeline_part1 = function(res) {
  message(paste0("Clustering ", res, "..."))
  message("\tFinding ground truth...")

  # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
  this.obj <- find_clusters(
                       obj,
                       reduction = reduction,
                       assay = assay,
                       resolution = res,
                       npcs = npcs,
                       batch.num = batch_num,
                       min.dist = min_dist,
                       n.neighbors = n_neighbors
  )
  this.str = paste0("X", batch_num, "_", min_dist, "_", n_neighbors, "_res.", res)
  clusters <- this.obj[[glue::glue(this.str)]]

  # Now perform iterative, sub-sampled clusters
  results <- multiple_cluster(
                              this.obj,
                              n = 100,
                              size = 0.8,
                              npcs = npcs,
                              res = res,
                              reduction = reduction,
                              assay = assay,
                              min.dist = min_dist,
                              n.neighbors = n_neighbors
  )
  rm(this.obj) # memory cleanup
  saveRDS(results, paste0(results_path, "results_", this.str, ".rds"))
  saveRDS(clusters,paste0(results_path, "clusters_", this.str, ".rds"))
}

pipeline_part2 = function(res) {
  this.str = paste0("X", batch_num, "_", min_dist, "_", n_neighbors, "_res.", res)
  results = readRDS(paste0(results_path, "results_", this.str, ".rds"))

  # Now calculate the co-clustering frequencies
  message(paste0("Tallying ", res, "..."))
  columns <- colnames(dplyr::select(results, -cell))
  mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
  i <- 1 # Counter
  for (col in columns) {
    message(paste0("\tRound ", i, "..."))
    mtchs <- Reduce("+", list(
      mtchs,
      find_matches(col, df = results)
    ))
    i <- i + 1
  }

  message(paste0("Scoring ", res, "..."))
  mtchs <- dplyr::mutate_all(
    dplyr::as_tibble(mtchs),
    function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
  )
  saveRDS(mtchs, paste0(results_path, "mtchs_", this.str, ".rds"))
}

pipeline_part3 = function(res) {
  this.str = paste0("X", batch_num, "_", min_dist, "_", n_neighbors, "_res.", res)
  mtchs    = readRDS(paste0(results_path, "mtchs_", this.str, ".rds"))
  clusters = readRDS(paste0(results_path, "clusters_", this.str, ".rds"))
  
  # Now calculate silhouette scores
  message(paste0("Silhouette ", res, "..."))
  sil <- cluster::silhouette(
    x = as.numeric(as.character(unlist(clusters))),
    dmatrix = (1 - as.matrix(mtchs))
  )
  saveRDS(sil, paste0(results_path, "silhouette_", this.str, ".rds"))

  # Finally, calculate grouped metrics
  message(paste0("Grouping ", res, "..."))
  grp <- group_scores(mtchs, unlist(clusters))
  saveRDS(grp, paste0(results_path, "frequency_grouped_", this.str, ".rds"))
  sil <- group_sil(sil, res)
  saveRDS(sil, paste0(results_path, "silhouette_grouped_", this.str, ".rds"))
  write.csv(clusters, paste0(results_path, "clusters_", this.str, ".csv"))
}

message(" - - - Clustering (George) - - - ")
mclapply(resolutions, function(x) pipeline_part1(x), mc.cores = num_cores)
message(" - - - Done Clustering (George) - - - ")

message(" - - - Tallying/Scoring (George) - - - ")
mclapply(resolutions, function(x) pipeline_part2(x), mc.cores = num_cores)
message(" - - - Done Tallying/Scoring (George) - - - ")

message(" - - - Silhoutte (George) - - - ")
mclapply(resolutions, function(x) pipeline_part3(x), mc.cores = num_cores)
message(" - - - Done Silhoutte (George) - - - ")

