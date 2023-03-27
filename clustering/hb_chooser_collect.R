#*******************************************************************************
# Clustering Robustness ========================================================
#*******************************************************************************
# All hard-coded paths in this block
# git_dir_path: path to git directory
# obj_path: unclustered Seurat Object
# results_path: path to store results in 
git_dir_path = "/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/chooser/my/cluster.stability/"
obj_path = "/storage/home/hcoda1/6/ggruenhagen3/scratch/st/data/st_c2b2_010622.rds"
results_path = paste0("/storage/home/hcoda1/6/ggruenhagen3/scratch/brain/chooser/chooser_spatial/examples/")
output_path = results_path

source(paste0(git_dir_path, "clustering/chooseR_code/pipeline.R"))
batch_num = 1
results_path = paste0(results_path, "batch_data_", batch_num, "/")
min_dists = c(0.001, 0.1, 0.2, 0.3, 0.4, 0.5)
n_neighbors = c(5, 10, 20, 30, 40, 50)
resolutions = seq(0.2, 2.0, by = 0.2)
scores = list()
for (min_dist in min_dists) {
  for (n_neighbor in n_neighbors) {
    for (res in resolutions) {
      this.str = paste0("X", batch_num, "_", min_dist, "_", n_neighbor, "_res.", res)
      this.sil =  readRDS(paste0(results_path, "silhouette_grouped_", this.str, ".rds"))
      this.sil$min_dist = min_dist
      this.sil$n_neighbor = n_neighbor
      scores[[length(scores)+1]] = this.sil
    }
  }
}

library("dplyr")
scores <- dplyr::bind_rows(scores) %>%
  dplyr::group_by(min_dist, n_neighbor, res) %>%
  dplyr::mutate("n_clusters" = dplyr::n()) %>%
  dplyr::ungroup()

unique(data.frame(scores)[,c("min_dist", "n_neighbor", "res", "n_clusters")])
range(unique(data.frame(scores)[,c("min_dist", "n_neighbor", "res", "n_clusters")])$n_clusters)

meds <- scores %>%
  dplyr::group_by(min_dist, n_neighbor, res) %>%
  dplyr::summarise(
    "boot" = list(boot_median(avg_sil)),
    "n_clusters" = mean(n_clusters)
  ) %>%
  tidyr::unnest_wider(boot)

threshold <- max(meds$low_med)
meds %>% dplyr::filter(med >= threshold) %>% dplyr::arrange(n_clusters) %>% tail(n = 1)
meds$id = paste(meds$min_dist, meds$n_neighbor, meds$res, sep = "_")
scores$id = paste(scores$min_dist, scores$n_neighbor, scores$res, sep = "_")
meds = meds %>% dplyr::arrange(n_clusters)
choice <- as.character(meds %>% dplyr::filter(med >= threshold) %>% dplyr::arrange(n_clusters) %>% tail(n = 1) %>% dplyr::pull(id))
print("The best clustering paramaters are:")
print(choice)

p.meds = data.frame(meds)
p.meds$id = factor(p.meds$id, levels = unique(p.meds$id))
ggplot(p.meds, aes(id, med)) + geom_crossbar(aes(ymin = low_med, ymax = high_med), fill = "grey", size = 0.25) +
                               geom_hline(aes(yintercept = threshold), colour = "blue") +
                               geom_vline(aes(xintercept = choice), colour = "red") +
                               scale_x_discrete("Min Distance _ N Neighbors _ Resolution") +
                               scale_y_continuous("Silhouette Score", expand = c(0, 0), limits = c(-1, 1), breaks = seq(-1, 1, 0.25), oob = scales::squish) +
                               cowplot::theme_minimal_hgrid() +
                               theme(axis.title = element_text(size = 8), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 45, vjust = 1, hjust = 1), axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), axis.ticks = element_line(colour = "black"),)
ggsave(filename = paste0(output_path, "chooser_near_optimal.png"), dpi = 300, height = 3.5, width = 7, units = "in")

