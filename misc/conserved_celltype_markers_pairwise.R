# Helper Functions =============================================================
geneOverlap = function(df1, df2, cluster1, cluster2, hgnc1, hgnc2, return_genes = F, return_ovlp = F) {
  df1 = df1[which(df1$cluster == cluster1),]
  df2 = df2[which(df2$cluster == cluster2),]
  
  df1$hgnc = hgnc1$hgnc[match(df1$gene, hgnc1$gene)]
  df2$hgnc = hgnc2$hgnc[match(df2$gene, hgnc2$gene)]
  
  df1_genes = df1$hgnc[which(!is.na(df1$hgnc))]
  df2_genes = df2$hgnc[which(!is.na(df2$hgnc))]
  
  ovlp = unique(df1_genes[which(df1_genes %in% df2_genes)])
  smallest_cluster = min(length(df1_genes), length(df2_genes))
  pct = length(ovlp) / smallest_cluster
  pct = ifelse(length(ovlp) == 0, 0, pct)
  
  if (return_ovlp) {
    return(ifelse(length(ovlp) == 0, 0, ovlp))
  }
  
  if (return_genes) {
    df_ortho = merge(df1[which(!is.na(df1$hgnc)),], df2[which(!is.na(df1$hgnc)),], by = "hgnc", suffixes = c("_1", "_2"))
    df_ortho = df_ortho[order(df_ortho$p_val_adj_1, df_ortho$p_val_adj_2),]
    df_ortho = df_ortho[which(!duplicated(df_ortho$gene_2)),]
    df_ortho = df_ortho[order(df_ortho$p_val_adj_2, df_ortho$p_val_adj_1),]
    df_ortho = df_ortho[which(!duplicated(df_ortho$gene_1)),]
    df_ortho = df_ortho[,which(!colnames(df_ortho) %in% c("X_1", "X_2"))]
    return(df_ortho)
  }
  return(pct)
}

# Gene Orthologs ===============================================================
vert = sort(c("mouse", "axolotl", "bird", "turtle", "cichlid"))
ortho = list()
for (v1 in vert[1:length(vert)]) {
  if (v1 == "mouse") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/mouse_human_many.csv"))
    mz_mm_gene_map = mz_mm_gene_map[,c("gene", "hgnc")]
  } else if (v1 == "turtle") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/turtle_ortho_final.csv"))
    mz_mm_gene_map$hgnc = mz_mm_gene_map$hgnc.many
  } else if (v1 == "cichlid") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/gene_info_4.csv"))
    mz_mm_gene_map$hgnc = mz_mm_gene_map$human_reasonable
    mz_mm_gene_map$gene = mz_mm_gene_map$seurat_name
  } else if (v1 == "bird") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/bird_ortho_final.csv"))
    mz_mm_gene_map$hgnc = mz_mm_gene_map$hgnc.many
  } else if (v1 == "axolotl") {
    mz_mm_gene_map = read.csv(paste0("~/scratch/bcs/data/", v1, "_egg_many.csv"))  
  }
  mz_mm_gene_map = mz_mm_gene_map[,c("gene", "hgnc")]
  mz_mm_gene_map$hgnc[which(mz_mm_gene_map$hgnc == "")] = NA
  ortho[[v1]] = mz_mm_gene_map
}
ortho[["zeisel"]] = ortho[["mouse"]]
ortho[["saunders"]] = ortho[["mouse"]]
ortho[["oritz"]] = ortho[["mouse"]]

# Main =========================================================================
# Inputs
mz.dataset = "bb"
mm.dataset = "zeisel"
mouse.dataset = mm.dataset

# Load Datsets and DEGs
mz.hgnc = ortho[["cichlid"]]
mm.hgnc = ortho[[mm.dataset]]
mz_deg = read.csv("~/scratch/brain/results/bb53_deg_012323.csv")
if (mz.dataset == "st") {
  mz_deg = read.csv("~/scratch/bcs/results/c2b2_brain_structure_deg_042023.csv")
}
mz_deg$cluster = stringr::str_replace(as.character(as.vector(mz_deg$cluster)), "Astro", "RG")
mz_deg$cluster = stringr::str_replace(as.character(as.vector(mz_deg$cluster)), "OB gc", "OB-gc")
mz_deg$cluster = stringr::str_replace(as.character(as.vector(mz_deg$cluster)), "OB gml", "OB-gml")

message("Loading mouse DEGs")
mouse.deg.path = list()
mouse.deg.path[["axolotl"]]    = read.csv("~/scratch/bcs/results/axolotl_cluster_markers_022123.csv")
mouse.deg.path[["bird"]]    = read.csv("~/scratch/bcs/results/bird_markers_041923.csv")
mouse.deg.path[["oritz"]]    = read.csv("~/scratch/bcs/results/oritzg_markers_042023.csv")
mouse.deg.path[["saunders"]] = read.csv("~/scratch/bcs/results/saunders_cluster_region_subcluster_020723.csv")
mouse.deg.path[["turtle"]]   = read.csv("~/scratch/bcs/results/turtle_cluster_markers_041923.csv")
mouse.deg.path[["zeisel"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_082823.csv")
mouse.deg = mouse.deg.path[[mouse.dataset]]
mm_deg = mouse.deg

# Identify the percentage of overlapping genes
mm_clusters = unique(mm_deg$cluster)
mz_clusters = unique(mz_deg$cluster)
ovlp_df = expand.grid(mm_clusters, mz_clusters)
colnames(ovlp_df) = c("mm_cluster", "mz_cluster")
ovlp_df$mz_cluster = stringr::str_replace(as.character(as.vector(ovlp_df$mz_cluster)), "Astro", "RG")
ovlp_df$id = paste0(ovlp_df$mz_cluster, "_", ovlp_df$mm_cluster)
ovlp_df$pct_ovlp  = unlist(mclapply(1:nrow(ovlp_df), function(x) geneOverlap(mz_deg, mm_deg, as.vector(ovlp_df$mz_cluster[x]), as.vector(ovlp_df$mm_cluster[x]), mz.hgnc, mm.hgnc, return_genes = F), mc.cores = 20))

# Categorize by whether the celltype pairs are significant
hit_sup = read.csv(paste0("~/scratch/bcs/results/", mz.dataset, "_", mm.dataset, "_sig_hits.csv"))
hit_sup$X = NULL
hit_sup$id = paste0(hit_sup$mz_name, "_", hit_sup$mm_name)
ovlp_df[,"p0"] = hit_sup[match(ovlp_df$id, hit_sup$id), "p0"]
ovlp_df[which(is.na(ovlp_df[,"p0"])),"p0"] = FALSE

# Save the genes for celltype pairs that are significant
ovlp_df_p0 = ovlp_df[which(ovlp_df$p0),]
res = mclapply(1:nrow(ovlp_df_p0), function(x) geneOverlap(mz_deg, mm_deg, as.vector(ovlp_df_p0$mz_cluster[x]), as.vector(ovlp_df_p0$mm_cluster[x]), mz.hgnc, mm.hgnc, return_genes = T), mc.cores = 20)
df_ortho = do.call('rbind', res)
write.csv(df_ortho, paste0("~/scratch/bcs/results/", mz.dataset, "_", mm.dataset, "_hit_genes.csv"))

# Define the color palette
if (grepl("tasic", mm.dataset)) {
  mm_col = "mm_cluster"
  col.pal = rev(brewer.pal(11, "PRGn")[1:6])
} else if ( mm.dataset == "saunders" ) {
  mm_col = "mm_tissue_cluster"
  col.pal = rev(brewer.pal(11, "PiYG")[1:6])
} else if ( mm.dataset == "saunders_sub" ) {
  mm_col = "mm_tissue_subcluster"
  col.pal = rev(brewer.pal(11, "PiYG")[1:6])
} else if (grepl("oritz", mm.dataset)) {
  mm_col = "mm_g_parent"
  col.pal = brewer.pal(9, "Blues")
} else if (grepl("zeisel", mm.dataset)) {
  mm_col = "mm_ClusterName"
  col.pal = brewer.pal(9, "Greens")
} else if (grepl("turtle", mm.dataset)) {
  mm_col = "cp_detail"
  col.pal = brewer.pal(11, "BrBG")[6:11]
} else if (grepl("axolotl", mm.dataset)) {
  mm_col = "am_cluster"
  col.pal = rev(c(magma(100)[40:100], magma(100)[80:100]))
} else if (grepl("bird", mm.dataset)) {
  mm_col = "tg_cluster_orig2"
  col.pal = rev(brewer.pal(11, "PuOr")[1:6])
}

# Plot
ovlp_df$pct_ovlp = ovlp_df$pct_ovlp * 100
ggplot(ovlp_df, aes(x = pct_ovlp, fill = p0, color = p0)) + geom_density(alpha = 0.75) + scale_color_manual(values = c(colorRampPalette(colors = col.pal)(100)[50], colorRampPalette(colors = col.pal)(100)[100])) + scale_fill_manual(values = c(colorRampPalette(colors = col.pal)(100)[30], colorRampPalette(colors = col.pal)(100)[80])) + theme_classic() + scale_x_continuous(expand = c(0,0), name = "") + scale_y_continuous(expand = expansion(mult = c(0, .05)), name = "") + theme(axis.text = element_text(size = 10))
ggsave(paste0("~/scratch/bcs/results/", mz.dataset, "_", mm.dataset, "_pct_ovlp", ".pdf"), width = 5, height = 4)
system(paste0("rclone copy ~/scratch/bcs/results/", mz.dataset, "_", mm.dataset, "_pct_ovlp", ".pdf dropbox:BioSci-Streelman/George/Brain/spatial/analysis/samap/"))

# Calculate the mean and standard deviation
se <- function(x) sd(x)/sqrt(length(x))
this_p0 = ovlp_df$pct_ovlp[which(ovlp_df$p0)]
this_not_p0 = ovlp_df$pct_ovlp[which(!ovlp_df$p0)]
t.test(this_p0, this_not_p0)$p.value
t.test(this_p0, this_not_p0)$statistic

# Mean number of conserved marker genes
mean(data.frame(df_ortho %>% group_by(id) %>% dplyr::summarize(mean = n()))[,2])
se(data.frame(df_ortho %>% group_by(id) %>% dplyr::summarize(mean = n()))[,2])
