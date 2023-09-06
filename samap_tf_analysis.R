# *** Transcription Factor Analysis *** #
# Load Libraries
library("plyr")
library("dplyr")
library("ggplot2")

# Load Gene Lists: TF, ligands, and receptors
s5_cat = read.csv("~/scratch/msc/S5_goi_by_cat_and_other.csv")
s5_cat = s5_cat[which(!duplicated(s5_cat$human)),]
tfs_turtle = data.table::fread("~/scratch/msc/tosches_tf_list.txt", data.table = F)
tfs_turtle$X=NULL
tfs_turtle = tfs_turtle[which(!tfs_turtle[,2] %in% s5_cat$human),]

# Turtle -> Mouse ==============================================================
# Load Genes
df = read.csv("~/scratch/bcs/results/turtle_zeisel_samap_genes_clean.csv")
is_sig_df = read.csv(paste0("~/scratch/bcs/results/turtle_zeisel_samap_genes_num_w_score.csv"))
df = df[which(df$id %in% is_sig_df$id[which(is_sig_df$sig)]),]
colnames(df) = c("X", "id", "cluster_1", "cluster_2", "gene_1", "gene_2", "hgnc_1", "hgnc_2")
# df = df[which(df$hgnc_1 == df$hgnc_2),]
df$hgnc = df$hgnc_2

# Find number and percentage of the gene lists that overlap the driver genes
df_sum = unique(df[,c("cluster_1", "cluster_2")])
df_sum$num_genes = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])), mc.cores = 20))
df_sum$num_tf = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$hgnc[which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])] %in% tfs_turtle[,2])), mc.cores = 20))
df_sum$num_lig = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$hgnc[which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])] %in% s5_cat$human[which(s5_cat$category == "ligand")])), mc.cores = 20))
df_sum$num_rec = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$hgnc[which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])] %in% s5_cat$human[which(s5_cat$category == "receptor")])), mc.cores = 20))
df_sum$num_other = df_sum$num_genes - df_sum$num_tf - df_sum$num_lig - df_sum$num_rec
df_sum$pct_tf  = (df_sum$num_tf  / df_sum$num_genes) * 100
df_sum$pct_lig = (df_sum$num_lig / df_sum$num_genes) * 100
df_sum$pct_rec = (df_sum$num_rec / df_sum$num_genes) * 100
df_sum$pct_other = (df_sum$num_other / df_sum$num_genes) * 100

# Group celltypes by brain region
cat1_ct = c("e07", "e08", "e14")
cat2_ct = c("e01", "e02")
neo_ct = c("TEGLU4", "TEGLU9", "TEGLU10")
df_sum$cat = "other"
df_sum$cat[which(df_sum$cluster_1 %in% cat1_ct & df_sum$cluster_2 %in% neo_ct)] = "aDC"
df_sum$cat[which(df_sum$cluster_1 %in% cat2_ct & df_sum$cluster_2 %in% neo_ct)] = "aDVR"
df_sum$cat2 = paste0(df_sum$cluster_1, ";", df_sum$cluster_2)
df_sum$cat2[which(df_sum$cat == "other")] = "other"

# Calculate mean by brain region
df_p_tf    = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_tf),  sd_num = sd(num_tf),  mean_pct = mean(pct_tf),  sd_pct = sd(pct_tf))  %>% mutate(gene_cat = "tf")
df_p_lig   = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_lig), sd_num = sd(num_lig), mean_pct = mean(pct_lig), sd_pct = sd(pct_lig)) %>% mutate(gene_cat = "lig")
df_p_rec   = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_rec), sd_num = sd(num_rec), mean_pct = mean(pct_rec), sd_pct = sd(pct_rec)) %>% mutate(gene_cat = "rec")
df_p_other = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_other), sd_num = sd(num_other), mean_pct = mean(pct_other), sd_pct = sd(pct_other)) %>% mutate(gene_cat = "other")
df_p = rbind(df_p_tf, df_p_lig, df_p_rec, df_p_other)

# Calculate stats
df_p$total_dif = df_p$mean_num_genes - df_p$mean_num
df_p[,c("p", "stat")] = 1
df_p$p[which(df_p$gene_cat == "tf")]    = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "tf"),    c("mean_num", "total_dif")]))$p.value
df_p$p[which(df_p$gene_cat == "lig")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "lig"),   c("mean_num", "total_dif")]))$p.value
df_p$p[which(df_p$gene_cat == "rec")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "rec"),   c("mean_num", "total_dif")]))$p.value
df_p$p[which(df_p$gene_cat == "other")] = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "other"), c("mean_num", "total_dif")]))$p.value
df_p$stat[which(df_p$gene_cat == "tf")]    = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "tf"),    c("mean_num", "total_dif")]))$statistic
df_p$stat[which(df_p$gene_cat == "lig")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "lig"),   c("mean_num", "total_dif")]))$statistic
df_p$stat[which(df_p$gene_cat == "rec")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "rec"),   c("mean_num", "total_dif")]))$statistic
df_p$stat[which(df_p$gene_cat == "other")] = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "other"), c("mean_num", "total_dif")]))$statistic
df_p$gene_cat = factor(df_p$gene_cat, levels = c("tf", "lig", "rec", "other"))
df_p_turtle = df_p
df_p_turtle$species = "t_to_m"

# Cichlid -> Mouse =============================================================
df = read.csv("~/scratch/bcs/results/bb_zeisel_samap_genes_clean.csv")
is_sig_df = read.csv(paste0("~/scratch/bcs/results/bb_zeisel_samap_genes_num_w_score.csv"))
df = df[which(df$id %in% is_sig_df$id[which(is_sig_df$sig)]),]
colnames(df) = c("X", "id", "cluster_1", "cluster_2", "gene_1", "gene_2", "hgnc_1", "hgnc_2")
df$hgnc = df$hgnc_2

# Find number and percentage of the gene lists that overlap the driver genes
df_sum = unique(df[,c("cluster_1", "cluster_2")])
df_sum$num_genes = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])), mc.cores = 20))
df_sum$num_tf = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$hgnc[which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])] %in% tfs_turtle[,2])), mc.cores = 20))
df_sum$num_lig = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$hgnc[which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])] %in% s5_cat$human[which(s5_cat$category == "ligand")])), mc.cores = 20))
df_sum$num_rec = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$hgnc[which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])] %in% s5_cat$human[which(s5_cat$category == "receptor")])), mc.cores = 20))
df_sum$num_other = df_sum$num_genes - df_sum$num_tf - df_sum$num_lig - df_sum$num_rec
df_sum$pct_tf  = (df_sum$num_tf  / df_sum$num_genes) * 100
df_sum$pct_lig = (df_sum$num_lig / df_sum$num_genes) * 100
df_sum$pct_rec = (df_sum$num_rec / df_sum$num_genes) * 100
df_sum$pct_other = (df_sum$num_other / df_sum$num_genes) * 100

# Group celltypes by brain region
cat1_ct = c("8.1_Glut", "12_Glut")
neo_ct = c("TEGLU6", "TEGLU9")
df_sum$cat = "other"
df_sum$cat[which(df_sum$cluster_1 %in% cat1_ct & df_sum$cluster_2 %in% neo_ct)] = "Dl-g"
df_sum$cat2 = paste0(df_sum$cluster_1, ";", df_sum$cluster_2)
df_sum$cat2[which(df_sum$cat == "other")] = "other"
df_p_tf    = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_tf),  sd_num = sd(num_tf),  mean_pct = mean(pct_tf),  sd_pct = sd(pct_tf))  %>% mutate(gene_cat = "tf")
df_p_lig   = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_lig), sd_num = sd(num_lig), mean_pct = mean(pct_lig), sd_pct = sd(pct_lig)) %>% mutate(gene_cat = "lig")
df_p_rec   = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_rec), sd_num = sd(num_rec), mean_pct = mean(pct_rec), sd_pct = sd(pct_rec)) %>% mutate(gene_cat = "rec")
df_p_other = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_other), sd_num = sd(num_other), mean_pct = mean(pct_other), sd_pct = sd(pct_other)) %>% mutate(gene_cat = "other")
df_p = rbind(df_p_tf, df_p_lig, df_p_rec, df_p_other)

# Calculate stats
df_p$total_dif = df_p$mean_num_genes - df_p$mean_num
df_p[,c("p", "stat")] = 1
df_p$p[which(df_p$gene_cat == "tf")]    = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "tf"),    c("mean_num", "total_dif")]))$p.value
df_p$p[which(df_p$gene_cat == "lig")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "lig"),   c("mean_num", "total_dif")]))$p.value
df_p$p[which(df_p$gene_cat == "rec")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "rec"),   c("mean_num", "total_dif")]))$p.value
df_p$p[which(df_p$gene_cat == "other")] = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "other"), c("mean_num", "total_dif")]))$p.value
df_p$stat[which(df_p$gene_cat == "tf")]    = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "tf"),    c("mean_num", "total_dif")]))$statistic
df_p$stat[which(df_p$gene_cat == "lig")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "lig"),   c("mean_num", "total_dif")]))$statistic
df_p$stat[which(df_p$gene_cat == "rec")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "rec"),   c("mean_num", "total_dif")]))$statistic
df_p$stat[which(df_p$gene_cat == "other")] = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "other"), c("mean_num", "total_dif")]))$statistic
df_p$gene_cat = factor(df_p$gene_cat, levels = c("tf", "lig", "rec", "other"))
df_p_cichlid = df_p
df_p_cichlid$species = "c_to_m"

# Cichlid -> Turtle ============================================================
df = read.csv("~/scratch/bcs/results/bb_turtle_samap_genes_clean.csv")
is_sig_df = read.csv(paste0("~/scratch/bcs/results/bb_turtle_samap_genes_num_w_score.csv"))
df = df[which(df$id %in% is_sig_df$id[which(is_sig_df$sig)]),]
colnames(df) = c("X", "id", "cluster_1", "cluster_2", "gene_1", "gene_2", "hgnc_1", "hgnc_2")
df$hgnc = df$hgnc_1

# Find number and percentage of the gene lists that overlap the driver genes
df_sum = unique(df[,c("cluster_1", "cluster_2")])
df_sum$num_genes = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])), mc.cores = 20))
df_sum$num_tf = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$hgnc[which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])] %in% tfs_turtle[,2])), mc.cores = 20))
df_sum$num_lig = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$hgnc[which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])] %in% s5_cat$human[which(s5_cat$category == "ligand")])), mc.cores = 20))
df_sum$num_rec = unlist(mclapply(1:nrow(df_sum), function(x) length(which(df$hgnc[which(df$cluster_1 == df_sum$cluster_1[x] & df$cluster_2 == df_sum$cluster_2[x])] %in% s5_cat$human[which(s5_cat$category == "receptor")])), mc.cores = 20))
df_sum$num_other = df_sum$num_genes - df_sum$num_tf - df_sum$num_lig - df_sum$num_rec
df_sum$pct_tf  = (df_sum$num_tf  / df_sum$num_genes) * 100
df_sum$pct_lig = (df_sum$num_lig / df_sum$num_genes) * 100
df_sum$pct_rec = (df_sum$num_rec / df_sum$num_genes) * 100
df_sum$pct_other = (df_sum$num_other / df_sum$num_genes) * 100

# Group celltypes by brain region
df_sum$cat = "other"
df_sum$cat[which(df_sum$cluster_1 %in% c("8.2_Glut") & df_sum$cluster_2 %in% c("e16"))] = "Dl-g DC"
df_sum$cat[which(df_sum$cluster_1 %in% c("8.8_Glut") & df_sum$cluster_2 %in% c("e16"))] = "Dl-g DC"
df_sum$cat[which(df_sum$cluster_1 %in% c("12_Glut") & df_sum$cluster_2 %in% c("e14"))] = "Dl-g DC"
df_sum$cat[which(df_sum$cluster_1 %in% c("12_Glut") & df_sum$cluster_2 %in% c("e01"))] = "Dl-g DVR"
df_sum$cat2 = paste0(df_sum$cluster_1, ";", df_sum$cluster_2)
df_sum$cat2[which(df_sum$cat == "other")] = "other"
df_p_tf    = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_tf),  sd_num = sd(num_tf),  mean_pct = mean(pct_tf),  sd_pct = sd(pct_tf))  %>% mutate(gene_cat = "tf")
df_p_lig   = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_lig), sd_num = sd(num_lig), mean_pct = mean(pct_lig), sd_pct = sd(pct_lig)) %>% mutate(gene_cat = "lig")
df_p_rec   = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_rec), sd_num = sd(num_rec), mean_pct = mean(pct_rec), sd_pct = sd(pct_rec)) %>% mutate(gene_cat = "rec")
df_p_other = df_sum %>% group_by(cat2) %>% dplyr::summarize(mean_num_genes = mean(num_genes), sd_num_genes = sd(num_genes), mean_num = mean(num_other), sd_num = sd(num_other), mean_pct = mean(pct_other), sd_pct = sd(pct_other)) %>% mutate(gene_cat = "other")
df_p = rbind(df_p_tf, df_p_lig, df_p_rec, df_p_other)

# Calculate stats
df_p$total_dif = df_p$mean_num_genes - df_p$mean_num
df_p[,c("p", "stat")] = 1
df_p$p[which(df_p$gene_cat == "tf")]    = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "tf"),    c("mean_num", "total_dif")]))$p.value
df_p$p[which(df_p$gene_cat == "lig")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "lig"),   c("mean_num", "total_dif")]))$p.value
df_p$p[which(df_p$gene_cat == "rec")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "rec"),   c("mean_num", "total_dif")]))$p.value
df_p$p[which(df_p$gene_cat == "other")] = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "other"), c("mean_num", "total_dif")]))$p.value
df_p$stat[which(df_p$gene_cat == "tf")]    = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "tf"),    c("mean_num", "total_dif")]))$statistic
df_p$stat[which(df_p$gene_cat == "lig")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "lig"),   c("mean_num", "total_dif")]))$statistic
df_p$stat[which(df_p$gene_cat == "rec")]   = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "rec"),   c("mean_num", "total_dif")]))$statistic
df_p$stat[which(df_p$gene_cat == "other")] = chisq.test(as.matrix(df_p[which(df_p$gene_cat == "other"), c("mean_num", "total_dif")]))$statistic
df_p$gene_cat = factor(df_p$gene_cat, levels = c("tf", "lig", "rec", "other"))
df_p_c_to_t = df_p
df_p_c_to_t$species = "c_to_t"

# Outputs ======================================================================
# Write stats to a csv file
df_p_all = rbind(df_p_turtle, df_p_cichlid, df_p_c_to_t)
write.csv(df_p_all, "~/scratch/bcs/results/bb_turtle_zeisel_composition_stats3.csv")

# Visualize the results
fname = "~/scratch/st/results/bb_turtle_zeisel_genecat_num3.pdf"
ggplot(df_p_all, aes(x = cat2, y = mean_num, fill = gene_cat)) + geom_bar(position="stack", stat = "identity") + theme_classic() + scale_y_continuous(expand= expansion(mult = c(0, .1)), name = "") + xlab("") + facet_grid(~ species, scales = "free_x", space = "free_x")
ggsave(fname, width = 11, height = 3)
paste0("rclone copy ", fname, " dropbox:BioSci-Streelman/George/Brain/spatial/analysis/samap/zeisel/")
