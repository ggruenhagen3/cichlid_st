# Input ========================================================================
# Read Input
# mouse.dataset = "zei_down"
# isGlut = F
# isGABA = F
# isNN   = F
# isSub1 = F
# isBB   = F
message(paste0("Running correlation comparison of the mouse Zeisel dataset vs the downsampled mouse Zeisel dataset"))

# Load Libraries
message("Loading Libraries")
suppressMessages(source("~/scratch/st/st_scripts/st_f.R"))
suppressMessages(source("~/scratch/bcs/bcs_scripts/bcs_f.R"))
library("ggh4x")

# Mouse Object
message("Loading mouse Zeisel object")
mouse = readRDS("~/scratch/bcs/data/l5_tel_all_sct.rds")

# Set Identity of Mouse Object
message("Setting mouse identity")
Idents(mouse) =  mouse$ClusterName

# Mouse DEG
message("Loading mouse DEGs")
mouse.deg  = read.csv("~/scratch/bcs/results/l5_cluster_markers_082823.csv")

# Downsampled mouse object
message("Loading downsampled mouse Zeisel object")
mz = readRDS("~/scratch/bcs/data/l5_tel_all_down_norm_022823.rds")

# Set Identity of Downsampled Mouse Object
message("Setting downsampled mouse identity")
Idents(mz) =  mz$ClusterName

# Downsampled mouse DEG
mz.deg = read.csv("~/scratch/bcs/results/zei_down_deg_082923.csv")

# Find Correlations ============================================================
# Find Common Gene Set
message("Finding correlations")
common.gene.set = sort(unique(toupper(mouse.deg$gene)))
common.gene.set = common.gene.set[which(common.gene.set %in% mz.deg$gene)]

# Cichlid Average Expression
mz.assay.to.use = "SCT"
mz.avg.exp = AverageExpression(mz, features = common.gene.set, assays = mz.assay.to.use, slot = "data")[[1]]
mz.avg.exp.norm = log(mz.avg.exp+1) + 0.1
mz.avg.exp.norm = mz.avg.exp.norm / rowMeans(mz.avg.exp.norm)

# Mouse Average Expression
mouse.avg.exp = AverageExpression(mouse, features = common.gene.set, assays = "SCT", slot = "data")[[1]]
mouse.avg.exp.norm = log(mouse.avg.exp+1) + 0.1
mouse.avg.exp.norm = mouse.avg.exp.norm / rowMeans(mouse.avg.exp.norm)

# Correlation
mz.mouse.cor = cor(mz.avg.exp.norm, mouse.avg.exp.norm, method = "spearman")

# Permutation Function =========================================================
permCor = function(old.mat) {
  new.mat.list = lapply(1:nrow(old.mat), function(x) sample(old.mat[x,]))
  new.mat = do.call('rbind', new.mat.list)
  perm.cor = cor(new.mat, mouse.avg.exp.norm, method = "spearman")
  perm.cor.melt = reshape2::melt(perm.cor)
  return(perm.cor.melt[,3])
}

# Permutations =================================================================
message("Performing permutations")
n.perms = 10000
mz.mouse.cor.melt = reshape2::melt(mz.mouse.cor)
orig.gene.labels = unname(rownames(mz.avg.exp.norm))
perm.cor.list = mclapply(1:n.perms, function(x) permCor(mz.avg.exp.norm), mc.cores = 20)
perm.cor.mat = do.call('cbind', perm.cor.list)
mz.mouse.cor.melt$num_perm_greater = unlist(lapply(1:nrow(mz.mouse.cor.melt), function(x) length(which(perm.cor.mat[x,] > mz.mouse.cor.melt$value[x]))))
mz.mouse.cor.melt$p.perm = (mz.mouse.cor.melt$num_perm_greater) / n.perms
mz.mouse.cor.melt$bon.perm = p.adjust(mz.mouse.cor.melt$p.perm, method = "bonferroni")

# Spearman Correlation Test ====================================================
message("Performing spearman correlation test")
mz.mouse.cor.p = matrix(NA, nrow = nrow(mz.mouse.cor), ncol = ncol(mz.mouse.cor), dimnames = list(rownames(mz.mouse.cor), colnames(mz.mouse.cor)))
for (mz.clust in colnames(mz.avg.exp.norm)) {
  for (mouse.clust in colnames(mouse.avg.exp.norm)) {
    mz.mouse.cor.p[mz.clust, mouse.clust] = cor.test(mz.avg.exp.norm[,mz.clust], mouse.avg.exp.norm[,mouse.clust], method = "spearman", alternative = "greater")$p.value
  }
}
mz.mouse.cor.bon = matrix(p.adjust(mz.mouse.cor.p, method = "bonferroni"), nrow = nrow(mz.mouse.cor), ncol = ncol(mz.mouse.cor), dimnames = list(rownames(mz.mouse.cor), colnames(mz.mouse.cor)))

# Output: Plots and CSV ========================================================
message("Plotting")
mz.mouse.cor.maxed.out.melt = reshape2::melt(mz.mouse.cor)
colnames(mz.mouse.cor.maxed.out.melt) = c("mz.cluster", "mouse.cluster", "cor")
mz.mouse.cor.maxed.out.melt[, c("mouse.region", "mouse.celltype")] = reshape2::colsplit(mz.mouse.cor.maxed.out.melt$mouse.cluster, "_", c('1', '2'))
mz.mouse.cor.maxed.out.melt = mz.mouse.cor.maxed.out.melt[,c("mz.cluster", "mouse.cluster", "mouse.region", "mouse.celltype", "cor")]
maxed.num=1
mz.mouse.cor.maxed.out.melt$cor.maxed = mz.mouse.cor.maxed.out.melt$cor
mz.mouse.cor.maxed.out.melt$cor.maxed[which(mz.mouse.cor.maxed.out.melt$cor >  maxed.num)] =  maxed.num
mz.mouse.cor.maxed.out.melt$cor.maxed[which(mz.mouse.cor.maxed.out.melt$cor < -maxed.num)] = -maxed.num
mz.mouse.cor.maxed.out.melt$cor.test.p = reshape2::melt(mz.mouse.cor.p)[,3]
mz.mouse.cor.maxed.out.melt$cor.test.bon = reshape2::melt(mz.mouse.cor.bon)[,3]
mz.mouse.cor.maxed.out.melt$perm.test.p = mz.mouse.cor.melt$p.perm
mz.mouse.cor.maxed.out.melt$perm.test.bon = mz.mouse.cor.melt$bon.perm
mz.mouse.cor.maxed.out.melt$all.sig = mz.mouse.cor.maxed.out.melt$perm.test.p < 0.05
mz.mouse.cor.maxed.out.melt_no_order = mz.mouse.cor.maxed.out.melt
mz.order  = hclust(dist(mz.mouse.cor), method = "complete")
mz.mouse.cor.maxed.out.melt$mz.cluster = factor(mz.mouse.cor.maxed.out.melt$mz.cluster, levels = mz.order$labels[mz.order$order])
mouse.order = hclust(dist(t(mz.mouse.cor)), method = "complete")
mz.mouse.cor.maxed.out.melt$mouse.cluster = factor(mz.mouse.cor.maxed.out.melt$mouse.cluster, levels = mouse.order$labels[mouse.order$order])
ggplot(mz.mouse.cor.maxed.out.melt, aes(x = mouse.cluster, y = mz.cluster, fill = cor.maxed)) + geom_raster() + geom_point(data = mz.mouse.cor.maxed.out.melt[which(mz.mouse.cor.maxed.out.melt$all.sig),], size = 1.2, color = "gray60") + scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), n.breaks = 6, limits = c(-maxed.num, maxed.num)) + coord_fixed() + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 10), axis.line=element_blank()) + force_panelsizes(rows = unit(nrow(mz.mouse.cor)/8, "in"), cols = unit(ncol(mz.mouse.cor)/8, "in"))


out_name = paste0("~/scratch/bcs/results/zei_down_zeisel_cluster_cor")
buffer.width = 2
ggsave(paste0(out_name, ".pdf"), width = (ncol(mz.mouse.cor)/5) + buffer.width, height = (nrow(mz.mouse.cor)/5) + buffer.width)
write.csv(mz.mouse.cor.maxed.out.melt, paste0(out_name, ".csv"))
message(paste0("rclone copy ", paste0(out_name, ".pdf"), " dropbox:BioSci-Streelman/George/Brain/spatial/analysis/bcsm/", mouse.dataset))
system( paste0("rclone copy ", paste0(out_name, ".pdf"), " dropbox:BioSci-Streelman/George/Brain/spatial/analysis/bcsm/", mouse.dataset))
system( paste0("rclone copy ", paste0(out_name, ".csv"), " dropbox:BioSci-Streelman/George/Brain/spatial/analysis/bcsm/", mouse.dataset))
message("Done")

mz.mouse.cor.maxed.out.melt_no_order$cor.maxed[which(mz.mouse.cor.maxed.out.melt_no_order$cor < 0)] = 0
ggplot(mz.mouse.cor.maxed.out.melt_no_order, aes(x = mouse.cluster, y = mz.cluster, fill = cor.maxed)) + geom_raster() + scale_fill_viridis(limits = c(0, 1)) + coord_fixed() + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(vjust = 0.5, hjust = 1, size = 10), axis.line=element_blank()) + force_panelsizes(rows = unit(nrow(mz.mouse.cor)/8, "in"), cols = unit(ncol(mz.mouse.cor)/8, "in"))
ggsave(paste0(out_name, "_viridis.pdf"), width = (ncol(mz.mouse.cor)/5) + buffer.width, height = (nrow(mz.mouse.cor)/5) + buffer.width)
