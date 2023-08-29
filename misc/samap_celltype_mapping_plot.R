# Load Libraries
library("stringr")
library("ggplot2")
library("Seurat")
library("SeuratObject")
library("Matrix")
library("RColorBrewer")
library("viridis")
library("reshape2")
library("data.table")
library("colourvalues")
library("colorspace")
library("patchwork")
library("dplyr")
library("parallel")
library("grid")
library("ggh4x")

# Inputs
mz.dataset = "bb"
mm.dataset = "zeisel2"
isGE = F  # plot just the ganglionic eminence celltypes?

# Load celltype mapping after SAMap
mzmm   = as.matrix(read.csv(paste0("~/research/st/results/samap_res/", mz.dataset, "_", mm.dataset, "_mapping_mine4.csv"),   row.names = 1))
mzmm.p = as.matrix(read.csv(paste0("~/research/st/results/samap_res/", mz.dataset, "_", mm.dataset, "_mapping_mine_p4.csv"), row.names = 1))

# Correct special character issues, if necessary
if (mz.dataset == "st") {
  colnames(mzmm) = str_replace_all(colnames(mzmm), "\\.", "-")
  colnames(mzmm) = plyr::revalue(colnames(mzmm), c("Dc-1-2" = "Dc-1/2"))
  
  colnames(mzmm.p) = str_replace_all(colnames(mzmm.p), "\\.", "-")
  colnames(mzmm.p) = plyr::revalue(colnames(mzmm.p), c("Dc-1-2" = "Dc-1/2"))
} else {
  colnames(mzmm)[which(startsWith(colnames(mzmm), "X"))] = str_sub(colnames(mzmm)[which(startsWith(colnames(mzmm), "X"))], 2, 50)
  colnames(mzmm) = str_replace(colnames(mzmm), "Astro", "RG")
  colnames(mzmm)[which(startsWith(colnames(mzmm), "mz_"))] = str_sub(colnames(mzmm)[which(startsWith(colnames(mzmm), "mz_"))], 4, 50)
  rownames(mzmm)[which(startsWith(rownames(mzmm), "mm_"))] = str_sub(rownames(mzmm)[which(startsWith(rownames(mzmm), "mm_"))], 4, 50)
  colnames(mzmm) = plyr::revalue(colnames(mzmm), c("8.9_Glut"="8-9_Glut", "8.9_Glut.1"="8.9_Glut", "8.9_Glut.1"="8.9_Glut", "15.1_GABA.Glut"="15.1_GABA/Glut", "15.5_GABA.Glut"="15.5_GABA/Glut"))

  colnames(mzmm.p)[which(startsWith(colnames(mzmm.p), "X"))] = str_sub(colnames(mzmm.p)[which(startsWith(colnames(mzmm.p), "X"))], 2, 50)
  colnames(mzmm.p) = str_replace(colnames(mzmm.p), "Astro", "RG")
  colnames(mzmm.p)[which(startsWith(colnames(mzmm.p), "mz_"))] = str_sub(colnames(mzmm.p)[which(startsWith(colnames(mzmm.p), "mz_"))], 4, 50)
  rownames(mzmm.p)[which(startsWith(rownames(mzmm.p), "mm_"))] = str_sub(rownames(mzmm.p)[which(startsWith(rownames(mzmm.p), "mm_"))], 4, 50)
  colnames(mzmm.p) = plyr::revalue(colnames(mzmm.p), c("8.9_Glut"="8-9_Glut", "8.9_Glut.1"="8.9_Glut", "8.9_Glut.1"="8.9_Glut", "15.1_GABA.Glut"="15.1_GABA/Glut", "15.5_GABA.Glut"="15.5_GABA/Glut"))
}

# Apply a cutoff so values with high scores don't skew the color legend
mzmm[which(mzmm > quantile(mzmm, 0.992))] = quantile(mzmm, 0.992)

# Reshape the matrix of scores to a dataframe
mzmm.melt = reshape2::melt(mzmm)
mzmm.melt = mzmm.melt[which(!is.na(mzmm.melt$value) & mzmm.melt$Var2 != ""),]
colnames(mzmm.melt) = c("mm_name", "mz_name", "Score")
mzmm.melt$id = paste0(mzmm.melt$mm_name, "_", mzmm.melt$mz_name)

# Reshape the matrix of p values to a dataframe
mzmm.p.melt = reshape2::melt(mzmm.p)
colnames(mzmm.p.melt) = c("mm_name", "mz_name", "p")
mzmm.p.melt$id = paste0(mzmm.p.melt$mm_name, "_", mzmm.p.melt$mz_name)

# Add p-value's from the p-value dataframe to the dataframe with the scores
mzmm.melt$p_perm = mzmm.p.melt$p[match(mzmm.melt$id, mzmm.p.melt$id)]
mzmm.melt$bh_perm = p.adjust(mzmm.melt$p_perm, method = "BH")
mzmm.melt$bh_sig = mzmm.melt$bh_perm < 0.05
mzmm.melt$p_sig  = mzmm.melt$p_perm < 0.05
mzmm.melt$p0     = mzmm.melt$p_perm == 0
mzmm.melt$bon_perm = p.adjust(mzmm.melt$p_perm, method = "bonferroni")
mzmm.melt$p0 = mzmm.melt$bon_perm < 0.05

# Color Pallette for plot
col.pal = viridis(100)
mouse.order = ""
if (grepl("tasic", mm.dataset)) {
  col.pal = rev(brewer.pal(11, "PRGn")[1:6])
} else if (grepl("saunders", mm.dataset)) {
  col.pal = rev(brewer.pal(11, "PRGn")[1:6])
} else if (grepl("oritz", mm.dataset)) {
  col.pal = brewer.pal(9, "Reds")
} else if (mm.dataset == "zeisel") {
  col.pal = brewer.pal(9, "Greens")
  mouse.order = read.csv("~/research/st/results/zeisel_cell-type_order_070323.csv")[,1]
} else if (mm.dataset == "zeisel2") {
  col.pal = brewer.pal(9, "Greens")
  mouse.order = read.csv("~/research/st/results/zeisel_cell-type_order_new_081523.csv")[,1]
} else if (grepl("turtle", mm.dataset)) {
  col.pal = brewer.pal(11, "BrBG")[6:11]
  mouse.order = read.csv("~/research/st/results/turtle_cell-type_order_070423.csv")[,1]
} else if (grepl("axolotl", mm.dataset)) {
  col.pal = rev(brewer.pal(11, "PiYG")[1:6])
  mouse.order = read.csv("~/research/st/results/axolotl_cell-type_order_070423.csv")[,1]
} else if (grepl("bird", mm.dataset)) {
  col.pal = rev(brewer.pal(11, "PuOr")[1:6])
  mouse.order = read.csv("~/research/st/results/bird_cell-type_order_070423.csv")[,1]
}

# Order the axis labels
if ( isGE ) {
  mzmm.melt$mm_name = as.vector(mzmm.melt$mm_name)
  mz_lge = c("4.1_GABA", "4.2_GABA", "4.3_GABA", "4.4_GABA", "4.5_GABA", "4.6_GABA", "4.7_GABA", "4.8_GABA", "7_GABA"); mz_cge = c("15.1_GABA/Glut", "15.2_GABA", "15.4_GABA"); mz_mge = c("15.3_GABA", "15.5_GABA", "6_GABA");
  if (grepl("tasic", mm.dataset)) {
    mm_lge = c(); mm_cge = unique(mzmm.melt$mm_name[which(startsWith(as.vector(mzmm.melt$mm_name), "Vip"))]); mm_mge = unique(mzmm.melt$mm_name[which(startsWith(as.vector(mzmm.melt$mm_name), "Sst") | startsWith(as.vector(mzmm.melt$mm_name), "Pvalb"))]);    
  } else if (grepl("zeisel", mm.dataset)) {
    mm_lge = c("MSN1", "MSN2", "MSN3", "MSN4", "MSN5"); mm_cge = c("TEINH1", "TEINH4", "TEINH12"); mm_mge = c("TEINH16", "TEINH17", "TEINH18", "TEINH19", "TEINH21");  
  } else if (grepl("turtle", mm.dataset)) {
    mm_lge = c("i01", "i04", "i05", "i06"); mm_cge = c("i14", "i15", "i16", "i17", "i18"); mm_mge = c("i07", "i08", "i09", "i10", "i11", "i12", "i13");
  } else if (grepl("axolotl", mm.dataset)) {
    mm_lge = paste0("GABA", c("1", "3", "9", "10", "11", "12", "22", "25", "29")); mm_cge = c(); mm_mge = paste0("GABA", c("2", "6", "17", "20"));
  } else if (grepl("bird", mm.dataset)) {
    mm_lge = c("MSN1", "MSN2", "MSN3", "MSN4"); mm_cge = c("GABA-5-1", "GABA-5-2", "GABA-5-3"); mm_mge = c("GABA-2", "GABA-3", "GABA-4", "GABA-6");
  }
  
  mzmm.melt = mzmm.melt[which(mzmm.melt$mz_name %in% c(mz_lge, mz_cge, mz_mge)),]
  mzmm.melt = mzmm.melt[which(mzmm.melt$mm_name %in% c(mm_lge, mm_cge, mm_mge)),]
  mzmm.melt$mz_name = factor(as.vector(mzmm.melt$mz_name), levels = c(mz_lge, mz_cge, mz_mge))
  mzmm.melt$mm_name = factor(as.vector(mzmm.melt$mm_name), levels = c(mm_lge, mm_cge, mm_mge))
  mzmm = reshape2::acast(mzmm.melt, mm_name ~ mz_name, value.var = "bh_perm")
} else {
  if (mz.dataset == "st") {
    mz.order = c("OB-gc", "OB-gml", "Vs", "Vd-c", "Vd-r", "Vx", "Vi", "Vl", "Vv", "Vc", "Dm-1", "Dl-d", "Dl-vv", "Dl-v", "Dl-g", "Dp", "Dm-2c", "Dm-2r", "Dc-3", "Dd", "vVZ", "VZ", "ON", "Dm-3", "Dc-1/2", "Dc-4", "Dc-5")  
    mm.order = c("MOB", "OLF", "ACB", "PAL", "LSX", "CEA", "sAMY", "NLOT", "MEA", "sAMY-other", "CP", "STRv-other", "CA3", "DG", "CA1", "CA2", "HIP-other", "VIS", "ORB", "PTLp", "AUD", "ILA", "ACA", "MO", "PERI", "PL", "TEa", "AON", "VS", "fiber tracts", "DP", "RHP", "ECT", "AI", "RSP", "BMA", "PA", "EP", "LA", "FS", "BLA", "U_CTX", "CLA", "GU", "SS", "VISC", "PIR", "TT", "PAA", "TR")
    mzmm.melt$mz_name = factor(as.vector(mzmm.melt$mz_name), levels = mz.order)
    mzmm.melt$mm_name = factor(mzmm.melt$mm_name, levels = mm.order)
  } else if (mouse.order != "") {
    mz.order = c("1.1_RG", "1.2_RG", "2.2_Oligo", "2.1_OPC", "1.3_MG", "3_Peri", "9.5_Glut", "15.5_GABA/Glut", "15.7_Glut", "14_Glut", "15.6_Glut", "5.1_GABA", "5.2_GABA", "6_GABA", "15.3_GABA", "15.4_GABA", "15.1_GABA/Glut", "15.2_GABA", "4.8_GABA", "4.4_GABA", "4.3_GABA", "4.6_GABA", "4.5_GABA", "4.2_GABA", "4.1_GABA", "4.7_GABA", "7_GABA", "8.1_Glut", "8.2_Glut", "8.8_Glut", "8.11_Glut", "12_Glut", "8.10_Glut", "8.7_Glut", "11.3_Glut", "11.2_Glut", "10.2_Glut", "10.1_Glut", "11.1_Glut", "8.4_Glut", "8.3_Glut", "8.9_Glut", "8.5_Glut", "8.6_Glut", "8-9_Glut", "9.8_Glut", "13_Glut", "9.6_Glut", "9.7_Glut", "9.4_Glut", "9.1_Glut", "9.2_Glut", "9.3_Glut")
    mzmm.melt$mz_name = factor(as.vector(mzmm.melt$mz_name), levels = mz.order)
    mzmm.melt$mm_name = factor(mzmm.melt$mm_name, levels = mouse.order)
  } else {
    mz.order = c("1.1_RG", "1.2_RG", "2.2_Oligo", "2.1_OPC", "1.3_MG", "3_Peri", "9.5_Glut", "15.5_GABA/Glut", "15.7_Glut", "14_Glut", "15.6_Glut", "5.1_GABA", "5.2_GABA", "6_GABA", "15.3_GABA", "15.4_GABA", "15.1_GABA/Glut", "15.2_GABA", "4.8_GABA", "4.4_GABA", "4.3_GABA", "4.6_GABA", "4.5_GABA", "4.2_GABA", "4.1_GABA", "4.7_GABA", "7_GABA", "8.1_Glut", "8.2_Glut", "8.8_Glut", "8.11_Glut", "12_Glut", "8.10_Glut", "8.7_Glut", "11.3_Glut", "11.2_Glut", "10.2_Glut", "10.1_Glut", "11.1_Glut", "8.4_Glut", "8.3_Glut", "8.9_Glut", "8.5_Glut", "8.6_Glut", "8-9_Glut", "9.8_Glut", "13_Glut", "9.6_Glut", "9.7_Glut", "9.4_Glut", "9.1_Glut", "9.2_Glut", "9.3_Glut")
    mzmm.melt$mz_name = factor(as.vector(mzmm.melt$mz_name), levels = mz.order)
    mouse.order = hclust(dist(mzmm), method = "complete")
    mzmm.melt$mm_name = factor(mzmm.melt$mm_name, levels = mouse.order$labels[mouse.order$order])
  }
}

# Custom string for the filename of the plot
ge_str = ifelse(isGE, "_ge3", "")
plot_fname = paste0("~/research/st/results/", mz.dataset, "_", mm.dataset,  ge_str, "_mine4.pdf")

# Plot the celltype mapping
pdf(plot_fname, width = (nrow(mzmm)/5) + 2, height = (ncol(mzmm)/5) + 2)
ggplot(mzmm.melt, aes(x = mm_name, y = mz_name, fill = Score)) + geom_raster() + scale_fill_gradientn(colors = col.pal, breaks = c(min(mzmm.melt$Score), max(mzmm.melt$Score))) + theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size = 10), axis.line=element_blank()) + scale_x_discrete(expand=c(0,0), name="") + scale_y_discrete(expand=c(0,0), name="") + coord_fixed() + force_panelsizes(cols = unit(nrow(mzmm)/8, "in"), rows = unit(ncol(mzmm)/8, "in")) +  geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1.4, color = "black") + geom_point(data = mzmm.melt[which(mzmm.melt$p0),], size = 1, color = "white")
dev.off()

# Save the dataframe with the scores and the p-values
write.csv(mzmm.melt[,c("mz_name", "mm_name", "Score", "p_perm")], paste0("~/research/st/results/supplemental/samap/", mz.dataset, "_", mm.dataset, ".csv"))
