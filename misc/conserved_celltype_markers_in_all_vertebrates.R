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
library("pheatmap")
library("scales")
library("viridis")

# Load markers
deg = list()
vert = sort(c("mouse", "axolotl", "bird", "turtle", "cichlid"))
deg[["mouse"]]   = read.csv("~/scratch/bcs/results/l5_cluster_markers_082823.csv")
deg[["axolotl"]] = read.csv("~/scratch/bcs/results/axolotl_cluster_markers_022123.csv")
deg[["turtle"]]  = read.csv("~/scratch/bcs/results/turtle_cluster_markers_041923.csv")
deg[["bird"]]    = read.csv("~/scratch/bcs/results/bird_markers_042823.csv")
deg[["cichlid"]] = read.csv("~/scratch/brain/results/bb53_deg_012323.csv")
deg[["cichlid"]]$cluster = stringr::str_replace(as.character(as.vector(deg[["cichlid"]]$cluster)), "Astro", "RG")

# Gene Orthologs
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
length(which(ortho[["mouse"]]$hgnc %in% ortho[["axolotl"]]$hgnc & ortho[["mouse"]]$hgnc %in% ortho[["turtle"]]$hgnc & ortho[["mouse"]]$hgnc %in% ortho[["bird"]]$hgnc & ortho[["mouse"]]$hgnc %in% ortho[["cichlid"]]$hgnc))

# Remove duplicate orthologs
ortho2 = list()
for (v1 in vert) {
  if (v1 == "cichlid") {
    ortho2[[v1]] = ortho[[v1]][which(ortho[[v1]]$gene %in% rownames(v[[v1]]@assays$RNA@counts)),]
  } else {
    ortho2[[v1]] = ortho[[v1]][which(ortho[[v1]]$gene %in% rownames(v[[v1]]@assays$SCT@counts)),]
  }
  ortho2[[v1]] = ortho2[[v1]][which(!duplicated(paste0(ortho2[[v1]]$gene, "_", ortho2[[v1]]$hgnc))),]
  ortho2[[v1]]$dup = duplicated(ortho2[[v1]]$gene) | duplicated(ortho2[[v1]]$gene, fromLast = TRUE)
  ortho2[[v1]]$dup[which(ortho2[[v1]]$dup & ortho2[[v1]]$hgnc == toupper(ortho2[[v1]]$gene))] = FALSE
  ortho2[[v1]] = ortho2[[v1]][which(!ortho2[[v1]]$dup),]
}
common_gene = ortho2[["mouse"]]$hgnc[which(ortho2[["mouse"]]$hgnc %in% ortho2[["axolotl"]]$hgnc & ortho2[["mouse"]]$hgnc %in% ortho2[["turtle"]]$hgnc & ortho2[["mouse"]]$hgnc %in% ortho2[["bird"]]$hgnc & ortho2[["mouse"]]$hgnc %in% ortho2[["cichlid"]]$hgnc)]
length(common_gene)

# Conserved Celltypes
common_name = c("OLIG", "OB", "MSN", "MGE1", "MGE2", "DGNBL", "CA3")
common = list()
common[[common_name[1]]] = data.frame(name = common_name[1], celltype = c("MFOL1", "OLIG15", "tsOlig", "Oligo", "2.2_Oligo"),
                                      species = c("mouse",  "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[2]]] = data.frame(name = common_name[2], celltype = c("OBDOP2", "GABA3", "i01", "GABA-1-1", "5.2_GABA"),
                                      species = c("mouse",  "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[3]]] = data.frame(name = common_name[3], celltype = c("MSN1", "GABA11", "i05", "MSN3", "4.1_GABA"),
                                      species = c("mouse",  "axolotl",  "turtle", "bird", "cichlid"))
common[[common_name[4]]] = data.frame(name = common_name[4], celltype = c("TEINH17", "GABA2", "i07", "GABA-3", "6_GABA"),
                                      species = c("mouse", "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[5]]] = data.frame(name = common_name[5], celltype = c("TEINH21", "GABA17", "i08", "GABA-2", "15.3_GABA"),
                                      species = c("mouse", "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[6]]] = data.frame(name = common_name[6], celltype = c("DGNBL1", "NB1", "tsNPCs", "Pre-2", "9.5_Glut"),
                                      species = c("mouse",  "axolotl", "turtle", "bird", "cichlid"))
common[[common_name[7]]] = data.frame(name = common_name[7], celltype = c("TEGLU23", "GLUT7", "e34", "HVC_Glut-3", "8.9_Glut"),
                                      species = c("mouse", "axolotl",  "turtle", "bird", "cichlid"))
common_celltype = do.call('rbind', common)

# Find the common DEGs for all vertebrates
vert = c("mouse", "bird", "turtle", "axolotl", "cichlid")
common_celltype_genes = data.frame()
common_celltype_genes4 = data.frame()
for (name_i in common_name) {
  cluster1 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == "cichlid")]
  for (v1 in vert[which(vert != "cichlid")]) {
    print(v1)
    cluster2 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == v1)]
    
    df1 = deg[["cichlid"]]
    df2 = deg[[v1]]
    df1 = df1[which(df1$cluster %in% cluster1),]
    df2 = df2[which(df2$cluster %in% cluster2),]
    
    df1$hgnc = ortho[["cichlid"]]$hgnc[match(df1$gene, ortho[["cichlid"]]$gene)]
    df2$hgnc = ortho[[v1]]$hgnc[match(df2$gene, ortho[[v1]]$gene)]
    
    ovlp = unique(df1$hgnc[which(df1$hgnc %in% df2$hgnc & !is.na(df1$hgnc) & df1$hgnc != "")])
    if (length(ovlp) > 0) {
      common_celltype_genes = rbind(common_celltype_genes, data.frame(name = name_i, species = v1, gene = ovlp))
    }
  }
  gene_rep = table(common_celltype_genes$gene[which(common_celltype_genes$name == name_i)])
  print(table(gene_rep))
  if (length(which(gene_rep == 4)) > 0) { common_celltype_genes4 = rbind(common_celltype_genes4, data.frame(name = name_i, gene = names(gene_rep)[which(gene_rep == 4)])) }
}
write.csv(common_celltype_genes4, "~/scratch/bcs/results/vert_cons_genes.csv")

# Find the common DEGs for pairwise combinations of species
common_pairwise = data.frame()
for (name_i in common_name) {
  for (v1 in vert) {
    cluster1 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == v1)]
    for (v2 in vert[1:(which(vert == v1)-1)]) {
      cluster2 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == v2)]
      
      df1 = deg[[v1]]
      df2 = deg[[v2]]
      df1 = df1[which(df1$cluster %in% cluster1),]
      df2 = df2[which(df2$cluster %in% cluster2),]
      
      df1$hgnc = ortho[[v1]]$hgnc[match(df1$gene, ortho[[v1]]$gene)]
      df2$hgnc = ortho[[v2]]$hgnc[match(df2$gene, ortho[[v2]]$gene)]
      
      ovlp = unique(df1$hgnc[which(df1$hgnc %in% df2$hgnc & !is.na(df1$hgnc))])
      min_num = min(length(df1$hgnc[which(!is.na(df1$hgnc))]), length(df2$hgnc[which(!is.na(df2$hgnc))]))
      pct = length(ovlp) / min_num
      common_pairwise = rbind(common_pairwise, data.frame(v1 = v1, v2 = v2, name = name_i, pct = pct, num = length(ovlp)))
    }
  }
}
common_pairwise = common_pairwise[which(common_pairwise$v1 != common_pairwise$v2),]
common_pairwise$v1 = factor(common_pairwise$v1, vert)
common_pairwise$v2 = factor(common_pairwise$v2, vert[1:(length(vert)-1)])
pdf(paste0("~/scratch/bcs/results/vert_ovlp_pct5.pdf"), width = 8, height = 10.67)
ggplot(common_pairwise, aes(x = v1, y = v2, size = pmax(pmin(num, 250),0), color = pct)) + geom_point() + facet_wrap(~ name, ncol = 2) + theme_classic() + theme(axis.text.x = element_text(size = 10)) + coord_fixed() + scale_color_viridis() + xlab("") + ylab("")
dev.off()

# Visualize the stats for the common DEGs across all vertebrates
common_celltype_genes_stats = data.frame()
for (name_i in common_name) {
  for (v1 in vert) {
    cluster1 = common_celltype$celltype[which(common_celltype$name == name_i & common_celltype$species == v1)]
    
    df1 = deg[[v1]]
    df1 = df1[which(df1$cluster %in% cluster1),]
    df1$hgnc = ortho[[v1]]$hgnc[match(df1$gene, ortho[[v1]]$gene)]
    
    high_cons_genes = common_celltype_genes4$gene[which(common_celltype_genes4$name == name_i)]
    df1 = df1[which(df1$hgnc %in% high_cons_genes),]
    df1 = df1[order(df1$p_val_adj, decreasing = F),]
    df1 = df1[which(!duplicated(df1$hgnc)),]
    common_celltype_genes_stats = rbind(common_celltype_genes_stats, data.frame(name = name_i, species = v1, gene = df1$hgnc, p_val_adj = df1$p_val_adj, log2FC = df1$avg_log2FC))
  }
}
neg_log_bon_thresh = 300
common_celltype_genes_stats$neg_log_bon = -log10(common_celltype_genes_stats$p_val_adj)
common_celltype_genes_stats$neg_log_bon[which(common_celltype_genes_stats$neg_log_bon > neg_log_bon_thresh)] = neg_log_bon_thresh
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "cichlid")] = "goldenrod1"
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "mouse")] = "#006d2c"
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "axolotl")] = "#bb3665"
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "turtle")] = "#004034"
common_celltype_genes_stats$species[which(common_celltype_genes_stats$species == "bird")] = "#CD700E"
pdf(paste0("~/scratch/bcs/results/vert_ovlp_genes_9422.pdf"), width = 8, height = 10.67)
ggplot(common_celltype_genes_stats, aes(x = neg_log_bon, y = gene, color = species, size = log2FC)) + geom_point() + scale_color_identity() + theme_bw() + facet_grid(name ~ ., scales = "free_y", space = "free_y") + theme(strip.text.x = element_blank()) + scale_x_continuous(limits=c(0, neg_log_bon_thresh+10))
dev.off()


# Common Gene Patterns =========================================================
# Find genes that exhibits similar patterns of gene expression, but may not strictly be DEGs

# Load the datasets
v = list()
v[["cichlid"]] = readRDS("~/scratch/brain/data/bb_demux_102021.rds")
v[["mouse"]]   = readRDS("~/scratch/bcs/data/l5_tel_all_sct.rds")
v[["axolotl"]] = readRDS("~/scratch/bcs/data/axolotl_norm_w_cluster.rds")
v[["turtle"]]  = readRDS("~/scratch/bcs/data/turtle_norm_cp_detail.rds")
v[["bird"]]    = readRDS("~/scratch/bcs/data/bird_smallx_033023_norm.rds")

convert53 = read.csv("~/scratch/st/data/convert53.csv")
v[["cichlid"]]$good_names53 = factor(convert53$new[match(v[["cichlid"]]$seurat_clusters, convert53$old)], levels = rev(convert53$new))

Idents(v[["cichlid"]]) = "good_names53"
Idents(v[["mouse"]])   = v[["mouse"]]$ClusterName
Idents(v[["axolotl"]]) = v[["axolotl"]]$cluster
Idents(v[["turtle"]])  = v[["turtle"]]$cp_detail
Idents(v[["bird"]])    = v[["bird"]]$cluster_orig2

# Find the foldchange of genes in celltypes
findSpeciesFC = function(i) {
  c1 = common_celltype$name[i]
  v1 = common_celltype$species[i]
  this_ident = common_celltype$celltype[i]
  print(paste0('common celltype = ', c1, ', species = ', v1, ", species celltype = ", this_ident))
  
  this_genes = ortho2[[v1]]$gene[which(ortho2[[v1]]$hgnc %in% common_gene)]
  this_fc = FoldChange(v[[v1]], features = this_genes, ident.1 = this_ident)
  this_fc$gene = rownames(this_fc)
  this_fc$hgnc = ortho2[[v1]]$hgnc[match(this_fc$gene, ortho2[[v1]]$gene)]
  this_fc$name = c1
  this_fc$species = v1
  this_fc$ident = this_ident
  print(paste0('*DONE* common celltype = ', c1, ', species = ', v1, ", species celltype = ", this_ident))
  return(this_fc)
}

v_common_fc = parallel::mclapply(1:nrow(common_celltype), function(x) findSpeciesFC(x), mc.cores = 1)
fc_df = do.call('rbind', v_common_fc)

hgnc_counts = fc_df[which(fc_df$name == "OLIG"),]
hgnc_counts = table(hgnc_counts$hgnc, hgnc_counts$species)
one_hgnc   = rownames(hgnc_counts)[which(rowSums(hgnc_counts) == 5)] 
multi_hgnc = rownames(hgnc_counts)[which(rowSums(hgnc_counts) != 5 & rowSums(hgnc_counts > 0) == 5)]
one_common = data.frame()
for (c1 in common_name) {
  this_one = data.frame(hgnc = c(one_hgnc, multi_hgnc), isOne = c(rep("TRUE", length(one_hgnc)), rep("FALSE", length(multi_hgnc))))
  for (v1 in vert) {
    this_df = fc_df[which(fc_df$name == c1 & fc_df$species == v1 & fc_df$hgnc %in% c(one_hgnc, multi_hgnc)),]
    multi_df = this_df[order(this_df$avg_log2FC, decreasing = T),]
    multi_df = multi_df[which(multi_df$hgnc %in% c(one_hgnc, multi_hgnc)),]
    multi_df = multi_df[which(!duplicated(multi_df$hgnc)),]
    multi_df = multi_df[match(c(one_hgnc, multi_hgnc), multi_df$hgnc),]
    this_one = cbind(this_one, multi_df$avg_log2FC)
  }
  colnames(this_one) = c("hgnc", "isOne", vert)
  this_one$name = c1
  one_common = rbind(one_common, this_one)
}
one_common$same_sign = rowSums(one_common[,vert] > 0)
table(one_common$same_sign, one_common$name)
one_common2 = one_common[which(one_common$same_sign == 5 & rowSums(one_common[,vert] > 0.1) == 5),]
one_common2$isDEG = one_common2$hgnc %in% common_celltype_genes4$gene
write.csv(one_common2, "~/scratch/bcs/results/vert_cons_genes_loose_090423.csv")

# Plot the foldchange
all_hit = one_common2
rownames(all_hit) = paste0("G", 1:nrow(all_hit))
this_other = -all_hit[,vert]
colnames(this_other) = paste0(colnames(this_other), "_other")
all_hit = cbind(all_hit[,vert], this_other)
all_hit = as.matrix(all_hit)
fc_thresh = 1
all_hit[which(all_hit >  fc_thresh)] =  fc_thresh
all_hit[which(all_hit < -fc_thresh)] = -fc_thresh
gene_annot = data.frame(annot = one_common2$name, row.names = rownames(all_hit))
pheatmap::pheatmap(all_hit, annotation_row = gene_annot, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, filename = paste0("~/scratch/st/results/loose_cons_genes_hit.pdf"))

hit_genes = one_common2$hgnc
hit_combos = paste0(one_common2$hgnc, "_", one_common2$name)
for (v1 in vert) {
  this_df = one_common[which(one_common$hgnc %in% hit_genes), c("hgnc", "name", v1)]
  this_mat = reshape2::acast(this_df, hgnc ~ name, value.var = v1)
  this_mat2 = this_mat[hit_genes,]
  rownames(this_mat2) = 1:nrow(this_mat2)
  this_df2 = reshape2::melt(this_mat2)
  this_df2$Var1 = factor(this_df2$Var1, levels = 1:nrow(this_mat2))
  print(ggplot(this_df2, aes(y = Var1, x = Var2, fill = value)) + geom_raster() + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), limits = c(-2, 2), oob = scales::squish ) + scale_y_discrete(labels = hit_genes, name = "", expand = c(0,0)) + scale_x_discrete(name = "", expand = c(0,0)) + theme_classic() + theme(axis.line=element_blank(), axis.text.y = element_text(face = ifelse(one_common2$hgnc %in% common_celltype_genes4$gene, "bold", "plain"))))
  ggsave(paste0("~/scratch/st/results/loose_cons_genes_hit_", v1, "_090423.pdf"), width = 5, height = 16)
}
