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

# Load ST data
all_merge_hi = qs::qread(paste0("~/scratch/st/data/all_merge_hi.qs"))

# Load Cell2location results
c2l_mean = read.csv(paste0("~/scratch/st/results/cell2location_c2b2_spatial53_output_means.csv")); rownames(c2l_mean) = c2l_mean$X; c2l_mean$X = NULL; colnames(c2l_mean) = as.character(0:52)

# Convert the cluster numbers into cell type names
convert53 = read.csv("~/scratch/st/results/convert53.csv") # dataframe containing a column for cluster number and cell type names
colnames(c2l_mean) = convert53$new[match(colnames(c2l_mean), convert53$old)]; c2l_mean = c2l_mean[, convert53$new]

# Round the number of cells. Cells that round to 0, give them 1 count for the cell type closest to 1.
zero.cell.st = which(rowSums(round(c2l_mean)) == 0); for (i in zero.cell.st) { c2l_mean[i,which.max(c2l_mean[i,])] = 1 }
c2l_mean = round(c2l_mean)
c2l_mean_pct = c2l_mean / rowSums(c2l_mean)

# Get the celltype with the most cells per spot
st.celltype.char = unlist(lapply(1:nrow(c2l_mean), function(x) colnames(c2l_mean)[which.max(c2l_mean[x,])]))

# Plot the composition of ST clusters by cell types
c2l.to.cluster = expand.grid(cluster = levels(all_merge$cluster), celltype = colnames(c2l_mean))
sum.cells.per.cluster = unlist(lapply(levels(all_merge$cluster), function(x) sum(c2l_mean[colnames(all_merge)[which(all_merge$cluster == x)],]) ))
names(sum.cells.per.cluster) = levels(all_merge$cluster)
c2l.to.cluster$cluster.sum.num.cells = sum.cells.per.cluster[c2l.to.cluster$cluster]
c2l.to.cluster$num = unlist(lapply(1:nrow(c2l.to.cluster), function(x) sum(c2l_mean[colnames(all_merge)[which(all_merge$cluster == c2l.to.cluster$cluster[x])], c2l.to.cluster$celltype[x]]) ))
c2l.to.cluster$pct.of.cluster = c2l.to.cluster$num / c2l.to.cluster$cluster.sum.num.cells
pdf(paste0( "stcluster_to_celltype.pdf"), width = 8, height = 5)
ggplot(c2l.to.cluster, aes(x = celltype, y = cluster, fill = pct.of.cluster)) + geom_raster() + coord_fixed() + scale_x_discrete(expand = c(0,0), name = "") + scale_y_discrete(expand = c(0,0), name = "") + theme_classic() + scale_fill_viridis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.ticks = element_blank(), axis.line = element_blank())
dev.off()

# Plot the composition of ST brain regions by cell types
c2l.to.region = expand.grid(region = levels(all_merge$struct), celltype = colnames(c2l_mean))
sum.cells.per.region = unlist(lapply(levels(all_merge$struct), function(x) sum(c2l_mean[colnames(all_merge)[which(all_merge$struct == x)],]) ))
names(sum.cells.per.region) = levels(all_merge$struct)
c2l.to.region$region.sum.num.cells = sum.cells.per.region[c2l.to.region$region]
c2l.to.region$num = unlist(lapply(1:nrow(c2l.to.region), function(x) sum(c2l_mean[colnames(all_merge)[which(all_merge$struct == c2l.to.region$region[x])], c2l.to.region$celltype[x]]) ))
c2l.to.region$pct.of.region = c2l.to.region$num / c2l.to.region$region.sum.num.cells
c2l.to.region$maxed = c2l.to.region$pct.of.region * 100
c2l.to.region$maxed[which(c2l.to.region$maxed > 50)] = 50
c2l.to.region$col = convert53$col[match(c2l.to.region$celltype, convert53$new)]
col2 = scales::hue_pal()(nrow(convert53))
pdf(paste0( "stregion_to_celltype.pdf"), width = 6, height = 8)
ggplot(c2l.to.region, aes(x = region, y = celltype, fill = maxed)) + geom_raster() + coord_fixed() + scale_x_discrete(expand = c(0,0), name = "") + scale_y_discrete(expand = c(0,0), name = "") + theme_classic() + scale_fill_viridis() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), axis.text.y = element_text(size=10, color = col2), axis.line = element_blank()) + force_panelsizes(rows = unit(length(colnames(c2l_mean))/8, "in"), cols = unit(length(levels(all_merge$struct))/8, "in"))
dev.off()


