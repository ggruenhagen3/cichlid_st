# Load Libraries
library("Seurat")
library("SeuratObject")
library("SeuratDisk")
library("ggplot2")

# Load Data
Convert("~/scratch/bcs/data/bb_zeisel.h5ad", dest = "h5seurat", overwrite = TRUE)
merged = LoadH5Seurat("~/scratch/bcs/data/bb_zeisel.h5seurat", meta.data = FALSE, misc = FALSE) # cichlid + mouse data merged
mz = readRDS("~/scratch/brain/data/bb_demux_102021.rds") # cichlid data
merged$species = plyr::revalue(as.character(colnames(merged) %in% colnames(mz)), replace = c("TRUE" = "mz", "FALSE" = "mm"))
merged$species = factor(merged$species, levels = c("mm", "mz"))

# UMAP plot of cichlid + mouse data
Idents(merged) = merged$species
pdf("~/scratch/bcs/results/mz_mm_zei_umap.pdf", width = 10, height = 10)
print(DimPlot(merged, raster = F, order = T) + coord_fixed() + theme_bw() + theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Greens")[6], "goldenrod1")))
dev.off()
