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

# Load Zeisel Dataset
zei = readRDS("~/scratch/bcs/data/l5_tel_all.rds")

# Randomly Remove 10% of Counts
zei.counts = zei@assays$RNA@counts
zei.num.to.remove = round(zei$nCount_RNA * 0.10)
idx.to.rm = mclapply(1:ncol(zei.counts), function(x) { pos.idx = which(zei.counts[,x] > 0); rm.idx = sample(pos.idx, zei.num.to.remove[x]); return(1:nrow(zei.counts) %in% rm.idx) }, mc.cores = 24 )
idx.rm.mat = do.call('cbind', idx.to.rm)
class(idx.rm.mat) = "numeric"
zei.counts.down = zei.counts - idx.rm.mat
down.zei = CreateSeuratObject(counts = zei.counts.down, meta.data = zei@meta.data)

# Save Output
saveRDS(down.zei, "~/scratch/bcs/data/l5_tel_all_down_022823.rds")

