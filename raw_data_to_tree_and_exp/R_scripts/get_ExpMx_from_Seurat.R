suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(argparse)
})


## accept argument from command line
parser <- ArgumentParser(description = "get oneCelll node exp matrix and all tree node exp matrix")
parser$add_argument("-i", "--inSeuratObj", type="character", help="PATH of raw Seurat, Rds")
parser$add_argument("-a", "--alleleInfo", type="character", help="PATH of alleleinfo file path")
parser$add_argument("--outRawSeuratObj", type="character", help="PATH of all raw Seurat object Rds")
parser$add_argument("--outOverlapRaw", type="character", help="PATH of ovarall raw exp matrix Rds")

args <- parser$parse_args()

### run
# read in allele info 
alleleInfo <- read.csv(args$alleleInfo, sep = '\t', stringsAsFactors = F)

all.raw.seuratObj <- readRDS(args$inSeuratObj)

### all raw counts matrix
all.raw.counts <- all.raw.seuratObj@assays$RNA@counts
colnames(all.raw.counts) <- gsub("*_|-1", "", colnames(all.raw.counts))
all.raw.counts <- all.raw.counts[Matrix::rowSums(all.raw.counts) != 0,]
new.raw.seurat.obj <- CreateSeuratObject(counts = all.raw.counts)
new.raw.seurat.obj <- NormalizeData(new.raw.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

new.raw.exp.mx <- as.matrix(new.raw.seurat.obj@assays$RNA@data)

### overlap raw counts matrix
exp.mx <- new.raw.exp.mx[, alleleInfo$BC]

## output results
saveRDS(new.raw.seurat.obj, file = args$outRawSeuratObj)
saveRDS(exp.mx, file = args$outOverlapRaw)




