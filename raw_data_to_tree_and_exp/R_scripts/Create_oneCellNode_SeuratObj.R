suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(argparse)
})

## accept argument from command line
parser <- ArgumentParser(description = "get oneCelll node exp matrix and all tree node exp matrix")
parser$add_argument("-i", "--inSeuratObj", type="character", help="PATH of raw Seurat, Rds")
parser$add_argument("-a", "--alleleInfo", type="character", help="PATH of alleleinfo file path")
parser$add_argument("--outOneCellNodeSeuratObj", type="character", help="PATH of all raw Seurat object Rds")

args <- parser$parse_args()

### run
# read in allele info 
alleleInfo <- read.csv(args$alleleInfo, sep = '\t', stringsAsFactors = F)

all.raw.seuratObj <- readRDS(args$inSeuratObj)
### all raw counts matrix
all.raw.counts <- all.raw.seuratObj@assays$RNA@counts
colnames(all.raw.counts) <- gsub("*_|-1", "", colnames(all.raw.counts))
all.raw.counts <- all.raw.counts[Matrix::rowSums(all.raw.counts) != 0,]

### overlap raw counts matrix
oneCellNode.BC <- alleleInfo %>% filter(., cellNum == 1) %>% `[`("BC") %>% unlist()

keep.counts <- all.raw.counts[, oneCellNode.BC]
keep.counts <- keep.counts[Matrix::rowSums(keep.counts) != 0, ]

colnames(keep.counts) <- alleleInfo$NodeName[match(colnames(keep.counts), alleleInfo$BC)]

keep.raw.seurat.obj <- CreateSeuratObject(counts = keep.counts)

keep.raw.seurat.obj <- NormalizeData(keep.raw.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
keep.raw.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(keep.raw.seurat.obj, pattern = "^MT-")
## output results
saveRDS(keep.raw.seurat.obj, file = args$outOneCellNodeSeuratObj)
