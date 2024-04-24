library(Seurat)
library(dplyr)
library(argparse)
## import some function
source("~/project/simulatetion_Tree/scripts/use_function.R")
## some function
source("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/scripts/R_scripts/DELTA_analysis_use_function.R")

## accept argument from command line
parser <- ArgumentParser(description = "Analysis raw data, get the material for further analysis.")
parser$add_argument("--sc10xPath", type="character", help="PATH of 10x cellranger output feature bc matrix dir")
parser$add_argument("--alleleInfo", type="character", help="PATH of alleleinfo file")
parser$add_argument("--outSeurat", type="character", help="PATH of output Seurat object RData")
parser$add_argument("--outCounts", type="character", help="PATH of output node raw counts RData")
parser$add_argument("--outNorExp", type="character", help="PATH of output node normalize expression RData")
parser$add_argument("--outBinExp", type="character", help="PATH of output one cell node binary expression RData, use to DELTA find program genes")

args <- parser$parse_args()


#### load 10x data and get the one cell node raw counts and normalization expression
# 1. load 10x data and create Seurat object
sc293data <- Read10X(args$sc10xPath)
seuratObj <- CreateSeuratObject(sc293data, min.cells = 3, min.features = 200, project = "10x293T")
# filter mitochondrial
seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
seuratObj <- subset(seuratObj, subset = percent.mt < 5)

seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)

# 2. fine highly variable genes
seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)

# 3. get the exp data and save Seurat object
saveRDS(seuratObj, args$outSeurat)
seuratObj.hv.genes <- VariableFeatures(seuratObj)
seuratObj.norExpData <- as.matrix(GetAssayData(seuratObj[["RNA"]], slot = "data"))
seuratObj.rawCounts <- as.matrix(GetAssayData(seuratObj[["RNA"]], slot = "counts"))


# 4. get one cell node norExp and raw conuts
## norExpe
allele.info <- read.csv(args$alleleInfo, header = T, sep = '\t', stringsAsFactors = F)
# get oneCellNode normalize expression
oneCellNode.info <- filter(allele.info, cellNum == 1)
oneCellNode.norExpData <- seuratObj.norExpData[, oneCellNode.info$NodeName]
# get all node nor expression and output results
Node.norExpData <- seuratObj.norExpData[,]
colnames(Node.norExpData) <- allele.info$NodeName
saveRDS(Node.norExpData, args$outNorExp)

## raw counts
#oneCellNode.rawCounts <- seuratObj.rawCounts[, oneCellNode.info$BC]
Node.rawCounts <- seuratObj.rawCounts[,]
colnames(oneCellNode.rawCounts) <- oneCellNode.info$NodeName
saveRDS(oneCellNode.rawCounts, args$outCounts)

#### get the binary expression of one cell node, use to DELTA find program genes
t.oneCellNode.norExpData <- t(oneCellNode.norExpData)
t.oneCellNodeHV.norExpData <- t.oneCellNode.norExpData[,seuratObj.hv.genes]

t.remove0 <- apply(t.oneCellNodeHV.norExpData, 2, function(col) if(!all(col == 0)) return(col))
t.remove0 <- t.remove0[lengths(t.remove0) != 0]
t.remove0 <- as.data.frame(t.remove0)

binExp <- lapply(seq(ncol(t.remove0)), function(x){
  sub.df <- t.remove0[x]
  gene <- colnames(sub.df)
  node.names <- rownames(sub.df)
  if (!all(as.vector(unlist(sub.df)) == 0)){
    newExpBin <- ifelse(as.vector(unlist(sub.df)) > 0, 1, 0)
    sub.re <- data.frame(exp = newExpBin)
    colnames(sub.re) <- gene
    rownames(sub.re) <- node.names
    return(sub.re)
  } 
}) %>% bind_cols()

rownames(binExp) <- rownames(t.remove0)
saveRDS(binExp, args$outBinExp)









