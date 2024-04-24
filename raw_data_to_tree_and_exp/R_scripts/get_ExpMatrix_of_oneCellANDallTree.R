# load pkgs
library(argparse)
library(dplyr)

## accept argument from command line
parser <- ArgumentParser(description = "get oneCelll node exp matrix and all tree node exp matrix")
parser$add_argument("-i", "--inMatrix", type="character", help="PATH of cell exp matrix extract from Seurat, Rds")
parser$add_argument("-a", "--alleleInfo", type="character", help="PATH of alleleinfo file path")
parser$add_argument("--outOneCell", type="character", help="PATH of output oneCell node exp matrix Rds")
parser$add_argument("--outAllTree", type="character", help="PATH of output allTree node exp matrix Rds")

args <- parser$parse_args()

# run
exp.mx <- readRDS(args$inMatrix)

# read in allele info 
alleleInfo <- read.csv(args$alleleInfo, sep = '\t', stringsAsFactors = F)
# substitu the cell barcode to node name
bc_node.df <- alleleInfo[,c("NodeName", "BC")]
row.names(bc_node.df) <- bc_node.df$BC
colnames(exp.mx) <- bc_node.df[colnames(exp.mx), 1]
exp.mx <- exp.mx[, !is.na(colnames(exp.mx))]

####   get oneCell exp matrix
oneCell.exp <- exp.mx[, grepl("_1$", colnames(exp.mx))]

####  get all node exp matrix
#  get multi cell node exp matrix, the multi cell node expression is the mean of all cells in this node
multiCell.exp <- exp.mx[, !grepl("_1$", colnames(exp.mx))]
multiCellNode <- alleleInfo$NodeName %>% `[`(!grepl("_1$", alleleInfo$NodeName)) %>% unique()

mulcellExp <- lapply(seq(length(multiCellNode)), function(i){
  nodeName <- multiCellNode[i]
  pattern <- paste0("^", nodeName)
  nodeExp.df <- multiCell.exp[, grepl(pattern, colnames(multiCell.exp))]
  nodeExp.mean <- rowMeans(nodeExp.df) %>% data.frame()
  colnames(nodeExp.mean) <- nodeName
  rownames(nodeExp.mean) <- rownames(nodeExp.df)
  return(nodeExp.mean)
}) %>% bind_cols()
rownames(mulcellExp) <- rownames(multiCell.exp)
rm(multiCell.exp)
#  all tree node exp matrix (multicell node exp use their mean exp represent)
wholeNodeExp <- cbind(mulcellExp, oneCell.exp)
oneCell.exp <- oneCell.exp[rowSums(oneCell.exp) != 0, ]
wholeNodeExp <- wholeNodeExp[rowSums(wholeNodeExp) != 0,]

# save out results
saveRDS(oneCell.exp, args$outOneCell)
saveRDS(wholeNodeExp, args$outAllTree)









