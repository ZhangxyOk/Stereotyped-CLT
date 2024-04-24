#.libPaths(c("/mnt/data/home/lzz/R/x86_64-redhat-linux-gnu-library/3.6", "/usr/lib64/R/library", "/usr/share/R/library" ))
library(argparse)
library(ggtree)
library(ggplot2)
library(stringr)
library(dplyr)
library(ggstance)
#library(ggnewscale)

plotTreeEdit <- function(treeFile, EditPos, alleleInfo, plotfre=FALSE, noName=FALSE){
  tree <- read.tree(treeFile)
  editPos <- read.csv(EditPos, header = T, sep = '\t')
  Node.info.df <- read.csv(alleleInfo, header = T, sep = '\t')
  Node.info.df$NodeName <- str_trim(Node.info.df$NodeName)
  ## split INSERSION vs DELETION and NONE
  ## 
  if (noName){
    plotTree <- ggtree(tree, size=0.1, branch.length = "none") + geom_treescale(width=3, fontsize = 0, linesize = 0, color="white")
  } else {
    plotTree <- ggtree(tree, size=0.1, branch.length = "none") + geom_tiplab(size=0.8) + geom_treescale(width=3, fontsize = 0, linesize = 0, color="white")
  }
  # plot NONE and DELETION
  plotDN <- facet_plot(plotTree , panel = "Barcode", data=editPos, mapping=aes(x=NULL, y=NULL,xmin=start, xmax=end, ymin = y - 1, ymax = y + 1, fill=event),color="white",
                       size=0,geom = geom_rect, show.legend=TRUE) + scale_fill_manual(values = c("NONE"="#DEDEDE","DELETION"="red", "INSERTION"="blue")) + 
    theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), legend.text = element_text(size=10))
  ## plot Frequency
  if (!plotfre) return(plotDN)
  else{
    plotFre <- facet_plot(plotDN, panel = "Frequency", geom = geom_barh, data=Node.info.df, mapping=aes(x=Frequency),width = 0.4,fill="black", stat='identity')
    return(plotFre)
  }
}

## accept argument from command line
parser <- ArgumentParser(description = "plot a tree edit, and save as pdf.")
parser$add_argument("--tree", type="character", help="PATH of nwk format tree file.")
parser$add_argument("--ai", type="character", help="PATH of alleleinfo file.")
parser$add_argument("--pl", type="character", help="PATH of plot barcode edit file.")
parser$add_argument("--pdf", type="character", help="PATH of out put PDF file.")
parser$add_argument("--plotFre", action="store_true", default=FALSE, help="plot frequency.")
parser$add_argument("--noName", action="store_true", default=FALSE, help="don't plot node name")

args <- parser$parse_args()


treeplot <- plotTreeEdit(treeFile = args$tree, EditPos = args$pl, alleleInfo = args$ai, plotfre = args$plotFre, noName = args$noName)


pdf(args$pdf)
treeplot
dev.off()


## test













