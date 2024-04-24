#.libPaths(c("/mnt/data/home/lzz/R/x86_64-redhat-linux-gnu-library/3.6", "/usr/lib64/R/library", "/usr/share/R/library" ))
library(argparse)
library(ggtree)
library(ggplot2)
library(stringr)
library(dplyr)
library(ggstance)
library(ggnewscale)

plotTreeEdit <- function(treeFile, EditPos, alleleInfo, plotfre=FALSE, noName=FALSE){
  tree <- read.tree(treeFile)
  editPos <- read.csv(EditPos, header = T, sep = '\t')
  Node.info.df <- read.csv(alleleInfo, header = T, sep = '\t')
  Node.info.df$NodeName <- str_trim(Node.info.df$NodeName)
  ## split INSERSION vs DELETION and NONE
  inser.df <- editPos %>% filter(., event == "INSERTION")
  deleAndNone.df <- editPos %>% filter(., event == "DELETION" | event == "NONE")
  ## 
  if (noName){
    plotTree <- ggtree(tree, size=0.1, branch.length = "none") + geom_treescale(width=3, fontsize = 0, linesize = 0, color="white")
  } else {
    plotTree <- ggtree(tree, size=0.1, branch.length = "none") + geom_tiplab(size=0.8) + geom_treescale(width=3, fontsize = 0, linesize = 0, color="white")
  }
  # plot NONE and DELETION
  plotDN <- facet_plot(plotTree , panel = "Barcode", data=deleAndNone.df, mapping=aes(x=NULL, y=NULL,xmin=start, xmax=end, ymin = y - 0.5, ymax = y + 0.5, fill=event),color="white",
                       size=0,geom = geom_rect, show.legend=TRUE) + scale_fill_manual(values = c("NONE"="#DEDEDE","DELETION"="red")) + 
    theme(legend.position = "bottom", legend.key.size = unit(0.5, "cm"), legend.text = element_text(size=10))
  # plot INSERSION
  plotDN <- plotDN + new_scale_fill()
  plotIN <- facet_plot(plotDN, panel = "Barcode", data=inser.df, mapping = aes(x=NULL, y=NULL,xmin=start, xmax=end, ymin = y - 0.5, ymax = y + 0.5, fill=editLen), 
                       size=0,geom = geom_rect) +
    scale_fill_gradient(low="green", high="blue", name="INSER LEN")
  ## plot Frequency
  if (!plotfre) return(plotIN)
  else{
    plotFre <- facet_plot(plotIN, panel = "Frequency", geom = geom_barh, data=Node.info.df, mapping=aes(x=cellNum),width = 0.4,fill="black", stat='identity')
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

# if (args$plotFre == "TRUE" | args$plotFre == "T"){
#   treeplot <- plotTreeEdit(treeFile = args$tree, EditPos = args$pl, alleleInfo = args$ai, plotfre = TRUE)
# } else if (args$plotFre == "FALSE" | args$plotFre == "F") {
#   treeplot <- plotTreeEdit(treeFile = args$tree, EditPos = args$pl, alleleInfo = args$ai, plotfre = FALSE)
# } else {
#   stop("plotFre argument must be one of TRUE T FALSE F.")
# }
treeplot <- plotTreeEdit(treeFile = args$tree, EditPos = args$pl, alleleInfo = args$ai, plotfre = args$plotFre, noName = args$noName)


pdf(args$pdf)
treeplot
dev.off()


# ## 
# testBUG <- plotTreeEdit(treeFile = "~/project/2019-10-22-PacBio.SampleE/results/m54061_191012_053836/m54061_191012_053836.union.nwk",
#                         EditPos = "~/project/2019-10-22-PacBio.SampleE/results/m54061_191012_053836/m54061_191012_053836.union.Plots.txt",
#                         alleleInfo = "~/project/2019-10-22-PacBio.SampleE/results/m54061_191012_053836/m54061_191012_053836.union.AllelesInfo.txt")
# 
# 














