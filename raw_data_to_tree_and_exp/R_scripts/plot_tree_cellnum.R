suppressMessages({
  library(argparse)
  library(ggtree)
  library(dplyr)
  library(ggplot2)
  library(treeio)
  library(stringr)
  library(scales)
})


plotTreeAndCellNum <- function(tree, allele.info.df){
  # get node cell count
  node.cell.count <- allele.info.df %>% dplyr::group_by(NodeName) %>% dplyr::summarise(fre=n())
  plot.node.cell.count <- lapply(seq(nrow(node.cell.count)), function(x){
    sub.df <- node.cell.count[x, ]
    fre <- sub.df$fre[1]
    color <- ifelse(fre == 1, 1, ifelse(fre >= 10, 10, 5))
    return(cbind(sub.df, color))
  }) %>% bind_rows(.)
  rownames(plot.node.cell.count) <- plot.node.cell.count$NodeName
  rownames(plot.node.cell.count) <- str_trim(rownames(plot.node.cell.count))
  # plot tree and count
  p <- ggtree(tree, branch.length="none" ,layout = "circular", size=0.1) + xlim(-10, NA) + scale_color_manual(values=c("black", "red")) 
  p <- rotate_tree(p, 80)
  plot.node.cell.count <- plot.node.cell.count[-1]
  plot.node.cell.count$color <- as.factor(plot.node.cell.count$color)
  plot.node.cell.count <- plot.node.cell.count[-1]
  p <- gheatmap(p, plot.node.cell.count, offset=.1, width=.1, colnames_angle=120, colnames_offset_y = .25, colnames = FALSE) + xlim(-25, NA) +
    scale_fill_discrete(name="cell number", labels=c("1 cell", "2~9 cells", ">=10 cells"),limits=c("1", "5", "10")) 
  return(p)
}


## accept argument from command line
parser <- ArgumentParser(description = "plot a tree and node cell numbers, and save as pdf.")
parser$add_argument("--ai", type="character", help="PATH of alleleinfo file.")
parser$add_argument("--tree", type="character", help="PATH of tree file.")
parser$add_argument("--pdf", type="character", help="PATH of out put PDF file.")

useArgs <- parser$parse_args()

### run 
tree.obj <- read.tree(useArgs$tree)
alleleinfo.df <- read.csv(useArgs$ai, sep = '\t', stringsAsFactors = F)

out.plot <- 
  plotTreeAndCellNum(tree = tree.obj,
                     allele.info.df = alleleinfo.df)

pdf(useArgs$pdf)
out.plot + theme(legend.position = c(0.5, 0.5))
dev.off()




