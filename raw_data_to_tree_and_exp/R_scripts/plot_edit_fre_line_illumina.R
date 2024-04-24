suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(ggpubr)
  library(argparse)
})



plot_edit_fre <- function(edit_fre_path){
  ### input the annotate file
  #plot_annotate <- read.table("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/material/new_annotate.txt", header = T, sep = '\t')
  plot_annotate <- read.table("~/project/Reconstruct_lineage/reference/plot_annotate.txt", header = T, sep = '\t')
  linker <- plot_annotate[plot_annotate$type == "linker",]
  barcode <- plot_annotate[plot_annotate$type == "barcode",]
  pam <- plot_annotate[plot_annotate$type == "pam",]
  ### plot edit percent
  EditFre <- read.table(edit_fre_path, header = T, sep = '\t')
  ## just extract the ecit percent
  subEditFre <- EditFre[,c(1,3,5)]
  aa <-subEditFre %>%
    melt(id.vars = ("Position"))
  p <- 
    aa %>%
    ggplot(.,aes(x=Position,y=value*100,color = variable))+ labs(x="Position", y= "edit/all")+
    annotate("rect", xmin = linker$start, xmax = linker$end+1, ymin = 0, ymax = round(max(aa$value)*100) + 2, alpha = linker$alpha, fill = "grey") +
    annotate("rect", xmin = barcode$start, xmax = barcode$end+1, ymin = 0, ymax = round(max(aa$value)*100) + 2, alpha = barcode$alpha, fill = "grey") +
    annotate("rect", xmin = pam$start, xmax = pam$end+1, ymin = 0, ymax = round(max(aa$value)*100) + 2, alpha = pam$alpha, fill = "grey") +
    geom_line() + scale_color_manual(values = c("InserPercent" = "blue", "DeletPercent" = "red")) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  return(p)
}



## accept argument from command line
parser <- ArgumentParser(description = "plot edit frequency, and save as pdf.")
parser$add_argument("--EditFreq", type="character", help="PATH of edit frequency file.")
parser$add_argument("--pdf", type="character", help="PATH of out put PDF file.")

useArgs <- parser$parse_args()

### run 

line.p <- 
  plot_edit_fre(useArgs$EditFreq)


pdf(useArgs$pdf, width = 14)
line.p
dev.off()









