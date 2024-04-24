suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(ggpubr)
  library(ggstance)
  library(argparse)
})


plotAllEdit <- function(plot_path){
  plot_df <- read.csv(plot_path, sep = '\t', stringsAsFactors = F)
  ### input the annotate file
  plot_annotate <- read.table("~/project/Reconstruct_lineage/reference/plot_annotate.txt", header = T, sep = '\t')
  linker <- plot_annotate[plot_annotate$type == "linker",]
  barcode <- plot_annotate[plot_annotate$type == "barcode",]
  pam <- plot_annotate[plot_annotate$type == "pam",]
  # 3. get top 100 allele plots and frequency
  plot_edit_df <- plot_df 
  # 4. scale y axis
  all.record <- length(unique(plot_edit_df$y))
  plot_edit_df$y <- as.factor(plot_edit_df$y)
  plot_edit_df$start <- as.numeric(plot_edit_df$start)
  plot_edit_df$end <- as.numeric(plot_edit_df$end)
  # plot 
  edit_plot <- ggplot() +
    annotate("rect", xmin = linker$start, xmax = linker$end+1, ymin = 0.5, ymax = all.record + 0.5, alpha = linker$alpha, fill = "grey") +
    annotate("rect", xmin = barcode$start, xmax = barcode$end+1, ymin = 0.5, ymax = all.record + 0.5, alpha = barcode$alpha, fill = "white") +
    annotate("rect", xmin = pam$start, xmax = pam$end+1, ymin = 0.5, ymax = all.record + 0.5, alpha = pam$alpha, fill = "grey") +
    geom_rect(data = plot_edit_df, mapping = aes(x=NULL, y=NULL, xmin=start, xmax=end, ymin=as.numeric(y)-0.5, ymax=as.numeric(y)+0.5, fill=event), alpha=0.7) +
    scale_fill_manual(values = c("INSERTION"="blue","NONE"="#DEDEDE","DELETION"="red")) +
    labs(x=NULL, y=NULL) + 
    scale_y_continuous(breaks = seq(1,length(levels(plot_edit_df$y))), labels = levels(plot_edit_df$y), expand = c(0,0)) + 
    scale_x_continuous(expand = c(0,0)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          axis.line = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          plot.margin = unit(c(0,0,0.8,0), "lines"))
  
  return(edit_plot)
}

## accept argument from command line
parser <- ArgumentParser(description = "plot all edit, and save as pdf.")
parser$add_argument("--pl", type="character", help="PATH of plot barcode edit file.")
parser$add_argument("--pdf", type="character", help="PATH of out put PDF file.")

useArgs <- parser$parse_args()

# run

all_plot <- plotAllEdit(useArgs$pl)

pdf(useArgs$pdf, height = 56)
print(all_plot)
dev.off()