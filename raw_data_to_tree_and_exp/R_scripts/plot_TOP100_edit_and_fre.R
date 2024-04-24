suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(ggpubr)
  library(ggstance)
  library(argparse)
})

plotT100Edit <- function(plot_path, allele_path, sample.name){
  plot_df <- read.csv(plot_path, sep = '\t', stringsAsFactors = F)
  alleleInfo <- read.csv(allele_path, sep = '\t', stringsAsFactors = F)
  ### input the annotate file
  plot_annotate <- read.table("/mnt/data/home/lzz/project/2019-8-22-PacBio.SampleA/material/new_annotate.txt", header = T, sep = '\t')
  linker <- plot_annotate[plot_annotate$type == "linker",]
  barcode <- plot_annotate[plot_annotate$type == "barcode",]
  pam <- plot_annotate[plot_annotate$type == "pam",]
  # 1. get ratios
  # make alleleinfo unique
  if ("cellNum" %in% colnames(alleleInfo)){
    colnames(alleleInfo)[colnames(alleleInfo) == "cellNum"] <- "Frequency"
  }
  alleleInfo <- alleleInfo %>% group_by(!!!syms(paste0("target", seq(1,13)))) %>% do({
    df <- .
    df[1,]
  })
  
  alleleInfo <- alleleInfo[with(alleleInfo, order(-Frequency)),]
  allreads <- sum(alleleInfo$Frequency)
  alleleInfo$ratio <- alleleInfo$Frequency / allreads
  rownames(alleleInfo) <- c(1:nrow(alleleInfo))
  # 2. get top 100 allele name
  plot_edit_allele <- alleleInfo$NodeName
  plot_edit_allele <- plot_edit_allele[plot_edit_allele != "WT"]
  plot_edit_allele <- plot_edit_allele[1:100]
  # 3. get top 100 allele plots and frequency
  plot_edit_df <- plot_df %>% filter(y %in% plot_edit_allele)
  plot_freq <- filter(alleleInfo, NodeName %in% plot_edit_allele)
  # 4. scale y axis
  plot_freq <- plot_freq[with(plot_freq, order(Frequency)),]
  sort_level <- plot_freq$NodeName %>% unlist() %>% as.character()
  plot_freq$NodeName <- factor(plot_freq$NodeName, levels = plot_freq$NodeName)
  fre_scale <- plot_freq$NodeName
  plot_edit_df$y <- factor(plot_edit_df$y, levels = sort_level)
  levels(plot_edit_df$y) <- sort_level
  plot_edit_df$start <- as.numeric(plot_edit_df$start)
  plot_edit_df$end <- as.numeric(plot_edit_df$end)
  # plot 
  edit_plot <- ggplot() +
    annotate("rect", xmin = linker$start, xmax = linker$end+1, ymin = 0.5, ymax = length(plot_edit_allele) + 0.5, alpha = linker$alpha, fill = "grey") +
    annotate("rect", xmin = barcode$start, xmax = barcode$end+1, ymin = 0.5, ymax = length(plot_edit_allele) + 0.5, alpha = barcode$alpha, fill = "white") +
    annotate("rect", xmin = pam$start, xmax = pam$end+1, ymin = 0.5, ymax = length(plot_edit_allele) + 0.5, alpha = pam$alpha, fill = "grey") +
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
  
  freq_plot <- ggplot(plot_freq, mapping = aes(x=ratio, y=NodeName)) + geom_barh(stat="identity", width = 1) +
    labs(x=NULL, y=NULL) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_discrete(breaks = seq(1,length(fre_scale)), labels = fre_scale, expand = c(0,0.1)) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0,0.5,0,0), "lines"))
  
  combin_plot <- ggpubr::ggarrange(edit_plot, freq_plot, ncol=2, common.legend = TRUE, legend = "bottom", widths = c(2,0.7), vjust = 0)
  combin_plot <- annotate_figure(combin_plot, top = text_grob(sample.name,color = "black", face = "bold", size = 14),)
  return(combin_plot)
}

## accept argument from command line
parser <- ArgumentParser(description = "plot a tree edit, and save as pdf.")
parser$add_argument("--ai", type="character", help="PATH of alleleinfo file.")
parser$add_argument("--pl", type="character", help="PATH of plot barcode edit file.")
parser$add_argument("-n", "--name", type="character", help="sample name")
parser$add_argument("--pdf", type="character", help="PATH of out put PDF file.")

useArgs <- parser$parse_args()

# run

combn_plot <- plotT100Edit(useArgs$pl, useArgs$ai, useArgs$name)

pdf(useArgs$pdf, width = 8, height = 4)
print(combn_plot)
dev.off()

