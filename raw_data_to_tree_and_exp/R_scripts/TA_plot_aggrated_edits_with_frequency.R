suppressMessages({
  library(tidyverse);
  library(reshape2);
  library(ggpubr);
  library(ggstance);
  library(argparse)
})

### ================== define some functions =====================
plotEdit.fre <- function(plot_df, count_df){
  ### input the annotate file
  plot_annotate <- read.table("~/project/Reconstruct_lineage/reference/plot_annotate.txt", header = T, sep = '\t')
  linker <- plot_annotate[plot_annotate$type == "linker",]
  barcode <- plot_annotate[plot_annotate$type == "barcode",]
  pam <- plot_annotate[plot_annotate$type == "pam",]
  # 3. get top 100 allele plots and frequency
  plot_edit_df <- plot_df
  plot_freq <- count_df
  # 4. scale y axis
  sort_level <- plot_freq$sample %>% unlist() %>% as.character()
  plot_freq$sample <- factor(plot_freq$sample, levels = plot_freq$sample)
  fre_scale <- plot_freq$sample
  plot_edit_df$y <- factor(plot_edit_df$y, levels = sort_level)
  levels(plot_edit_df$y) <- sort_level
  plot_edit_df$start <- as.numeric(plot_edit_df$start)
  plot_edit_df$end <- as.numeric(plot_edit_df$end)
  # plot 
  edit_plot <- ggplot() +
    annotate("rect", xmin = linker$start, xmax = linker$end+1, ymin = 0.5, ymax = nrow(count_df) + 0.5, alpha = linker$alpha, fill = "grey") +
    annotate("rect", xmin = barcode$start, xmax = barcode$end+1, ymin = 0.5, ymax = nrow(count_df) + 0.5, alpha = barcode$alpha, fill = "white") +
    annotate("rect", xmin = pam$start, xmax = pam$end+1, ymin = 0.5, ymax = nrow(count_df) + 0.5, alpha = pam$alpha, fill = "grey") +
    geom_rect(data = plot_edit_df, mapping = aes(x=NULL, y=NULL, xmin=start, xmax=end, ymin=as.numeric(y)-0.3, ymax=as.numeric(y)+0.3, fill=event), alpha=0.7) +
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
          #axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          plot.margin = unit(c(0,0,0.8,0), "lines"))
  
  freq_plot <- ggplot(plot_freq, mapping = aes(x=counts, y=sample)) + geom_barh(stat="identity", width = 1) +
    labs(x=NULL, y=NULL) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_discrete(breaks = seq(1,length(fre_scale)), labels = fre_scale, expand = c(0,0.1)) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(0,0.5,0,0), "lines"))
  
  combin_plot <- ggpubr::ggarrange(edit_plot, freq_plot, ncol=2, common.legend = TRUE, legend = "bottom", widths = c(2,0.7), vjust = 0)
  #combin_plot <- annotate_figure(combin_plot, top = text_grob(sample.name,color = "black", face = "bold", size = 14),)
  return(combin_plot)
}


## =================== accept argument from command line ============================
parser <- ArgumentParser(description = "Analysis TA clone data, merge same edits and count.")
parser$add_argument("-e", "--events", type="character", help="PATH of raw events file")
parser$add_argument("-p", "--plots", type="character", help="PATH of plots positions file")
parser$add_argument("--split", type="character", default="TRUE", help="split by sample index or not, TRUE or FALSE")
parser$add_argument("-o", "--out", type="character", help="PATH of output plot pdf")

useArgs <- parser$parse_args()


## ============ run ============
raw.events.df <- read.csv(useArgs$events, sep = '\t', stringsAsFactors = F)
plot.df <- read.csv(useArgs$plots, sep = '\t', stringsAsFactors = F)

raw.events.df <- raw.events.df %>% filter(., !grepl("primers", name))

raw.events.df$sample.name <- str_split(raw.events.df$name, "-", simplify = T) %>% `[`(,1)

if (useArgs$split == "TRUE") {
  events.count <- 
    raw.events.df %>% group_by(sample.name, !!!syms(paste0("target", seq(13)))) %>% 
    dplyr::summarise(counts=n(), name=sample(name, 1))
} else {
  events.count <- 
    raw.events.df %>% group_by(!!!syms(paste0("target", seq(13)))) %>% 
    dplyr::summarise(counts=n(), name=sample(name, 1))
  
}


use.plot.df <- plot.df %>% filter(., y %in% events.count$name)

counts.df <- events.count[,c("name", "counts")]
colnames(counts.df)[1] <- "sample"
#events.count$name <- 

ta.edit.count.p <- 
  plotEdit.fre(plot_df = use.plot.df,
               count_df = counts.df)

pdf(useArgs$out)
ta.edit.count.p
dev.off()



