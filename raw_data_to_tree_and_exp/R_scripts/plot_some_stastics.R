suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(reshape2)
  library(ggpubr)
  library(argparse)
})


plot_sta_results <- function(stafile){
  # top 5 event plot
  top5_event <- filter(stafile, class == "Top5_event")
  top5_event$event <- factor(top5_event$event, levels = top5_event$event)
  top5_event_plot <- top5_event %>% ggplot() + geom_bar(aes(x=event, y=frequency), stat = "identity") + 
    labs(title = "Top 5 edit events", x="") + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  # top 5 long dele
  top5_longDele <- filter(stafile, class == "Top5_long_dele")
  top5_longDele$event <- factor(top5_longDele$event, levels = top5_longDele$event)
  top5_longDele_plot <- top5_longDele %>% ggplot() + geom_bar(aes(x=event, y=frequency), stat = "identity") + 
    labs(title = "Top 5 long deletion", x="") + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  # long dele type fre
  long_dele_fre <- filter(stafile, class == "Long_dele_event")
  long_dele_fre$event <- factor(long_dele_fre$event, levels = long_dele_fre$event)
  long_dele_fre_plot <- long_dele_fre %>% ggplot() + geom_bar(aes(x=event, y=frequency), stat = "identity") + 
    labs(title = "Long deletion class", x="") + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank())
  # allele componten
  allele_com <- filter(stafile, class == "Allele_of_leaf")
  allele_com_plot <- allele_com %>% ggplot() + geom_bar(aes(x=class, y=frequency, fill=event), stat = "identity", width = 0.6) + 
    labs(title = "Allele component",
         x = "",
         y = "frequency") + 
    scale_fill_manual(values = c("Allele_no_long_dele"="deepskyblue", "Allele_with_long_dele"="red"), labels=c("without long delet", "with long delet")) + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.title = element_blank())
##
  return(ggarrange(allele_com_plot, top5_event_plot, top5_longDele_plot, long_dele_fre_plot,ncol = 2, nrow = 2, widths = c(2,2,2,2),common.legend = F))
}



## accept argument from command line
parser <- ArgumentParser(description = "plot a stastics infos, and save as pdf.")
parser$add_argument("-i", "--inSta", type="character", help="PATH of stastics file.")
parser$add_argument("-o", "--out", type="character", help="PATH of out put PDF file.")

useArgs <- parser$parse_args()

# run
# test

#test <- "/mnt/data/home/lzz/project/2021-4-9_lung_scRNA-seq/results/PacBio/run_by_strand/A1-CBRAD5/A1-CBRAD5.Statistics.txt"

plot_stafile <- read.csv(useArgs$inSta, sep = '\t', stringsAsFactors = F)

p <- plot_sta_results(plot_stafile)




pdf(useArgs$out, width = 7, height = 7)
print(p)
dev.off()





