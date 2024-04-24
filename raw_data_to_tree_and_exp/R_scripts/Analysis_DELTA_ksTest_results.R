# load pkgs
library(argparse)
library(dplyr)

# load some self functions
source("~/project/2019-8-22-PacBio.SampleA/scripts/R_scripts/DELTA_analysis_use_function.R")
source("~/project/simulatetion_Tree/scripts/use_function.R")
## accept argument from command line
parser <- ArgumentParser(description = "Analysis DELTA ks test results.")
parser$add_argument("-i", "--inKSraw", type="character", help="PATH of DELTA ks test raw results RData")
parser$add_argument("-o", "--out", type="character", help="PATH of output results RData")

args <- parser$parse_args()
##
# 1. read in raw ks test result 
ks.test.Result <- readRDS(args$inKSraw)
# 2. adjust pvalue 
ks.test.padjust.df <- ks.test.Result %>% group_by(gene) %>% do(data.frame(p_adjust = p.adjust(.$p_value, method = "BH")))
# 3. calculate gm mean
ks.test.gm_mean <- ks.test.padjust.df %>% group_by(gene) %>% summarise(p_value.gm_mean = gm_mean(p_adjust))
# 4. get positive genes base on gm mean <= 0.05
ks.positive.pvalue.df <- filter(ks.test.gm_mean, p_value.gm_mean <= 0.05)
# 5. save results
resultList = list(raw.ks.test.result = ks.test.Result,
                  ks.test.adjust = ks.test.padjust.df,
                  ks.adjust.gm_mean = ks.test.gm_mean,
                  ks.positive = ks.positive.pvalue.df)

saveRDS(resultList, args$out)