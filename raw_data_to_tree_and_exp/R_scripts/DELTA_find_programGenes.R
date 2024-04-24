# load pkgs
library(argparse)
library(dplyr)
library(ggvita)
# load some self functions
source("~/project/2019-8-22-PacBio.SampleA/scripts/R_scripts/DELTA_analysis_use_function.R")
source("~/project/simulatetion_Tree/scripts/use_function.R")
## accept argument from command line
parser <- ArgumentParser(description = "Analysis program gene use DELTA.")
parser$add_argument("--binExp", type="character", help="PATH of onecell node binary exp RData, the binary exp is 0 if gene no express in cell; 1 if gene express in cell")
parser$add_argument("--lineage", type="character", help="PATH of one cell node no class lineage file")
parser$add_argument("--score", type="character", help="PATH of DELTA align score", 
                    default="~/project/2019-8-22-PacBio.SampleA/results/m54061_190813_092847/DELTA_results/newDELTA.score.txt")
parser$add_argument("--tarNum", type="integer", help="Number of DELTA local alignment results wants to output.")
parser$add_argument("--out", type="character", help="PATH of output results RData")

args <- parser$parse_args()
# load all score results
oneCellNode.Expdata <- readRDS(args$binExp)
# node lineage
oneCell.Lineage <- read.csv(args$lineage, sep = '\t', stringsAsFactors = F, colClasses = "character")
# DELTA score
newScore <- args$score
## some paramters
testPvalue = 20
pruneScore = 2
targetNum = args$tarNum

checkGenes <- names(oneCellNode.Expdata)

target.gene.pvalue.df <- lapply(seq(length(checkGenes)), function(i){
  ## get the gene alm
  gene <- checkGenes[i]
  gene.lineage.path = tempfile(pattern = paste0("gene.lineage", gene), fileext = ".alm") # tmp file
  getGen.alm(oneCell.Lineage, oneCellNode.Expdata, gene, gene.lineage.path)
  # align 
  gene.align.path <- tempfile(pattern = paste0("gene.align", gene), fileext = ".alml") # tmp file
  selfDELTA.align(almFile = gene.lineage.path, almlFile = gene.align.path, scoreFile = newScore, test = NULL, prune = pruneScore, target = targetNum)
  # get total score
  #raw.total.score <- getTotalScore(gene.align.path)
  # permutate
  raw.align <- readal.alml(gene.align.path)
  dedup.raw.score <- getScoreDF(raw.align) %>% `[`("Score") %>% unlist() %>% as.numeric()
  #score_string <- paste0(dedup.raw.score, collapse = "_")
  # get permutate score and t.test
  message("run on NO.", i, "gene, gene name: ", gene)
  ## permutate result dir
  #dir.create(paste0(test_dir, gene, "/"), recursive = T)
  test.gene.pvalue <- mclapply(seq(100), function(j){
    times <- j
    random.rawlineage <- tempfile(pattern = paste0(gene, "random.leaves", times), fileext = ".alm") # tmp file
    ##
    getRandomTree(rawlineagePath = gene.lineage.path, randomLinPath = random.rawlineage)
    ## self to self local align
    random.trim.align <- tempfile(pattern = paste0("random.local.align", times), fileext = ".alml") # tmp file
    selfDELTA.align(almFile = random.rawlineage, almlFile = random.trim.align, scoreFile = newScore, test = NULL, prune = pruneScore, target = targetNum)
    ## read local align results and do t.test
    permutate.align <- readal.alml(random.trim.align)
    # get raw score
    # get random score
    dedup.random.score <- getScoreDF(permutate.align) %>% `[`("Score") %>% unlist() %>% as.numeric()
    #random_score_string <- paste0(dedup.random.score, collapse = "_")
    ## t test
    p.value <- ks.test(dedup.raw.score, dedup.random.score, alternative = "less")$p.value
    unlink(c(random.rawlineage, random.trim.align), recursive = T)
    return(p.value)
  }, mc.cores = 20) %>% unlist()
  unlink(c(gene.lineage.path, gene.align.path), recursive = T)
  return(data.frame(gene = gene, p_value = test.gene.pvalue))
}) %>% bind_rows(.)

saveRDS(target.gene.pvalue.df, args$out)
message("run complete!")
