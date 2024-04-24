library(tidyverse);
library(ape);
library(parallel);
library(plyr);
library(tidyverse);
library(parallel);
library(amap); ## amap::dist is a multithreaded version of "dist"
library(ape);
library(treeio);
library(readr);
library(reshape2);
library(magrittr);
library(ggplot2)

needNtips <- readRDS("~/figs5/needNtips.Rds");

expectedNdiv.nTip <- rep(NA,2000);  ## the distance of every tip pair was defined by Ndiv
nRep = 1000;
for(x in c(4:2000)) {
  cat(paste0(x,"\n"));
  res <- mclapply(c(1:nRep),mc.cores=20,function(i){
    return(sum(cophenetic(rtree(x,br=1))) / choose(x,2) / 2 - 1); #  calculate the mean distance of two tips, the sum of distance divided by combination number(choose function) and 2 (cophenetic symmetical matrix)
  }) %>%
    unlist() %>%
    mean();   # each tree of specific size from 4 to 2000 were circularly calculated for 1000 times and mean lastly
  expectedNdiv.nTip[x] <- res;
}

expectedNdiv.nTip[1] <- 0;
expectedNdiv.nTip[2] <- 1;
expectedNdiv.nTip[3] <- 5/3; 

saveRDS(expectedNdiv.nTip,file="nTips2nDiv.Rds")

expectedNdiv.nTip<-readRDS("~/figs5/nTips2nDiv.Rds")%>%as.data.frame()

expectedNdiv.nTip<-expectedNdiv.nTip%>%mutate(.,size=seq(nrow(expectedNdiv.nTip)))%>%select(.,expDist.pair=.,size=size)

# functions

## Function to estimate correlation between relatedness and trait difference (inspired by Haseman-Elston Regression)
## This function is compatible for the traits of internal nodes 
## (such as cell type component of tips under the internal nodes)
##   dfTree - a tbl_tree object of the tree structure
##   dfFeat - a data.frame of the features (expression, PC, etc). 
##            Column names should be available in dfTree$label.
##            Features of nodes in dfTree but unavailable in dfFeat will be assigned NA
##   nBootstrap - number of bootstrap assessing SD of correlation. Set to <2 to disable.
##   nPermutation - number of permutation of nodes assessing significance of correlation. Set to <1 to disable.
##   nProc - number of threads to use
##   relatednessMethod - method to estimate relatedness coefficient, use "mrcaNtip" or "nSeparatingDivision"
##   distMethod - method to pairwise phenotypic distances between nodes, use "euclidean" or "pearson"
##   mrcaNtip2sepDiv - a vector indicating the correspondance between mrcaNtip and nSeparatingDivision, data from code named"nTips2nDiv.R"

heritabilityByHEregression <- function(dfTree,
                                       dfFeat,
                                       nBootstrap,
                                       nPermutation,
                                       nProc,
                                       relatednessMethod,
                                       distMethod,
                                       mrcaNtip2sepDiv){
  if(dfTree %>% filter(parent == node) %>% nrow() != 1) {
    stop("The dfTree should contain one and only one root !");
  } 
  if( !is(dfTree,"tbl_tree") ) {
    stop("Expecting a tidytree::tbl_tree for the dfTree");
  } 
  if( ! (distMethod %in% c("pearson","euclidean")) ) {
    stop("distMethod must be one of \"pearson\" or \"euclidean\"");
  } 
  if( !any(match(colnames(dfFeat),dfTree$label) ) ) {
    stop("At least one node of dfTree should be found in the columns of dfFeat");
  } 
  if(length(dfTree$label) != length(unique(dfTree$label)) ) {
    stop("Node labels in dfTree is not unique!")
  } 
  if((relatednessMethod == "nSeparatingDivision") & missing(mrcaNtip2sepDiv)) {
    stop("mrcaNtip2sepDiv is required for relatednessMethod of \"nSeparatingDivision\"");
  }
  
  allTipNodes <- setdiff(dfTree$node,dfTree$parent);
  
  ## SECTION 1: pairwise trait distance
  nNodeInTree <- length(dfTree$node); # dfTree$node: tip nodes and internal nodes
  nFeat <- nrow(dfFeat);
  dfFeat2 <- as.tibble(data.frame(matrix(NA,nrow=nFeat,ncol=nNodeInTree)));
  colnames(dfFeat2) <- as.character(dfTree$label);
  matchNodes <- match(dfTree$label,colnames(dfFeat));  ## find index 
  dfFeat2[,!is.na(matchNodes)] <- dfFeat[,na.omit(matchNodes)];
  ## this ensures that it works even for features of the internal nodes when dfFeat2 was compared with dfFeat,but internal node exp is NA
  
  matDist <- amap::Dist(t(dfFeat2),method=distMethod,nbproc=nProc) %>% as.matrix(); ## this is not squared
  if(distMethod == "euclidean") {
    matDist <- matDist ** 2;  
  }  
  
  ## SECTION 2 : calculate relatedness coefficient
  ##
  dfTipUnderNode <-                              # function dfTipUnderNode and dfNodeSize were so functional
    mclapply(allTipNodes,
             mc.cores = nProc,
             function(thisTip) {
               tidytree::ancestor(dfTree,thisTip) %>%
                 mutate(tip = thisTip) %>%
                 select(tip,node) %>%
                 rbind.fill(data.frame(tip=thisTip,node=thisTip)) %>%
                 return();
             }) %>%
    rbind.fill() %>%
    merge(dfTree %>% select(node,label), by.x="tip", by.y="node");
  dfNodeSize <- dfTipUnderNode %>%
    group_by(node) %>%
    dplyr::summarise(size = length(tip)) %>%
    merge(dfTree %>% select(node,label),by.x="node",by.y="node")
  
  phyloTree <- treeio::as.phylo(dfTree);   ## as.phylo might rename the nodes, so, to be cautious...
  phylo.node2label <- c(phyloTree$tip.label,phyloTree$node.label);
  names(phylo.node2label) <- as_tibble(phyloTree)$node %>% as.character();
  allMrca <- mrca(phyloTree,full=T); ## result from "full = T" always named the dimensions by node instead of label.mrca:find most common ancestor of two nodes
  allMrca <- matrix(nrow = nNodeInTree,phylo.node2label[as.character(allMrca)]);
  colnames(allMrca) <- phylo.node2label;
  rownames(allMrca) <- phylo.node2label;
  
  matRelat <- ## relatedness coeffcient
    matrix(nrow = nNodeInTree,
           dfNodeSize$size[match(allMrca,dfNodeSize$label)] -
             (dfNodeSize$size[match(rep(colnames(allMrca),each = nNodeInTree),dfNodeSize$label)] - 1) -
             (dfNodeSize$size[match(rep(rownames(allMrca),nNodeInTree),dfNodeSize$label)] - 1));
  ## for above,this is just based on the number of tips.(such as:root to root is -1045,comb-like leaf to leaf is 1047,range of matRelat is -1045~1047)
  # The nSeparatingDivision is based on this. See below
  matRelat[ ( allMrca == rep(colnames(allMrca),each = nNodeInTree) ) |
             ( allMrca == rep(rownames(allMrca),nNodeInTree) ) ] <- NA; 
  ## for above,should not compare a node with its ancestor: diagonal values are NA,and leaves or inter-nodes in the row or col of its' ancestor are NA. 
  # range of matRelat is 2~1047 after NA assignment 
  if(relatednessMethod == "mrcaNtip") {
    matRelat <- 1/ (matRelat - 1); ## One need (n-1) divisions to have n tips, "n-1" equal to the number of internalnodes 
    diag(matRelat) <- 1;
    dimnames(matRelat) <- dimnames(allMrca);
  } else if(relatednessMethod == "nSeparatingDivision") {
    matRelat <- mrcaNtip2sepDiv[matRelat] %>% matrix(nrow = nNodeInTree); 
    ## for above,mrcaNtip2sepDiv is the mean distance of two tips for random different size tree,the distance of tip pair was defined by Ndiv.
    # different size random tree is created for remove the error with tree division number with the true cell division number. 
    diag(matRelat) <- 0;
    matRelat <- 1 / (matRelat + 1);
    dimnames(matRelat) <- dimnames(allMrca);
  } else if(relatednessMethod == "2pow.nSeparatingDivision") {
    matRelat <- mrcaNtip2sepDiv[matRelat] %>% matrix(nrow = nNodeInTree);
    diag(matRelat) <- 0;
    matRelat <- 1 / 2 ** matRelat;
    dimnames(matRelat) <- dimnames(allMrca);
  } 
  dimnames(matRelat) <- dimnames(allMrca);
    
  ## Section 3 : Haseman-Elston Regression
  allNodeLabels <- colnames(allMrca);
  dfHEreg <- allMrca %>%
    reshape2::melt(as.is=T) %>%
    dplyr::rename(node1 = Var1,
                  node2 = Var2,
                  mrca = value) %>%
    cbind(data.frame(
      y = matDist[allNodeLabels,allNodeLabels] %>% as.vector(),
      x = matRelat[allNodeLabels,allNodeLabels] %>% as.vector(),
      stringsAsFactors = F)) %>%
    dplyr::filter(node1 < node2) %>%
    filter(!is.na(y)); 
  corObj <- cor.test(dfHEreg$x,dfHEreg$y,method="p");
  cor.obs <- corObj$estimate;
  cor.p <- corObj$p.value;
  
  myRet <- data.frame(spearman=cor.obs, spearman.p=cor.p);
  if(nBootstrap > 1) {
    bootstrap.cor <- seq(nBootstrap) %>%
      mclapply(
        mc.cores = nProc,
        function(i) {
          dfHEreg.boot <- dfHEreg %>%
            sample_n(nrow(dfHEreg),replace=T);
          return(cor(dfHEreg.boot$x,dfHEreg.boot$y,method="s"));
        }) %>%
      unlist();
    myRet$cor.se <- sd(bootstrap.cor,na.rm=T) / sqrt(nBootstrap);
  }
  
  if(nPermutation > 0) {
    permut.cor <- seq(nPermutation) %>%
      mclapply(
        mc.cores = nProc,
        function(i) {
          rangePerm <- which(apply(matDist,1,function(x){length(which(is.na(x)))}) < (ncol(matDist) - 1));
          permNodes <- allNodeLabels;
          permNodes[rangePerm] <- sample(permNodes[rangePerm]); ## only permute these nodes
          
          dfHEreg.perm <- allMrca %>%
            reshape2::melt(as.is=T) %>%
            dplyr::rename(node1 = Var1,
                          node2 = Var2,
                          mrca = value) %>%
            cbind(data.frame(
              y = matDist[allNodeLabels,allNodeLabels] %>% as.vector(),
              x = matRelat[permNodes,permNodes] %>% as.vector(),
              stringsAsFactors = F)) %>%
            dplyr::filter(node1 < node2) %>%
            filter(!is.na(y)); 
          return(cor(dfHEreg.perm$x,dfHEreg.perm$y,method="s"));
        }) %>%
      unlist();
    myRet$cor.perm.p <- rank(c(cor.obs,permut.cor))[1] / nPermutation; ## the alternative hypo is observed < permutated
  }
  
  return(myRet);
}

## The actual program starts here
##
allSamples <- c("A1.CBRAD5","G2.CBRAD5","G11.CBRAD5","GS.HESC");
perSample.totalCell <- c();
perSample.totalCell[allSamples] <- c(5000,4000,6450,8000);
perSample.numNode <- c();

## load PCAed transcriptomes (PCA done for all samples together, not individually)
load("~/figs5/all.nodeExp.pca.RData"); 
if(FALSE) {
  dim(all.nodeExp.pca);
  colnames(all.nodeExp.pca) %>% substr(1,6) %>% unique();
}
#rownames(all.nodeExp.pca);
## the "all.nodeExp.pca" includes PCAed transcriptomes (one PC per row) from all cells (one cell per column)

## load cell types
dfCellType <- readRDS("~/fig1/all_cbrad5_GS_hesc_tree_dataframe_modify_new.Rds");
dfCellType <- dfCellType %>% select(-1);
dfCellType <- dfCellType %>%
  mutate(sample = gsub("-",".",sample)) %>%
  mutate(to = gsub("-",".",to)) %>%
  mutate(from = gsub("-",".",from));
allCellType <- dfCellType$celltype %>% unique() %>% sort();
allCellType <- allCellType[allCellType != "inter"];

ntip2ndiv <- readRDS("~/figs5/nTips2nDiv.Rds"); ## the mean distance of two tips by 1000 times simulation for the random tree of size from 4 to 2000 

perSample.capturedCell <- c(1877,692,2395,4690);
names(perSample.capturedCell) <- allSamples
allRes <- list();
for(thisSample in allSamples) {
  thisTreedata <- get(paste0(toupper(thisSample),".treeDF")); # load tree data from nwk(read.tree()%>%as.tibble())
  class(thisTreedata) <- c("tbl_tree","tbl_df","tbl","data.frame");
  thisTreedata$label <- ifelse(thisTreedata$label == "1", thisTreedata$node,thisTreedata$label);
  thisTreedata$label <- ifelse(thisTreedata$label == "", "root",thisTreedata$label);
  thisTree.allTip <- setdiff(thisTreedata$node,thisTreedata$parent);
  thisTree.allTip <- thisTreedata$label[match(thisTree.allTip,thisTreedata$node)];
  thisFeatMat <- all.nodeExp.pca %>% select(starts_with(thisSample));
  colnames(thisFeatMat) <- gsub(paste0(thisSample,"_(.+)"),"\\1",perl=T,colnames(thisFeatMat));
  fracCellSampled <- perSample.capturedCell[thisSample] / perSample.totalCell[thisSample];  # 10x genomics cell capture rate
  
  try.x <- c("mrcaNtip","nSeparatingDivision","2pow.nSeparatingDivision");
  try.y <- c("euclidean","pearson");
  myRes <- expand.grid(x=try.x,y=try.y,stringsAsFactors = F);
  collectRes <- list();
  for(i in seq(nrow(myRes))) {
    cat(paste0("Analyzing sample ", thisSample, " using ", myRes[i,"x"], " as x and ", myRes[i,"y"], " as y.\n"));
    collectRes[[i]] <- heritabilityByHEregression(
      dfTree = thisTreedata,dfFeat = thisFeatMat,
      nBootstrap = 0,nPermutation=1000,nProc = 50,
      relatednessMethod = myRes[i,"x"], distMethod = myRes[i,"y"],
      mrcaNtip2sepDiv = ntip2ndiv) %>%
      mutate(x = myRes$x[i],y=myRes$y[i]);
  }
  allRes[[thisSample]] <- collectRes %>% 
    rbind.fill() %>%
    mutate(type = "transcriptome") %>%
    merge(myRes, by = c("x","y"));
  
  ## internal nodes (cell component)
  dfTip.long <-  
    mclapply(thisTree.allTip,
             mc.cores = 50,
             function(thisTip) {
               tidytree::ancestor(thisTreedata,thisTip) %>%
                 mutate(tip = thisTip) %>%
                 select(tip,node) %>%
                 rbind.fill(data.frame(tip=thisTip,node=thisTip)) %>% ## "node" are the ancestors of the "tip"
                 return();
             }) %>%
    rbind.fill();
  
  ### df listing depth and number of tips / cells (total and each type) under every internal node
  currentNode <- thisTreedata %>% filter(parent == node);
  currentDepth <- 0;
  currentRes <- data.frame(NULL);
  repeat{
    cat(paste0("Checking depth ",currentDepth,"\r"));
    currentRes <- currentNode %>%
      group_split(node) %>%
      mclapply(mc.cores=50,
               function(thisNode){
                 toTip <- thisNode %>%
                   mutate(depth = currentDepth,
                          size = dfTip.long %>% filter(node == thisNode$node[1]) %>% nrow());
                 toCell <- dfTip.long %>% 
                   filter(node == thisNode$node[1]) %>% 
                   merge(dfCellType %>% filter(sample == thisSample) %>% select("label","celltype"),by.x="tip",by.y="label");
                 nCell <- nrow(toCell);
                 toCellType <- sapply(allCellType,function(x){
                   return(length(which(toCell$celltype == x)))
                 }) %>% t() %>% data.frame() %>% 
                   mutate(intern = thisNode$node)
                 toTip %>% 
                   merge(toCellType, by.x="node",by.y="intern") %>%
                   mutate(nCell = nCell) %>%
                   return();
               }) %>%
      rbind.fill() %>%
      rbind(currentRes);
    currentNode <- thisTreedata %>%
      filter(parent != node) %>%
      filter(parent %in% currentNode$node);
    if(nrow(currentNode) == 0) break;
    currentDepth <- currentDepth + 1;
  }
  dfNode2depth_nTip <- currentRes;
  allCellTypeComp <- dfNode2depth_nTip[allCellType] / dfNode2depth_nTip$nCell;
  rownames(allCellTypeComp) <- dfNode2depth_nTip$node;
  matFullDist.celltypeComp <- amap::Dist(allCellTypeComp,method="euclidean",nbproc=50) %>% as.matrix; ## this is not squared
  
  myRes2 <- expand.grid(x=try.x,y=try.y,stringsAsFactors = F);
  collectRes2 <- list();
  for(i in seq(nrow(myRes2))) {
    cat(paste0("Analyzing sample ", thisSample, " using ", myRes[i,"x"], " as x and ", myRes[i,"y"], " as y.\n"));
    collectRes2[[i]] <- heritabilityByHEregression(
      dfTree = thisTreedata,dfFeat = t(allCellTypeComp),
      nBootstrap = 0,nPermutation=1000,nProc = 50,
      relatednessMethod = myRes2[i,"x"], distMethod = myRes2[i,"y"],
      mrcaNtip2sepDiv = ntip2ndiv) %>%
      mutate(x = myRes2$x[i],y=myRes2$y[i]);
  }
  allRes[[thisSample]] <- collectRes2 %>% 
    rbind.fill() %>%
    mutate(type = "cellTypeComp") %>%
    merge(myRes2, by = c("x","y")) %>%
    rbind.fill(allRes[[thisSample]]);
}

dfAllRes <- allSamples %>%
  lapply(function(ss){return(allRes[[ss]] %>% mutate(sample=ss))}) %>%
  rbind.fill();

dfAllRes %>% arrange(x,y);
## ends here

dfAllRes$type<-as.factor(dfAllRes$type)
dfAllRes$sample<-as.factor(dfAllRes$sample)
dfAllRes.euclidean<-dfAllRes[dfAllRes$y=="euclidean",] # euclidean or pearson
dfAllRes.euclidean.nSeparatingDivision<-dfAllRes.euclidean[dfAllRes.euclidean$x=="nSeparatingDivision",]
dfAllRes.euclidean.nSeparatingDivision$shape<-"circle"
dfAllRes.euclidean.nSeparatingDivision$shape[dfAllRes.euclidean.nSeparatingDivision$type=="transcriptome"]<-"triangle"
dfAllRes.euclidean.nSeparatingDivision$shape.type<-"filled"
dfAllRes.euclidean.nSeparatingDivision$shape.type[dfAllRes.euclidean.nSeparatingDivision$cor.perm.p>"0.05"]<-"hollow"

pdf("~/figs5/heritability.pdf")
ggplot(data = dfAllRes.euclidean.nSeparatingDivision, aes(x = sample, y = spearman))+
  geom_point(aes(color = sample, size = 20, shape = type))+
  theme_classic()+
  labs(x="lung differentiation samples", y = "correlation between relatedness and trait difference")
dev.off()

ggplot(dfAllRes.euclidean.nSeparatingDivision, aes(x=sample, y=spearman, shape=factor(shape), fill=factor(shape.type))) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21, 24),
                     name = "shapes!")+
  scale_fill_manual(values=c("black", NA))