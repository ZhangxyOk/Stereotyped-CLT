  library(cli);
  library(GEOquery);
  library(org.Hs.eg.db);
  library(dplyr);
  library(scCustomize)
  library(plyr);
  library(tidyverse);
  library(parallel);
  library(amap);
  library(dplyr);
  library(plyr);
  library(tidyverse);
  library(Seurat);
  library(ggplot2);
  library(ggtree);
  library(ape);
  library(viridis);
  library(scales);
  library(ggpubr);
  library(ggrastr);
  library(ggbreak);
  library(patchwork);
  library(enrichR);
  library(tidyverse);
  library(tidytext);
  library(RColorBrewer);
  library(scCustomize);
  library(rlist);
  library(viridis);
  library(magrittr);
  library(tidyr);
  library(RColorBrewer);
  library(ComplexHeatmap);
  library(reshape);
  library(ggrepel)})
  library(rstatix)
  
# fig 2d
#######function

## function to calculate variance between groups
##   D : distance matrix, i.e., as.matrix(dist)
##   grp : group indicators for each row of the distance matrix
calBetweenGrpVar <- function(D,grp) {
  g = length(unique(grp));
  n = dim(D)[1];
  G = sapply(unique(grp), function(x) as.integer(x == grp)); # Matrix of indicators
  Eij = G %*% t(G);
  TSS=0.5*sum(D^2,na.rm=T)/n; # total SS
  ng=diag(t(G) %*% G);
  WSS=sum(0.5*apply((Eij * D^2)%*%G, 2, sum,na.rm=T)/ng,na.rm=T); ## within group SS
  BSS=TSS-WSS; ## between group SS
  return(BSS/(g-1));
  ## it is important here that the traditional F ratio [ (BSS/(g-1)) / (WSS/(n-g))  ] is inappropriate here
  ## because the within group SS (WSS) traditionally means only measurement error,
  ## whereas the WSS here also includes further differentiation.
  ## in other words, we assume the error level is the same in all the groups
}

## function to detect differentiation given a control dataset
##   dfTree - a tbl_tree object of the tree structure
##   dfFeat - a data.frame of the features (expression, PC, etc)
##   ctrlFeat - a data.frame of the features from a control dataset
##   nPerm - number of permutation for assessing the significant level
##   distMethod - method to estimate distance, use euclidean or pearson
detectDiff <- function(dfTree,dfFeat,ctrlFeat,nPerm,distMethod) {
  if(dfTree %>% filter(parent == node) %>% nrow() != 1) {
    stop("The dfTree should contain one and only one root !");
  } else if( !is(dfTree,"tbl_tree") ) {
    stop("Expecting a tidytree::tbl_tree for the dfTree");
  }
  
  allTipNodes <- setdiff(dfTree$node,dfTree$parent); # find all tips,including onecelltip and multicelltip
  dfTipUnderNode <- lapply(allTipNodes,
                           function(thisTip) {
                             tidytree::ancestor(dfTree,thisTip) %>%
                               mutate(tip = thisTip) %>%
                               select(tip,node) %>%
                               rbind.fill(data.frame(tip=thisTip,node=thisTip)) %>%
                               return();
                           }) %>%
    rbind.fill() %>%
    merge(dfTree %>% select(node,label), by.x="tip", by.y="node")
  
  ## the full distance matrix
  matFullDist <- amap::Dist(t(dfFeat),method=distMethod) %>% as.matrix;
  nRowDist <- nrow(matFullDist);
  matFullDist.ctrl <- amap::Dist(t(ctrlFeat),method=distMethod) %>% as.matrix;
  nRowCtrl <- nrow(matFullDist.ctrl);
  
  ## do it
  dfAllF.perm <-
    lapply(seq(nPerm),
           function(thisPerm){
             matPermDist <- c();
             if(thisPerm == 1) { ## the observation
               matPermDist <- matFullDist;
             } else { ## the permutation
               permIdx <- sample(nRowCtrl,nRowDist,replace=T);
               matPermDist <- matFullDist.ctrl[permIdx,permIdx];
               ## It is important to note that one should use a pre-designated control dataset as the basis of permutation.
               ## This because differentiation samples are continuously differentiating, 
               ## permutation among them not only captures the errors but also the differentiation.
               ## But for our purpose we don't want the differentiation part being controlled off,
               ## so we should use HESC as the control dataset.
             }
             allF.thisPerm <- dfTree %>%
               group_by(parent) %>%
               do({
                 myFocalParent <- .; 
                 dfRelevantTips <- dfTipUnderNode %>% 
                   filter(node %in% myFocalParent$node);
                 relevantIdx <- match(dfRelevantTips$label,rownames(matFullDist)); ## match for index in the original distMat, but use it in the permutated
                 matRelDist <- matPermDist[relevantIdx,relevantIdx];
                 thisF <- calBetweenGrpVar(matRelDist,dfRelevantTips$node);
                 data.frame(f = ifelse(is.nan(thisF),0,thisF));
               });
             allF.thisPerm %>%
               mutate(perm = thisPerm) %>%
               return();
           }) %>% 
    rbind.fill();
  dfRes <- dfAllF.perm %>%
    group_by(parent) %>%
    arrange(perm) %>%
    do({
      myDf <- .;
      myDf$f[is.na(myDf$f)] <- 0;
      p.value = (sum(myDf$f[-1] >= myDf$f[1]) + 1)/(nPerm+1);
      data.frame(p.value = p.value,
                 obsF = myDf$f[1],
                 expF = mean(myDf$f[-1]),
                 sdF = sd(myDf$f[-1]));
    });
  dfTree %>%
    merge(dfRes, by.x="parent",by.y="parent", all.x=T) %>%
    return();
}
#  mdepth(dt-d) calculation based on depth(ds)
# function
subtree.depths <- function(tree){
  tree.tibble <- as_tibble(tree) %>% as.data.frame()
  all.subtrees <- subtrees(tree)
  distMatrix <- dist.nodes(tree)
  # build in function to calculate depth
  getDist <- function(leaf.label, node){
    n2 <- dplyr::filter(tree.tibble, label == leaf.label) %>% `[`("node") %>% unlist()
    dis <- distMatrix[node, n2]
    return(dis)
  }
  # get all subtrees depth
  all.subtrees.depth.df <- lapply(seq(length(all.subtrees)), function(s){
    subtree <- all.subtrees[[s]]
    subtree.root <- subtree$name
    subtree.depth <- Map(getDist, subtree$tip.label, subtree$name) %>% unlist() %>% as.vector() %>% max()
    return(data.frame(subtree.root = subtree.root,
                      depth = subtree.depth))
  }) %>% bind_rows()
  # return results
  return(all.subtrees.depth.df)
}
## data import
#or
tree.infos <- readRDS("~/fig1/all_cbrad5_GS_hesc_tree_dataframe_modify_new.Rds")
tree.internal.info<-tree.infos%>%'['(-grep("add",tree.infos$label),)%>%filter(.,type%in%c("inode","root"))
allsample.cordepth<-list()
for(thisSample in allSamples){
  this.sample.info<-filter(tree.internal.info,sample==gsub("\\.","\\-",thisSample))
  this.sample.cor.depth<-mclapply(seq(nrow(this.sample.info)),mc.cores = 50, function(i){
    this.node<-this.sample.info[i,]
    this.node$depth<-this.node$depth-1
    max.depth<-this.sample.info$height[this.sample.info$type=="root"]-1
    cor.depth<-((max.depth-this.node$height)+1)
    this.node<-this.node%>%mutate(.,cor.depth=cor.depth,max.depth=max.depth)
    return(this.node)})%>%rbind.fill()%>%select(.,label=label,sample=sample,node=to,depth=depth,height=height,max.depth=max.depth,cor.depth=cor.depth)
  allsample.cordepth[[thisSample]]<-this.sample.cor.depth}

allsample.cordepth<-do.call("rbind",allsample.cordepth)
allsample.cordepth$node<-gsub(".*_","",allsample.cordepth$node)
AllNode.cor.depth<-allsample.cordepth%>%select(.,sample=sample,node=node,depth=cor.depth)
AllNode.cor.depth$sample<-gsub("-",".",AllNode.cor.depth$sample)

AllNode.depth.bind<-merge(AllNode.depth,AllNode.cor.depth,c("sample","node"))
AllNode.depth.bind$depth<-rowMeans(select(AllNode.depth.bind,
                                          c(depth.x,depth.y)), na.rm = TRUE)
AllNode.depth.bind<-AllNode.depth.bind[,colnames(AllNode.depth.bind)%in%c("sample","node","depth","size")]
#or
AllNode.depth.bind<-AllNode.depth.bind[,colnames(AllNode.depth.bind)%in%c("sample","node","depth")]
colnames(AllNode.depth.bind)<-c("sample","subtree.root","Td")

# calculate
mdepth.df<-
  mclapply(seq(sample.name),mc.cores=50,function(i){
    tree<-read.tree(tree.path[i])
    inSubtrees <- subtrees(tree)
    subtree.depth<-subtree.depths(tree)
    sample.df<-sample.name[i]
    AllNode.depth.bind.df<-AllNode.depth.bind%>%filter(.,sample==gsub("-",".",sample.df))
    AllNode.depth.bind.df$mdepth<-max(subtree.depth$depth)-AllNode.depth.bind.df$Td
    mdepth.df<-AllNode.depth.bind.df%>%dplyr::select("mdepth","sample","subtree.root")
    subtree.depth.info<-merge(subtree.depth,mdepth.df,by="subtree.root")
      return(subtree.depth.info)
    })%>%bind_rows()
	
#write.csv(mdepth.df,"~/result/fig2/PERMANOVA/AllNode.mdepth.bind.csv")
#saveRDS(AllNode.depth.bind,"~/result/fig2/PERMANOVA/AllNode.depth.bind.Rds")

# fig 2e

## PCA for all samples together
##
all.cbra.gshesc.nodeExp<-readRDS("~/fig2/all.cbra.gshesc.nodeExp.Rds")
pca.tmp <- prcomp(t(all.cbra.gshesc.nodeExp));
keepPc <- length(which(cumsum(summary(pca.tmp)$importance[2,]) < 0.95)) + 1; ## only keeping 95% variance
all.nodeExp.pca <- t(pca.tmp$x[,1:keepPc]) %>% as_tibble(); ## gene expression matrix converted to main PCs'(95% variance) value matrix

cumsum(summary(pca.tmp)$importance[2,1:100])

# top 50 and 100 pca
top50.pc.nodeExp<-all.nodeExp.pca[1:50,]%>%as.data.frame()
top100.pc.nodeExp<-all.nodeExp.pca[1:100,]%>%as.data.frame()
nrow(top50.pc.nodeExp)
ncol(top100.pc.nodeExp)
allPCs<-c("top50","top100")

top.pc.dfAllResult<-mclapply(seq(allPCs),mc.cores=60,function(i){
thisPC<-allPCs[i]
allpc.nodeExp.pca<-get(paste0(thisPC,".pc.nodeExp"))
allSamples <- c("A1.CBRAD5","G2.CBRAD5","G11.CBRAD5","GS.HESC");
tree.allSample <- data.frame(NULL);
ctrl.nodeExp.pca <- allpc.nodeExp.pca %>% dplyr::select(starts_with("GS.HESC"));
allResult <-lapply(seq(allSamples),function(x){
  thisSample<-allSamples[x]
  tree.thisSample <- get(paste0(toupper(thisSample),".treeDF"));
  class(tree.thisSample) <- c("tbl_tree","tbl_df","tbl","data.frame"); ## make it a tree table
  tree.thisSample <- tree.thisSample %>% mutate(sample = thisSample);
  tree.allSample <- tree.allSample %>% rbind.fill(tree.thisSample);
  myExprMat <- allpc.nodeExp.pca %>% dplyr::select(starts_with(thisSample));
  colnames(myExprMat) <- gsub(paste0(thisSample,"_(.*)"),"\\1",perl=T,colnames(myExprMat))
  
  allResult.df<-
    detectDiff(dfTree=tree.thisSample,
               dfFeat=myExprMat,
               ctrlFeat=ctrl.nodeExp.pca,
               nPerm=1000,
               distMethod="pearson"); ## using "euclidean" gives similar results
    return(allResult.df)
})
   names(allResult)<- allSamples             

dfAllResult <- allSamples %>%
  lapply(function(x){
    allResult[[x]] %>% 
      mutate(sample = x) %>%
      return()
  }) %>% rbind.fill()
    return(dfAllResult)
})
names(top.pc.dfAllResult)<-paste0(allPCs, ".dfAllResult")
saveRDS(top.pc.dfAllResult,file = "~/fig2/update_top50_and_100_PC_PERMANOVA.pearson.Rds")
## Figure for CDF of permernova significnat node
# depth group according to E_L_cor_3.0.R mTreedist group
distance.root.df<-readRDS("~/fig2/AllNode.depth.bind.Rds")

allSamples<-c("A1.CBRAD5","G2.CBRAD5","G11.CBRAD5","GS.HESC")

## normalize Td from 0 to 1 and group

distance.root.group<-lapply(seq(allSamples),function(x){
  sample.df<-allSamples[x]
  df.distance<-distance.root.df%>%filter(grepl(sample.df,sample))
  df.distance$depth.group<-df.distance$depth/max(df.distance$depth)
  df<-df.distance
  zone<-seq(from=0,to=1,by=0.1)     
  df.new.1<-lapply(seq((length(zone)-1)),function(x){
    min<-zone[x]
    max<-zone[x+1]
    index<-which(df$depth.group>min&df$depth.group<=max)
    if(length(index)<1){
      df.new<-as.data.frame(matrix(nrow=1,ncol=4))
      colnames(df.new)<-colnames(df)
    }
    else
    {df.new<-df[index,]
    df.new$depth.group<-max}
    return(df.new)
  })%>%bind_rows()%>%na.omit()
  df.new.2<-df%>%filter(.,depth<=0)
  df.new.2$depth.group<-0
  df.new.all<-rbind(df.new.2,df.new.1)
  return(df.new.all)
})%>%bind_rows()
### Addtional code for PC permanova
## Add depth to permanova result of significantly changed internal node ratio of lineage tree for all samples for PC

top50.and.100.pc.dfAllResult<-readRDS("~/fig2/update_top50_and_100_PC_PERMANOVA.pearson.Rds")
top50.pc.dfAllResult<-top50.and.100.pc.dfAllResult[[1]]
top100.pc.dfAllResult<-top50.and.100.pc.dfAllResult[[2]]
this.pc.dfAllResult<-top100.pc.dfAllResult # or top50

PC.dfAllResult.Pernode<-list()
  for(thisSample in allSamples){
    dfAllResult.this<-this.pc.dfAllResult%>%filter(.,sample==thisSample)
    
    ###### OR Unique by internal node
    
    dfAllResult.this<-dfAllResult.this[!duplicated(dfAllResult.this$parent),] #length(which(resDiff.allgene$sample=="A1.CBRAD5"))
    
    dfAllResult.this<-dfAllResult.this[,-c(2,3,4,7,8,9)]
    colnames(dfAllResult.this)[1]<-"node"
    PC.dfAllResult.Pernode[[thisSample]]<-dfAllResult.this}
  
  PC.dfAllResult.Pernode <- allSamples %>%
    lapply(function(x){
      PC.dfAllResult.Pernode[[x]] %>% 
        return()
    }) %>% rbind.fill()

PC.dfAllResult.Pernode.depth<-right_join(PC.dfAllResult.Pernode,distance.root.group,by=c("sample","node"))%>%as.data.frame()  ## AllNode.depth.bind,AllNode.cor.depth,AllNode.depth

Focus.pathway.dfAllResult.Pernode.depth<-list(top50.PC.dfAllResult.Pernode.depth,top100.PC.dfAllResult.Pernode.depth)
names(Focus.pathway.dfAllResult.Pernode.depth)<-c("top50.pc","top100.pc")  

##relative depth by max depth y=significant diff node numer/all nodes number per sample

## run 
  
Focus.pathway.dfAllResult.Pernode.depth<-list(top50.PC.dfAllResult.Pernode.depth,top100.PC.dfAllResult.Pernode.depth)
names(Focus.pathway.dfAllResult.Pernode.depth)<-c("top50.pc","top100.pc")
  
  all.trans.CDF.prepare<-lapply(seq(Focus.pathway.dfAllResult.Pernode.depth), function(i){
    
    dfAlldepthResult<-Focus.pathway.dfAllResult.Pernode.depth[[i]]%>%as.data.frame()%>%
      dplyr::select(sample,node,p.value,depth.group)
    
    df.relativeDiffDepth<-dfAlldepthResult%>%dplyr::select(sample,depth.group,p.value)
    colnames(df.relativeDiffDepth)<-c("sample","relD","p.value")
    
    df.relativeDiffDepth.ratio <- 
      df.relativeDiffDepth %>% dplyr::group_by(sample) %>% do({
        df <- .
        all.df <- length(df$relD)
        re.cdf <- lapply(unique(df$relD), function(dp){
          dp.df <- filter(df, relD <= dp, p.value<0.05) # p.value<0.05 significant diff node number
          diff.ratio <- length(dp.df$relD) / all.df
          #diff.ratio <- sum(dp.df$p.value < 0.05) / nrow(df)
          sub.re.df <- data.frame(sample = df$sample[1],
                                  relD = dp,
                                  diff.ratio = diff.ratio)
        }) %>% bind_rows()
        re.cdf
      })
    pc.name<-names(Focus.pathway.dfAllResult.Pernode.depth)[[i]]
    df.relativeDiffDepth.ratio<-df.relativeDiffDepth.ratio%>%mutate(PC=pc.name)%>%as.data.frame()
    return(df.relativeDiffDepth.ratio)})%>%rbind.fill()

	all.trans.CDF.prepare %>% 
  ggplot(aes(x=relD, y=diff.ratio, color=sample)) +
  geom_point(shape=1,size=8) +
  geom_line(aes(group=sample, color=sample),size=1) +
  scale_color_manual(values=c("#F46D43","#4DAF4A","#74ADD1","#BDBDBD"))+
  theme_bw() + 
  theme_classic(base_size = 14)+
  labs(color="Sample", x="Relative Depth", y = "Differentail Fraction") +
  facet_wrap(~ PC, ncol = 2)+
  theme(strip.text.x = element_text(size=8),
        strip.background = element_rect(fill="white"),
        axis.title = element_text(size = 10))  
# fraction of cbrad5 samples 
df1<-all.trans.CDF.prepare%>%filter(sample=="A1.CBRAD5"&PC=="top100.pc")%>%dplyr::select(diff.ratio)%>%max()
df2<-all.trans.CDF.prepare%>%filter(sample=="G11.CBRAD5"&PC=="top100.pc")%>%dplyr::select(diff.ratio)%>%max()
df3<-all.trans.CDF.prepare%>%filter(sample=="G2.CBRAD5"&PC=="top100.pc")%>%dplyr::select(diff.ratio)%>%max()
mean(c(df1,df2,df3))

# fig 2a

top<- GSE83310_SeuratObj.markers %>% group_by(cluster) %>% top_n(n = length(avg_log2FC)*0.1, wt = avg_log2FC) ## choose top10% GEGs ranked by avg_log2FC and 
top.marker<-top%>%filter(.,cluster%in%c("neural GFP+","day 0 undifferentiated iPSC17","day 3 definitive endoderm",
                                        "day 6 anterior foregut endoderm","day 15 Nkx2-1-GFP-","day 15 Nkx2-1-GFP+"))%>%
  dplyr::select(gene)
  
DoHeatmap(GSE83310_SeuratObj, features = as.character(top.marker$gene),hjust=0, angle = 90,slot = "scale.data", draw.lines = FALSE, group.by= "ident",raster=FALSE)+
  scale_fill_gradientn(colors = rev(colorRampPalette(c("#F46D43", "#FFFFFF", "#3288BD"))(256))) +
  theme(text = element_text(size = 3))#+ NoLegend()

# fig 2b

de.0 <- top%>%subset(.,cluster=="day 0 undifferentiated iPSC17")%>%dplyr::select(gene)
de.1 <- top%>%subset(.,cluster=="day 3 definitive endoderm")%>%dplyr::select(gene)
de.2 <- top%>%subset(.,cluster=="day 6 anterior foregut endoderm")%>%dplyr::select(gene)
de.3 <- top%>%subset(.,cluster=="day 15 Nkx2-1-GFP+")%>%dplyr::select(gene)
de.4 <- top%>%subset(.,cluster=="day 15 Nkx2-1-GFP-")%>%dplyr::select(gene)
de.5 <- top%>%subset(.,cluster=="neural GFP+")%>%dplyr::select(gene)

dbs <- listEnrichrDbs()

dbs <- dbs[order(dbs$libraryName),]

websiteLive <- getOption("enrichR.live")

listEnrichrSites()

setEnrichrSite("Enrichr")

dbs_pw <- c("KEGG_2021_Human", "WikiPathways_2021_Human", "BioPlanet_2019","MSigDB_Hallmark_2020",  
            "Elsevier_Pathway_Collection","Reactome_2022","NCI-Nature_2016")
#de<-top.marker
# or
de<-de.5
enriched_pw<-enrichr(de$gene, dbs_pw)
enriched_pw.df<-do.call("rbind",enriched_pw)%>%filter(Adjusted.P.value<0.05)%>% #0.05 for top0.1 or 0.1 for top0.05
  arrange(Combined.Score,decreasing=TRUE)%>%filter(!grepl("mm",Term))%>%filter(!grepl("mouse",Term))%>%
  filter(!grepl("MOUSE",Term))%>%dplyr::select(Term,Adjusted.P.value,P.value,Overlap,Odds.Ratio,Genes,Combined.Score)

#enriched_pw.allDEGs<-enriched_pw.df
# or

#enriched_pw.day0<-enriched_pw.df%>%mutate(stage="day 0 undifferentiated iPSC17")%>%
  #enriched_pw.day3<-enriched_pw.df%>%mutate(stage="day 3 definitive endoderm")%>%
  #enriched_pw.day6<-enriched_pw.df%>%mutate(stage="day 6 anterior foregut endoderm")%>%
  #enriched_pw.day15_nkx2_1_pos<-enriched_pw.df%>%mutate(stage="day 15 Nkx2-1-GFP+")%>%
  #enriched_pw.day15_nkx2_1_neg<-enriched_pw.df%>%mutate(stage="day 15 Nkx2-1-GFP-")%>%
  enriched_pw.neural_nkx2_1_pos<-enriched_pw.df%>%mutate(stage="neural GFP+")%>%
  #arrange(Combined.Score,decreasing=TRUE)%>%
distinct(Genes, .keep_all = TRUE)

enriched_pw.allDEGs.top0.1<-rbind(enriched_pw.neural_nkx2_1_pos,enriched_pw.day0,enriched_pw.day3,
                             enriched_pw.day6,enriched_pw.day15_nkx2_1_neg,
                             enriched_pw.day15_nkx2_1_pos)%>%distinct(Term, .keep_all = TRUE)
							 
 # remove similar geneset pathways
stage<-unique(enriched_pw.allDEGs.top0.1$stage)

enriched_pw.allDEGs<-lapply(seq(stage), function(i){
  df.stage<-stage[i]
  df<-enriched_pw.allDEGs.top0.1%>%filter(stage==df.stage)
  df<-df%>%arrange(Adjusted.P.value,decreasing=FALSE)
  return(df)
})%>%bind_rows()

pw.gene.all<-lapply(seq(nrow(enriched_pw.allDEGs)),function(i){ # enriched_pw 
  df<-enriched_pw.allDEGs
  table<-df[i,]
  gene<-table$Genes%>%as.data.frame()
  colnames(gene)<-"gene"
  gene<-gene%>% 
    mutate(gene = strsplit(as.character(gene), ";")) %>% 
    unnest(gene)
  return(gene)
})
names(pw.gene.all)<-enriched_pw.allDEGs$Term
	     
pw.gene<-pw.gene.all
names(pw.gene)<-paste0(enriched_pw.allDEGs$stag,"_",names(pw.gene))

colnames(sce@meta.data)
sce@meta.data<-sce@meta.data[,-c(18:202)]

sce_mutate_score<-mclapply(seq(pw.gene),mc.cores = 50,function(x){
  dataset<-pw.gene[[x]]
  gene<-dataset$gene
  score.name<-names(pw.gene)[x]%>%as.character()
  sce.new <- try(Seurat::AddModuleScore(object = sce,features = list(gene),name = paste(score.name,"1",sep="")), silent = T)
  if(class(sce.new) == "try-error"){score=NA}else{score<-sce.new@meta.data[,length(sce.new@meta.data)]%>%as.numeric()}
  return(score)})%>%bind_cols()

colnames(sce_mutate_score)<-names(pw.gene)
rownames(sce_mutate_score)<-rownames(sce@meta.data)

na<-apply(sce_mutate_score, 2, function (x) all(is.na(x)))
na.pathway<-colnames(sce_mutate_score)[which(na)]
#[1] "neural GFP+_RNA Polymerase I Promoter Opening R-HSA-73728"              
#[2] "neural GFP+_Packaging Of Telomere Ends R-HSA-171306"                    
#[3] "neural GFP+_Systemic lupus erythematosus"                               
#[4] "neural GFP+_DNA Damage/Telomere Stress Induced Senescence R-HSA-2559586"
#[5] "neural GFP+_Nonhomologous End-Joining (NHEJ) R-HSA-5693571"             
#[6] "neural GFP+_HATs Acetylate Histones R-HSA-3214847" 
sce_mutate_score.new<-sce_mutate_score[,-which(na)]
enriched_pw.filter.top.0.1<-enriched_pw.allDEGs[-which(enriched_pw.allDEGs$Term%in%gsub(".*_","",na.pathway)),]

Features<-names(pw.gene)[-which(na)]%>%as.character()

Feature.dotplot<-DotPlot(sce,features = Features,cols = "RdYlBu",dot.scale=6)+RotatedAxis()+ # scaled the averaged expression of signal by dotplot
  ggplot2::coord_flip()+theme(axis.text=element_text(size=13))

Feature.avg.exp<-Feature.dotplot[[1]]

cluster_id<-sce@active.ident%>%as.data.frame()
colnames(cluster_id)<-"cluster"

celltype<-unique(Feature.avg.exp$id)
Feature.avg.exp.matrix<-lapply(seq(celltype), function(i){ # scaled the averaged expression of signal by dotplot
  pathway<-unique(Feature.avg.exp$features.plot)
  Feature.avg.exp.sub<-filter(Feature.avg.exp,id==celltype[i])
  df<-Feature.avg.exp.sub$avg.exp.scaled[match(pathway,Feature.avg.exp.sub$features.plot)]%>%as.numeric()%>%as.data.frame()%>%t()
  rownames(df)<-celltype[i]
  colnames(df)<-pathway
  df<-df%>%as.data.frame()
  return(df)
})%>%bind_rows()

is.na(Feature.avg.exp.matrix)

# cluster sequence
index<-match(c("C10","C9","C8","C7","C6","C5","C4","C3","C2","C1","R2","R1"),rownames(Feature.avg.exp.matrix))
Feature.avg.exp.matrix.new<-Feature.avg.exp.matrix[index,]
colnames(Feature.avg.exp.matrix.new)<-gsub("R-HSA.*","",colnames(Feature.avg.exp.matrix.new))
colnames(Feature.avg.exp.matrix.new)<-gsub("Homo sapiens.*","",colnames(Feature.avg.exp.matrix.new))
colnames(Feature.avg.exp.matrix.new)<-gsub("top0.1_","",colnames(Feature.avg.exp.matrix.new))

 Heatmap(as.matrix(Feature.avg.exp.matrix.new),
            column_title = "lung progenitor differetiation related pathway", row_title = "cell cluster",
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize=6),# Text size for row names
            col<- rev(colorRampPalette(c("#F46D43", "#FFFFFF", "#3288BD"))(256)),
            cluster_rows = FALSE,
            cluster_columns = FALSE
)

# fig 2f
## permanova for all pathway geneset

all.pathway.geneset<-pw.gene.all

all.pathway.nodeExp<-mclapply(seq(all.pathway.geneset),mc.cores = 60,function(i){
  this.pathway.geneset<-all.pathway.geneset[[i]]
  assign.name<-names(all.pathway.geneset)[i]
  index.all<-which(rownames(all.cbra.gshesc.nodeExp)%in%this.pathway.geneset$gene)
  this.pathway.nodeExp<-all.cbra.gshesc.nodeExp[index.all,]
})
names(all.pathway.nodeExp)<-names(all.pathway.geneset)

all.pathway.dfAllResult<-mclapply(seq(all.pathway.nodeExp),mc.cores=60,function(i){
  thisPathway<-names(all.pathway.nodeExp)[i]
  this.pathway.nodeExp<-all.pathway.nodeExp[[i]]
  allSamples <- c("A1.CBRAD5","G2.CBRAD5","G11.CBRAD5","GS.HESC");
  tree.allSample <- data.frame(NULL);
  ctrl.nodeExp.pca <- this.pathway.nodeExp %>%as.data.frame()%>% dplyr::select(starts_with("GS.HESC"));
  allResult <-lapply(seq(allSamples),function(x){
    thisSample<-allSamples[x]
    tree.thisSample <- get(paste0(toupper(thisSample),".treeDF"));
    class(tree.thisSample) <- c("tbl_tree","tbl_df","tbl","data.frame"); ## make it a tree table
    tree.thisSample <- tree.thisSample %>% mutate(sample = thisSample);
    tree.allSample <- tree.allSample %>% rbind.fill(tree.thisSample);
    myExprMat <- this.pathway.nodeExp %>%as.data.frame()%>%dplyr::select(starts_with(thisSample));
    colnames(myExprMat) <- gsub(paste0(thisSample,"_(.*)"),"\\1",perl=T,colnames(myExprMat))
    
    allResult.df<-
      detectDiff(dfTree=tree.thisSample,
                 dfFeat=myExprMat,
                 ctrlFeat=ctrl.nodeExp.pca,
                 nPerm=1000,
                 distMethod="pearson"); ## using "euclidean" gives similar results
    return(allResult.df)
    })
  names(allResult)<- allSamples
  
  dfAllResult <- allSamples %>%
    lapply(function(x){
      allResult[[x]] %>% 
        mutate(sample = x) %>%
        return()
    }) %>% rbind.fill()
  return(dfAllResult)
})

names(all.pathway.dfAllResult)<-names(all.pathway.nodeExp)

all.pathway.dfAllResult.df<-all.pathway.dfAllResult
# the string gsub are necessary
#names(all.pathway.dfAllResult.df)<-gsub("R-HSA","",names(all.pathway.dfAllResult.df))
#names(all.pathway.dfAllResult.df)<-gsub("Homo sapiens.*","",names(all.pathway.dfAllResult.df))
#names(all.pathway.dfAllResult.df)<-gsub("2-1-GFP-","2_1_GFP_neg",names(all.pathway.dfAllResult.df))
#names(all.pathway.dfAllResult.df)<-gsub("2-1-GFP","2_1_GFP",names(all.pathway.dfAllResult.df))
#names(all.pathway.dfAllResult.df)<-gsub("Beta-","Beta_",names(all.pathway.dfAllResult.df))
#names(all.pathway.dfAllResult.df)<-gsub("Pre-","pre_",names(all.pathway.dfAllResult.df))
#names(all.pathway.dfAllResult.df)<-gsub("HIF-","HIF_",names(all.pathway.dfAllResult.df))
#names(all.pathway.dfAllResult.df)<-gsub("-.*","",names(all.pathway.dfAllResult.df))
#names(all.pathway.dfAllResult.df)<-gsub("2_1_GFP_neg","2_1_GFP-",names(all.pathway.dfAllResult.df))
#names(all.pathway.dfAllResult.df)<-gsub("_1","",names(all.pathway.dfAllResult.df))

allSamples <- c("A1.CBRAD5","G2.CBRAD5","G11.CBRAD5","GS.HESC");

sig.diff.node.ratio.all<-lapply(seq(all.pathway.dfAllResult.pearson), function(i){
  this.pathway.dfAllResult<-all.pathway.dfAllResult.pearson[[i]]
  pathway.euclidean.dfAllResult.new<-list()
  for(thisSample in allSamples){
    dfAllResult.this<-this.pathway.dfAllResult%>%filter(.,sample==thisSample)
    dfAllResult.this<-dfAllResult.this[!duplicated(dfAllResult.this$parent),] #length(which(resDiff.allgene$sample=="A1.CBRAD5"))
    dfAllResult.this<-dfAllResult.this[,-c(2,3,4)]
    pathway.euclidean.dfAllResult.new[[thisSample]]<-dfAllResult.this}
  
  pathway.euclidean.dfAllResult.new <- allSamples %>%
    lapply(function(x){
      pathway.euclidean.dfAllResult.new[[x]] %>% 
        return()
    }) %>% 
    rbind.fill()
  
  ratio<-pathway.euclidean.dfAllResult.new %>%
    group_by(sample) %>%
    filter(!is.na(p.value)) %>%
    dplyr::summarise(num.sig = sum(p.value < 0.05),
                     frac.sig = mean(p.value < 0.05))%>%t()
  colnames<-ratio[1,]
  ratio<-ratio[3,]%>%as.numeric()
  ratio$pathway<-names(all.pathway.dfAllResult.pearson)[[i]]
  ratio<-as.data.frame(ratio)
  colnames(ratio)<-c(colnames,"pathway")
  return(ratio)
})%>% rbind.fill()

# filter.by.mean.fraction.100percentile.hesc.sample for 185 pathway geneset	  
hesc.ratio.100percentile<-max(sig.diff.node.ratio.all$GS.HESC)

sig.diff.node.ratio<-lapply(seq(all.pathway.dfAllResult.df), function(i){
  this.pathway.dfAllResult<-all.pathway.dfAllResult.df[[i]]
  pathway.euclidean.dfAllResult.new<-list()
  for(thisSample in allSamples){
    dfAllResult.this<-this.pathway.dfAllResult%>%filter(.,sample==thisSample)
    dfAllResult.this<-dfAllResult.this[!duplicated(dfAllResult.this$parent),] #length(which(resDiff.allgene$sample=="A1.CBRAD5"))
    dfAllResult.this<-dfAllResult.this[,-c(2,3,4)]
    pathway.euclidean.dfAllResult.new[[thisSample]]<-dfAllResult.this}
  
  pathway.euclidean.dfAllResult.new <- allSamples %>%
    lapply(function(x){
      pathway.euclidean.dfAllResult.new[[x]] %>% 
        return()
    }) %>% 
    rbind.fill()
  
  ratio<-pathway.euclidean.dfAllResult.new %>%
    group_by(sample) %>%
    filter(!is.na(p.value)) %>%
    dplyr::summarise(num.sig = sum(p.value < 0.05),
                     frac.sig = mean(p.value < 0.05))%>%t()
  colnames<-ratio[1,]
  ratio<-ratio[3,]%>%as.numeric()
  ratio$pathway<-names(all.pathway.dfAllResult.df)[[i]]
  ratio<-as.data.frame(ratio)
  colnames(ratio)<-c(colnames,"pathway")
  return(ratio)
})%>% rbind.fill()
sig.diff.node.ratio<-mutate(sig.diff.node.ratio, mean_col = rowMeans(dplyr::select(sig.diff.node.ratio, 
                                              c(A1.CBRAD5,G11.CBRAD5,G2.CBRAD5)), na.rm = TRUE))


sig.diff.pathway.100percentile.hesc.pear<-filter(sig.diff.node.ratio,mean_col>=hesc.ratio.100percentile)%>%'['(,5)

index<-which(sig.diff.node.ratio$mean_col>=hesc.ratio.100percentile)

Focus.pathway.dfAllResult<-all.pathway.dfAllResult.df[index]

Focus.pathway.dfAllResult.Pernode<-lapply(seq(Focus.pathway.dfAllResult), function(i){
  this.pathway.dfAllResult<-Focus.pathway.dfAllResult[[i]]
  Focus.pathway.dfAllResult.Pernode<-list()
  for(thisSample in allSamples){
    dfAllResult.this<-this.pathway.dfAllResult%>%filter(.,sample==thisSample)
    dfAllResult.this<-dfAllResult.this[!duplicated(dfAllResult.this$parent),] ## unique for internal node
    dfAllResult.this<-dfAllResult.this[,-c(2,3,4,7,8,9)]
    colnames(dfAllResult.this)[1]<-"node"
    Focus.pathway.dfAllResult.Pernode[[thisSample]]<-dfAllResult.this}
  
  Focus.pathway.dfAllResult.Pernode <- allSamples %>%
    lapply(function(x){
      Focus.pathway.dfAllResult.Pernode[[x]] %>% 
        return()
    }) %>% 
    rbind.fill()
})
names(Focus.pathway.dfAllResult.Pernode)<-names(Focus.pathway.dfAllResult)

distance.root.df<-readRDS("~/fig2/AllNode.depth.bind.Rds")

## normalize Td from 0 to 1 and group

distance.root.group<-lapply(seq(allSamples),function(x){
  sample.df<-allSamples[x]
  df.distance<-distance.root.df%>%filter(grepl(sample.df,sample))
  df.distance$depth.group<-df.distance$depth/max(df.distance$depth)
  df<-df.distance
  zone<-seq(from=0,to=1,by=0.1)     
  df.new.1<-lapply(seq((length(zone)-1)),function(x){
    min<-zone[x]
    max<-zone[x+1]
    index<-which(df$depth.group>min&df$depth.group<=max)
    if(length(index)<1){
      df.new<-as.data.frame(matrix(nrow=1,ncol=4))
      colnames(df.new)<-colnames(df)
    }
    else
    {df.new<-df[index,]
    df.new$depth.group<-max}
    return(df.new)
  })%>%bind_rows()%>%na.omit()
  df.new.2<-df%>%filter(.,depth<=0)
  df.new.2$depth.group<-0
  df.new.all<-rbind(df.new.2,df.new.1)
  return(df.new.all)
})%>%bind_rows()

### Add depth to permanova result of significantly changed internal node ratio of lineage tree for all samples for pathway

Focus.pathway.dfAllResult.Pernode.depth<-list()
Focus.pathway.dfAllResult.Pernode.depth<-lapply(seq(Focus.pathway.dfAllResult.Pernode),function(i){
  this.result.Pernode<-Focus.pathway.dfAllResult.Pernode[[i]]%>%as.data.frame()
  this.result.Pernode.depth<-right_join(this.result.Pernode,distance.root.group,by=c("sample","node"))%>%as.data.frame() 
  Focus.pathway.dfAllResult.Pernode.depth<-this.result.Pernode.depth%>%mutate(pathway=names(Focus.pathway.dfAllResult.Pernode)[i])
  return(Focus.pathway.dfAllResult.Pernode.depth)})
names(Focus.pathway.dfAllResult.Pernode.depth)<-names(Focus.pathway.dfAllResult.Pernode)

## CDF plot
all.pathway.CDF.prepare<-lapply(seq(Focus.pathway.dfAllResult.Pernode.depth), function(i){
  
  dfAlldepthResult<-Focus.pathway.dfAllResult.Pernode.depth[[i]]%>%as.data.frame()%>%
    dplyr::select(sample,node,p.value,depth.group)
  
  df.relativeDiffDepth<-dfAlldepthResult%>%dplyr::select(sample,depth.group,p.value)
  colnames(df.relativeDiffDepth)<-c("sample","relD","p.value")
  
  df.relativeDiffDepth.ratio <- 
    df.relativeDiffDepth %>% dplyr::group_by(sample) %>% do({
      df <- .
      all.df <- length(df$relD)
      re.cdf <- lapply(unique(df$relD), function(dp){
        dp.df <- filter(df, relD <= dp, p.value<0.05) # p.value<0.05 significant diff node number
        diff.ratio <- length(dp.df$relD) / all.df
        #diff.ratio <- sum(dp.df$p.value < 0.05) / nrow(df)
        sub.re.df <- data.frame(sample = df$sample[1],
                                relD = dp,
                                diff.ratio = diff.ratio)
      }) %>% bind_rows()
      re.cdf
    })
  pathway.name<-names(Focus.pathway.dfAllResult.Pernode.depth)[[i]]
  df.relativeDiffDepth.ratio<-df.relativeDiffDepth.ratio%>%mutate(pathway=pathway.name)%>%as.data.frame()
  return(df.relativeDiffDepth.ratio)})%>%rbind.fill()


all.pathway.CDF.prepare %>% 
  ggplot(aes(x=relD, y=diff.ratio, color=sample)) +
  geom_point() +
  geom_line(aes(group=sample, color=sample)) +
  #theme_bw() + 
  theme_classic(base_size = 14)+
  labs(color="Sample", x="Relative Depth", y = "Differentail Fraction") +
  facet_wrap(~ pathway, ncol = 6)+
  theme(strip.text.x = element_text(size=8),
        strip.background = element_rect(fill="white"),
        axis.title = element_text(size = 10))
		
#all.pathway.CDF.prepare%>%filter(grepl("top0.1_",pathway))%>%
  all.pathway.CDF.prepare%>%
  ggplot(aes(x=relD, y=diff.ratio, color=sample)) +
  geom_point(shape=1,size=8) +
  geom_line(aes(group=sample, color=sample),size=2) +
  geom_hline(aes(yintercept = 0.05), 
             linetype="dashed",
             size = 0.5, 
             color = "#585858")+
  scale_color_manual(values=c("#F46D43","#4DAF4A","#74ADD1","#BDBDBD"))+
  theme_bw() + 
  theme_classic(base_size = 14)+
  labs(color="Sample", x="Normalized Depth", y = "Cumulative fraction") +
  facet_wrap(~ pathway, ncol = 10,nrow=12)+
  theme(strip.text.x = element_text(size=8),
        strip.background = element_rect(fill="white"),
        axis.title = element_text(size = 10))
		
## fig 2g
## tracing divergence event for signals which ratio >=100% percentile hesc ratios

Focus.pathway.dfAllResult.Pernode.depth.rbind <- do.call("rbind", Focus.pathway.dfAllResult.Pernode.depth)

Focus.pathway.dfAllResult.Pernode.depth.rbind<-Focus.pathway.dfAllResult.Pernode.depth.rbind%>%filter(.,p.value<0.05)# filter significant permanova

# ##  to find the most early ancestor node of 
alldiffSamples<-allSamples[1:3]
dfTip.long$sample<-gsub("_.*","",dfTip.long$intern)
dfTip.long$sample<-gsub("-",".",dfTip.long$sample) 
dfTip.long$node<-gsub(".*_","",dfTip.long$intern)%>%as.numeric()

sig.depth.df<-Focus.pathway.dfAllResult.Pernode.depth.rbind%>%dplyr::select("node","sample","p.value","depth.group","pathway")
colnames(sig.depth.df)<-c("node","sample","p.value","depth","pathway")
allSigPathway <- sig.depth.df$pathway %>% unique()

## find initial time of signal divergence

## GS HESC plot as control use 10 percentile instead of mean as initial 

## find initial time of signal divergence
AllNode.depth.sub<-distance.root.group%>%dplyr::select("sample","node","depth.group")
colnames(AllNode.depth.sub)<-c("sample","node","depth")
# new code 
# key: divergnece initial
beginApp<-mclapply(seq(allSigPathway),mc.cores=60,function(i){
  this.signal.Result<-sig.depth.df%>%filter(pathway==allSigPathway[i])
  
  this.signal.leaf.ancestor.result<-lapply(seq(allSamples),function(f){
    dfTip.long.sub<-filter(dfTip.long,sample==allSamples[f])%>%filter(.,!grepl("add",intern))
    this.sample.Result<-this.signal.Result%>% filter(sample==allSamples[f])
    
    leave.ancestor<-lapply(seq(unique(dfTip.long.sub$tip)),function(x){
      this.leaf<-dfTip.long.sub[dfTip.long.sub$tip==unique(dfTip.long.sub$tip)[x],]%>%arrange(node)
      #this.leaf<-this.leaf[-which(this.leaf$node%in%root),]
      this.leaf.this.signal.Result<-this.sample.Result%>%filter(.,node%in%this.leaf$node)
      condition<-this.leaf$node%in%this.leaf.this.signal.Result$node
      if (length(which(condition=="FALSE"))>=1){ # negative >1
        if(length(which(condition=="FALSE"))<length(condition)){ # 1. discontinuous (negative + positive)
          if(condition[1]=="FALSE"){   # 1.1 discontinuous: root negative
            ancestor<-this.leaf$node[which(condition=="TRUE")]%>%as.numeric()%>%min()  ## nearest ancestor significant 
            this.leaf.ancestor.Result<-data.frame(
              node=ancestor,
              sample=allSamples[f],
              pathway=allSigPathway[i],
              depth=(filter(AllNode.depth.sub,sample==allSamples[f]&node==ancestor)%>%
                       dplyr::select(.,depth)%>%'['(1,1)))
          }else{this.leaf.ancestor.Result<-this.leaf.this.signal.Result%>%  # 1.2 discontinuous: root positive
            arrange(node)%>%'['(1,)%>%dplyr::select(depth,sample,node,pathway) # root
          }
        }else{ancestor<-this.leaf$node[which(condition=="FALSE")]%>%as.numeric()%>%max() # 2. all negative  
        this.leaf.ancestor.Result<-data.frame(  
          node=ancestor,
          sample=allSamples[f],
          pathway=allSigPathway[i],
          depth=(AllNode.depth.sub%>%filter(sample==allSamples[f]&node==ancestor)%>%  # infer the most late appearred time
                   dplyr::select(.,depth)%>%'['(1,1)+0.1))
        }}else{this.leaf.ancestor.Result<-this.leaf.this.signal.Result%>%       # 3. continuous (all positive)  
          arrange(node)%>%'['(1,)%>%dplyr::select(depth,sample,node,pathway)} # root
      this.leaf.ancestor.Result$node<-as.numeric(this.leaf.ancestor.Result$node)
      return(this.leaf.ancestor.Result)})%>%bind_rows()
    return(leave.ancestor)})%>%bind_rows()
  return(this.signal.leaf.ancestor.result)})

# Key: divergence end
beginDisapp<-mclapply(seq(allSigPathway),mc.cores=50,function(i){
  this.signal.Result<-sig.depth.df%>%filter(pathway==allSigPathway[i])
  
  this.signal.leaf.ancestor.result<-lapply(seq(allSamples),function(f){
    dfTip.long.sub<-filter(dfTip.long,sample==allSamples[f])%>%filter(.,!grepl("add",intern))
    this.sample.Result<-this.signal.Result%>% filter(sample==allSamples[f])
    
    leave.ancestor<-lapply(seq(unique(dfTip.long.sub$tip)),function(x){
      this.leaf<-dfTip.long.sub[dfTip.long.sub$tip==unique(dfTip.long.sub$tip)[x],]%>%arrange(node)
      #this.leaf<-this.leaf[-which(this.leaf$node%in%root),]
      this.leaf.this.signal.Result<-this.sample.Result%>%filter(.,node%in%this.leaf$node)
      condition<-this.leaf$node%in%this.leaf.this.signal.Result$node
      if (length(which(condition=="FALSE"))>=1){ # negative >1
        if(length(which(condition=="FALSE"))<length(condition)){ # 1. discontinuous (negative + positive)
          ancestor<-this.leaf$node[which(condition=="TRUE")]%>%as.numeric()%>%max() ## most late ancestor significant diverged in this pathway
          this.leaf.ancestor.Result<-data.frame(
            node=ancestor,
            sample=allSamples[f],
            pathway=allSigPathway[i],
            depth=(AllNode.depth.sub%>%filter(sample==allSamples[f]&node==ancestor)%>%
                     dplyr::select(.,depth)%>%'['(1,1)))
        }else{ancestor<-this.leaf$node[which(condition=="FALSE")]%>%as.numeric()%>%max() # 2. all negative
          this.leaf.ancestor.Result<-data.frame(  
          node=ancestor,
          sample=allSamples[f],
          pathway=allSigPathway[i],
          #depth=min(this.sample.Result$depth # infer the mean according events happened：most early disappear
          depth=(AllNode.depth.sub%>%filter(sample==allSamples[f]&node==ancestor)%>% # infer the most late appearred time
             dplyr::select(.,depth)%>%'['(1,1))+0.1)  
        }}else{ancestor<-this.leaf$node[which(condition=="TRUE")]%>%as.numeric()%>%max() # 3. continuous (all positive)
          this.leaf.ancestor.Result<-data.frame( 
          node=ancestor,
          sample=allSamples[f],
          pathway=allSigPathway[i],
          depth=(AllNode.depth.sub%>%filter(sample==allSamples[f]&node==ancestor)%>% # infer the most late appearred time
                   dplyr::select(depth)%>%'['(1,1))+0.1)} # most late ancestor significant divereged +0.1
      this.leaf.ancestor.Result$node<-as.numeric(this.leaf.ancestor.Result$node)
      return(this.leaf.ancestor.Result)})%>%
      bind_rows()
    return(leave.ancestor)})%>%
    bind_rows()
  return(this.signal.leaf.ancestor.result)})

fate.signal<-beginApp
# or
fate.signal<-beginDisapp

names(fate.signal)<-allSigPathway
fate.signal.data<-do.call("rbind",fate.signal)
fate.signal.data$depth<-gsub("Inf","1",fate.signal.data$depth)%>%as.numeric()
depth.initial<-fate.signal.data%>%na.omit()

# depth.initial.euclidean<-depth.initial
# depth.initial.pearson<-depth.initial

mean.depth.initial<-lapply(seq(allSigPathway),function(i){
  sub <- depth.initial%>% filter(pathway == allSigPathway[i])
  mean.initial<- sub%>% group_by(sample) %>% dplyr::summarise(obsD = mean(depth,na.rm=TRUE))%>%
    mutate(pathway=allSigPathway[i])
  return(mean.initial)
})%>%bind_rows()

plot.depth.initial<-mean.depth.initial%>%reshape2::dcast(pathway ~ sample, value.var = "obsD")%>%
  mutate(mean_col = rowMeans(dplyr::select(.,c(A1.CBRAD5,G11.CBRAD5,G2.CBRAD5)),na.rm = TRUE))


df<-plot.depth.initial #%>%filter(grepl("top0.1_",pathway))
df$pathway<-gsub("top0.1_","",df$pathway)

plot.app.pearson<-df
plot.Disapp.pearson<-df

#pathway.heatmap<-names(all.pathway.dfAllResult.df)
df<-plot.app.pearson#%>%filter(grepl("day 15",pathway))
df<-plot.Disapp.pearson#%>%filter(grepl("day 15",pathway))
df$pathway<-gsub("Nkx2_","Nkx2-1-",df$pathway)
# df$pathway<-gsub("Nkx2_1_1","Nkx2_1",df$pathway)
#df<-df[match(colnames(Feature.avg.exp.matrix.new),df$pathway),]%>%na.omit()
#df<-df[match(names(all.pathway.dfAllResult.df),df$pathway),]%>%na.omit()
df$pathwayNum<-seq(df$pathway)
df<-reshape::melt(as.data.frame(df),id=c("pathwayNum","pathway"))
df.2<-df%>%filter(.,variable=="mean_col"|variable=="GS.HESC") 
df.2$type<-"divergence disappear" # "divergence appear" or "divergence disappear"
df.initiation<-df.2
df.end<-df.2
df<-rbind(df.initiation,df.end)
df$state<-gsub(" _.*","",df$pathway)%>%as.character()

index<-c(which(grepl("day 3",df$pathway)),which(grepl("day 6",df$pathway)),which(grepl("neural",df$pathway)),which(grepl("day 15",df$pathway)))
#index<-c(which(grepl("day 3",df$pathway)),which(grepl("day 6",df$pathway)),which(grepl("day 15",df$pathway)))
df<-df[index,]
#Or 
df<-df

level<-c("neural GFP+","day 0 undifferentiated iPSC17","day 3 definitive endoderm",
         "day 6 anterior foregut endoderm","day 15 Nkx2-1-GFP-","day 15 Nkx2-1-GFP+")

df.new<-lapply(seq(level), function(i){
  state.df<-level[i]%>%as.character()
  df$variable<-as.character(df$variable)
  df.sub<-df%>%filter(variable=="mean_col"&type=="divergence disappear"&state==state.df)%>%
    arrange(as.numeric(value),decreasing=FALSE)
  df.sub.2<-df%>%filter(variable=="mean_col"&type=="divergence appear"&state==state.df)
  df.sub.2<-df.sub.2[match(df.sub$pathway,df.sub.2$pathway),]
  df.sub.3<-df%>%filter(variable=="GS.HESC"&type=="divergence disappear"&state==state.df)
  df.sub.3<-df.sub.3[match(df.sub$pathway,df.sub.3$pathway),]
  df.sub.4<-df%>%filter(variable=="GS.HESC"&type=="divergence appear"&state==state.df)
  df.sub.4<-df.sub.4[match(df.sub$pathway,df.sub.4$pathway),]
  df.sub.new<-rbind(df.sub,df.sub.2,df.sub.3,df.sub.4)%>%as.data.frame()
  return(df.sub.new)
})%>%bind_rows()

df.new.1<-df.new%>%filter(variable=="mean_col"&type=="divergence disappear")
df.new.1$pathwayNum<-seq(df.new.1$pathway)
df.new.2<-df.new%>%filter(variable=="mean_col"&type=="divergence appear")
df.new.2$pathwayNum<-seq(df.new.2$pathway)
df.new.3<-df.new%>%filter(variable=="GS.HESC"&type=="divergence disappear")
df.new.3$pathwayNum<-seq(df.new.3$pathway)
df.new.4<-df.new%>%filter(variable=="GS.HESC"&type=="divergence appear")
df.new.4$pathwayNum<-seq(df.new.4$pathway)

df.new<-rbind(df.new.1,df.new.2,df.new.3,df.new.4)

df<-df.new
#or
df<-df.new%>%filter(variable=="mean_col")

## boxplot and point plot

# pair-wised test for divergence end 
df.end.test<-df.new.1
df.end.test$state[grep("neural",df.end.test$state)]<-"neural"
df.end.test$state<-gsub("day 0 undifferentiated iPSC17","day0",df.end.test$state)
df.end.test$state<-gsub("day 3 definitive endoderm","endoderm",df.end.test$state)
df.end.test$state<-gsub("day 6 anterior foregut endoderm","endoderm",df.end.test$state)
df.end.test$state[grep("day 15",df.end.test$state)]<-"day15"

df.new.a<-df.end.test%>%filter(state=="day0")%>%
  arrange(as.numeric(value),decreasing=FALSE)%>%mutate(length=13)
df.new.a$pathwayNum<-seq(df.new.a$pathway)
df.new.b<-df.end.test%>%filter(state=="endoderm")%>%
  arrange(as.numeric(value),decreasing=FALSE)%>%mutate(length=21)
df.new.b$pathwayNum<-seq(df.new.b$pathway)+max(df.new.a$pathwayNum)
df.new.c<-df.end.test%>%filter(state=="neural")%>%
  arrange(as.numeric(value),decreasing=FALSE)%>%mutate(length=62)
df.new.c$pathwayNum<-seq(df.new.c$pathway)+max(df.new.b$pathwayNum)
df.new.d<-df.end.test%>%filter(state=="day15")%>%
  arrange(as.numeric(value),decreasing=FALSE)%>%mutate(length=111)
df.new.d$pathwayNum<-seq(df.new.d$pathway)+max(df.new.c$pathwayNum)

df.end.test<-rbind(df.new.a,df.new.b,df.new.c,df.new.d)

# t test
stat.test <- df.end.test %>%
  t_test(value ~ state, p.adjust.method = "bonferroni")
# Remove unnecessary columns and display the outputs
stat.test %>% select(-.y., -statistic, -df)

df.end.test$state <- factor(df.end.test$state, levels=unique(df.end.test$state))

my_comparisons <- list( c("day0","day15"),
                        c("day15","endoderm"), c("day15", "neural") )
						
# plot part1
ggplot(df.end.test, aes(x=state, y=value,group=state)) +
  stat_boxplot(geom = "errorbar", width=0.5, position = position_dodge(1))+
  geom_boxplot(aes(fill=state,alpha=0.5)) +
  stat_summary(fun.y=mean)+
  geom_point(position=position_dodge(width=0.75),size=4)+
  geom_signif(mapping=aes(x=as.factor(state),group=state),
              comparisons = my_comparisons,
              map_signif_level=F,
              step_increase = 0.05,
              test ="t.test")+
  ylim(0.2,0.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="differentiation state", y = "Normalized depth of divergence end")+
  theme(axis.text.x=element_text(size=48),axis.text.y=element_text(size=48),axis.title = element_text(size = 48))
# plot part2
df.end.test%>%ggplot() +
  stat_boxplot(aes(x=length, y=value,group=state),
               geom = "errorbar", width=0.5, position = position_dodge(1))+
  geom_boxplot(aes(x=length, y=value,,group=state,fill=state,alpha=0.5)) +
  geom_point(aes(x=pathwayNum, y=value,group=state,fill=state,size=4))+
  ylim(0.2,0.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="differentiation state", y = "Normalized depth of divergence end")+
  theme(axis.text.x=element_text(size=48),axis.text.y=element_text(size=48),axis.title = element_text(size = 48))
 # plot part3
 min.mean.sd.max<- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

ggplot(df.end.test, aes(x=state, y=value,group=state)) +
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot")+
  ylim(0.2,0.8)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="differentiation state", y = "Normalized depth of divergence end")+
  theme(axis.text.x=element_text(size=48),axis.text.y=element_text(size=48),axis.title = element_text(size = 48))

# fig 2c
library(tidyverse)
library(ape)
library(ggtree)
library(slingshot)
library(parallel)
library(monocle)
CBRAD5_Seuratobj<-readRDS("~/fig2/cbrad5_Seuratobj.Rds")
### step1
## data prepare before monocle
CBRAD5_data <- as(as.matrix(CBRAD5_Seuratobj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = CBRAD5_Seuratobj@meta.data)
fData <- data.frame(gene_short_name = row.names(CBRAD5_data), row.names = row.names(CBRAD5_data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(CBRAD5_data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

## hvgs selected by monocle
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)

# reduction
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree') 
# order
mycds <- orderCells(mycds,reverse = TRUE) 
# save for step 3 cds assignment
saveRDS(mycds,file = "~/fig2/G11_cbra_monocole.Rds")
# monocle State trajectory
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
## monocle Cluster trajectory
HESC.cbrad5.cluster<-A1_CBRAD5_Seuratobj@meta.data$HESC.cbrad5.cluster
plot2 <- plot_cell_trajectory(mycds, color_by = "HESC.cbrad5.cluster")
## monocle Pseudotime trajectory
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")

##monocle plot
plotc <- plot1|plot2|plot3

## pseudotime data
pData(mycds)

# ========================== step 2 save and load vars and functions ==============================
tree <- read.tree("~/fig1/G11-CBRAD5/G11-CBRAD5.nwk")
allele.info <- read.csv("~/fig1/G11-CBRAD5.AllelesInfo.txt",
                        sep = '\t', stringsAsFactors = F)
# get oneCell node pair
oneCellNode <- tree$tip.label[grepl("_1$", tree$tip.label)]
multiCellNode <- tree$tip.label[!grepl("_1$", tree$tip.label)]
oneCellTree <- drop.tip(tree, multiCellNode, trim.internal = FALSE)
oneCellNode.pair <- combn(oneCellNode, 2)
oneCellNode.pair <- data.frame(node1=oneCellNode.pair[1,], node2=oneCellNode.pair[2,], stringsAsFactors = F)

#oncellnode.tree <- drop.tip(tree, )
# get tree tibble
inSubtrees <- subtrees(oneCellTree)
distMatrix <- dist.nodes(oneCellTree)
treeTibble <- as_tibble(oneCellTree)

getDist <- function(nodeName, node){
  n2 <- dplyr::filter(treeTibble, label == nodeName) %>% `[`("node") %>% unlist()
  dis <- distMatrix[node, n2]
  return(dis)
}

#### functions ####
get.node.pdata <- function(raw.pdata, alleleinfo.df, cut.str = "G11-CBRAD5_"){
  rownames(raw.pdata) <- gsub(paste0(cut.str, "|-1"), "", rownames(raw.pdata))
  node.pdata <- raw.pdata[alleleinfo.df$BC %>% unlist(), , drop=F]
  node.pdata$NodeName <- alleleinfo.df$NodeName[match(rownames(node.pdata), alleleinfo.df$BC)]  # pacbio insection from 10X caputrued cells 
  return(node.pdata)
}

get.subtree.pseudotime.vars <- function(tree.obj, node.pdata.df){
  # create sub tree list and calculate distance
  subtree.list <- subtrees(tree.obj)
  distMatrix <- dist.nodes(tree.obj)
  treeTibble <- as_tibble(tree.obj)
  
  getDist <- function(nodeName, node){
    n2 <- dplyr::filter(treeTibble, label == nodeName) %>% `[`("node") %>% unlist()
    dis <- distMatrix[node, n2]
    return(dis)
  }
  
  subtree.pseudotime.vars <- 
    mclapply(seq(length(subtree.list)), function(s){
      # get tree depth 
      subtree <- subtree.list[[s]]
      s.name <- subtree$name
      treeDepth <- Map(getDist, subtree$tip.label, s.name) %>% unlist() %>% as.vector() %>% max()
      # calculate subtree pseudotime vars
      daughters <- subtree$tip.label %>% unlist()
      daughter.pseudotime <- node.pdata.df %>% filter(., NodeName %in% daughters) %>% `[`("Pseudotime") %>% unlist()
      sd.value <- sd(daughter.pseudotime)
      mean.value <- mean(daughter.pseudotime)
      return(data.frame(subtree = s.name,
                        subtree.depth = treeDepth,
                        subtree.leaves = length(daughters),
                        p.time.mean = mean.value,
                        p.time.sd = sd.value,
                        p.time.cv = sd.value / mean.value * 100,
                        stringsAsFactors = F))
    }, mc.cores = 30) %>% bind_rows()
  return(subtree.pseudotime.vars)
}

sim.subtree.ptime.vars <- function(node.pdata.df, tree.obj, ran.times=100){
  
  # create sub tree list and calculate distance
  subtree.list <- subtrees(tree.obj)
  distMatrix <- dist.nodes(tree.obj)
  treeTibble <- as_tibble(tree.obj)
  
  getDist <- function(nodeName, node){
    n2 <- dplyr::filter(treeTibble, label == nodeName) %>% `[`("node") %>% unlist()
    dis <- distMatrix[node, n2]
    return(dis)
  }
  
  get.subtree.pseudotime.vars.sim <- function(subtree.list, node.pdata.df){
    subtree.pseudotime.vars <- 
      mclapply(seq(length(subtree.list)), function(s){
        # get tree depth 
        subtree <- subtree.list[[s]]
        s.name <- subtree$name
        treeDepth <- Map(getDist, subtree$tip.label, s.name) %>% unlist() %>% as.vector() %>% max()
        # calculate subtree pseudotime vars
        daughters <- subtree$tip.label %>% unlist()
        daughter.pseudotime <- node.pdata.df %>% filter(., NodeName %in% daughters) %>% `[`("Pseudotime") %>% unlist()
        sd.value <- sd(daughter.pseudotime)
        mean.value <- mean(daughter.pseudotime)
        return(data.frame(subtree = s.name,
                          subtree.depth = treeDepth,
                          subtree.leaves = length(daughters),
                          p.time.mean = mean.value,
                          p.time.sd = sd.value,
                          p.time.cv = sd.value / mean.value * 100,
                          stringsAsFactors = F))
      }, mc.cores = 30) %>% bind_rows()
    return(subtree.pseudotime.vars)
  }
  
  #####
  random.ptime.results <- list()
  for (ran in seq(ran.times)) {
    sim.time <- sample(node.pdata.df$Pseudotime, size = nrow(node.pdata.df), replace = FALSE)
    sim.pdata <- node.pdata.df %>% dplyr::select(-Pseudotime)
    sim.pdata$Pseudotime <- sim.time
    sim.subtree.vars <- 
      get.subtree.pseudotime.vars.sim(subtree.list = subtree.list,
                                      node.pdata.df = sim.pdata)
    random.ptime.results[[ran]] <- sim.subtree.vars
  }
  return(random.ptime.results)
}

## ======================step3  monocle pseudotime ====================

cds <- readRDS("~/fig2/G11_cbrad5_monocole.Rds") ## for g11 a1 g2 from code of fig

mono.pdata <- cds@phenoData@data[, "Pseudotime", drop=F]

mono.node.pdata <- get.node.pdata(raw.pdata = mono.pdata,
                                  alleleinfo.df = allele.info,
                                  cut.str = "G11-CBRAD5_")

colnames(mono.node.pdata)[1] <- "Pseudotime"

mono.subtree.ptime.vars <- 
  get.subtree.pseudotime.vars(tree.obj = tree,
                              node.pdata.df = mono.node.pdata)
							  
# plot results
mono.subtree.ptime.vars %>% 
  ggplot() +
  geom_boxplot(aes(factor(subtree.depth), y=p.time.cv))

### random simmulation

mono.sim.subtree.ptime.vars <- 
  sim.subtree.ptime.vars(node.pdata.df = mono.node.pdata,
                         tree.obj = tree,
                         ran.times = 100)
# plot results
mono.sim.subtree.ptime.vars[[20]] %>% 
  ggplot() +
  geom_boxplot(aes(factor(subtree.depth), y=p.time.cv))

#### mean of 100 simulated tree
mono.sim.subtree.ptime.vars.new<-do.call("rbind", mono.sim.subtree.ptime.vars)
mono.sim.subtree.ptime.vars.mean<-aggregate(mono.sim.subtree.ptime.vars.new[,2:6], list(mono.sim.subtree.ptime.vars.new$subtree), FUN=mean)
colnames(mono.sim.subtree.ptime.vars.mean)[1]<-"subtree"

select.sim <- mono.sim.subtree.ptime.vars.mean
select.sim$type <- "simulated"

mono.plot.df <- mono.subtree.ptime.vars
mono.plot.df$type <- "real"

mono.plot.df <- bind_rows(mono.plot.df, select.sim)

mono.plot.df.real<-mono.plot.df[mono.plot.df$type=="real",]%>%`[`(,c(1,2,6,7))
mono.plot.df.simulated<-mono.plot.df[mono.plot.df$type=="simulated",]%>%`[`(,c(1,2,6,7))
mono.plot.df.new<-merge(mono.plot.df.real,mono.plot.df.simulated,by="subtree")
colnames(mono.plot.df.new)<-c("subtree","depth","V_obs","observe type","depth","V_exp","expected type")
mono.plot.df.new$group="G11 CBRAD5"

###function for p.time variance comparison of subtree between observed and expected tree###
p.time.var.obs.exp.conunt<-function(mono.plot.value){
  sum.poiont.obs.exp<-lapply(seq(unique(mono.plot.value$depth)),function(i){
  a<-mono.plot.value[which(mono.plot.value$depth==i),]
  b<-length(which(a$V_obs<a$V_exp))
  c<-nrow(a)
  d<-binom.test(b, c, 1/2,alternative = "greater")
  e<-i
  return(data.frame(smaller=b,biger=c-b,pvalue=d$p.value,depth=e))
  })%>% bind_rows(.)
  }

g11.sum.ptime.var<-p.time.var.obs.exp.conunt(g11_cbrad5_mono.plot.value)
  
g11.big.df<-g11.sum.ptime.var[,c("biger","pvalue","depth")]
g11.big.df$type<-"biger"
g11.big.df$sample<-"G11 CBRAD5"
colnames(g11.big.df)<-c("counts","pvalue","depth","type","sample")
g11.small.df<-g11.sum.ptime.var[,c("smaller","pvalue","depth")]
g11.small.df$type<-"smaller"
g11.small.df$sample<-"G11 CBRAD5"
colnames(g11.small.df)<-c("counts","pvalue","depth","type","sample")
g11.sum.ptime.var.plot<-rbind(g11.big.df,g11.small.df)

allcbrad5.sum.ptime.var.plot<-rbind(a1.sum.ptime.var.plot,g11.sum.ptime.var.plot,g2.sum.ptime.var.plot)
#sum(g11.sum.ptime.var$smaller)
#sum(g11.sum.ptime.var$smaller)+sum(g11.sum.ptime.var$biger)
# a1 cbrad5 variation of fate
binom.test(297, 389, 1/2,alternative = "greater") #p-value < 2.2e-16
# g2 cbrad5 variation of fate
binom.test(110, 146, 1/2,alternative = "greater") #p-value = 3.274e-10
# g11 cbrad5 variation of fate
binom.test(324, 441, 1/2,alternative = "greater") #p-value < 2.2e-16

mono.plot.value.merge<-rbind(a1_cbrad5_mono.plot.value,g2_cbrad5_mono.plot.value,g11_cbrad5_mono.plot.value)
mono.plot.value.merge$group<-factor(mono.plot.value.merge$group)

ggplot(mono.plot.value.merge,aes(x=V_exp,y=V_obs,color=group)) + 
  geom_point(size=8,alpha=0.7)+
  scale_color_manual(values = c("#F46D43","#4DAF4A","#74ADD1"))+
  stat_density2d(colour="#BDBDBD",bins = 10,size=1)+
  theme_bw()+
  theme(text = element_text(size =42),axis.text.x=element_text(size=42),axis.text.y=element_text(size=42))+ #legend.position = "top"
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  scale_x_continuous(expand = c(0,0),limits = c(0, 60))+
  scale_y_continuous(expand = c(0,0),limits = c(0, 100))+
  theme_classic()+
  xlab("p.times cv of simulated subtree") + ylab("p.times cv of real subtree")+
  #ggtitle("P.time variance comparison of subtree bewteen observed and expected tree")+
  theme(axis.text = element_text(size = 15))+
  theme(text = element_text(size = 20))
  