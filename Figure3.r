suppressMessages({
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
})

tip.long<-readRDS("~/fig1/all_cbrad5_gs_hesc_tree_modify_tip_long_new.Rds")
tip.long$sample<-gsub("_.*","",tip.long$intern)
tip.long$to<-gsub("-",".",tip.long$to)
colnames(tip.long)[1:2]<-c("NodeName","subtree")
tree.infos <- readRDS("~/fig1/all_cbrad5_GS_hesc_tree_dataframe_modify_new.Rds")
modify.all.cbra.gshesc.nodeExp<-readRDS("~/fig2/modify.all.cbra.gshesc.nodeExp.Rds")

gene.exp<-modify.all.cbra.gshesc.nodeExp

### cv calculate for all celltype in all sample (including A1,G11,G2 CBRAD5 and GS-HESC)

tip.long$celltype<-plyr::mapvalues(x = tip.long$celltype,
                                   from = c("3","0","7","6","10", "5","2","11","8","4","9","1"),
                                   to = c(paste0("C",seq(1,10)), "R1", "R2"))

								   
tree.infos$celltype<-plyr::mapvalues(x = tree.infos$celltype,
                                     from = c("3","0","7","6","10", "5","2","11","8","4","9","1"),
                                     to = c(paste0("C",seq(1,10)), "R1", "R2"))

### for hesc sample
# celltype: "R2","C1","C2",C3","C4","C5"

celltype<-"C2"
Allsample<-"GS-HESC"

## loop for hesc sample

#for (n in seq(Allsample)){
sample.df<-Allsample
#cv.intra.by.celltype<-lapply(seq(celltype),function(i){
df.cell.type<-celltype
df.tree.infos<-filter(tip.long,celltype==df.cell.type&sample==sample.df) 
df.gene.exp<-gene.exp[,colnames(gene.exp)%in%df.tree.infos$NodeName]
## filter gene based on cell number
gene.filter<-mclapply(seq(nrow(df.gene.exp)),mc.cores=40,function(x){
  df.gene<-rownames(df.gene.exp)[x]
  if (sum(df.gene.exp[x,]>0,na.rm = TRUE)>3){     ## raw filter cutoff is ">3"
    gene.list<-df.gene
  }else{gene.list<-NA}
  gene.list<-data.frame(gene.list)
  rownames(gene.list)<-df.gene
  colnames(gene.list)<-"gene"
  return(gene.list)
})%>%bind_rows()

gene.filter<-na.omit(gene.filter)
gene.exp.filter<-df.gene.exp[rownames(df.gene.exp)%in%rownames(gene.filter),]

### exp freq calculate
expression.frequency<-mclapply(seq(nrow(gene.exp.filter)),mc.cores=40,function(x){
  this.gene.exp.freq<-sum(gene.exp.filter[x,1:ncol(gene.exp.filter)]>0)
  this.gene.exp.freq<-data.frame(this.gene.exp.freq)
  rownames(this.gene.exp.freq)<-rownames(gene.exp.filter)[x]
  colnames(this.gene.exp.freq)<-"expression.frequency"
  return(this.gene.exp.freq)
})%>%bind_rows()

path.filename<-paste0("~/fig2/",sample.df,".expression.frequency.filter3.",df.cell.type,".Rds")
saveRDS(expression.frequency,path.filename)

## for cbrad5 sample
celltype<-c("R1","R2","C1","C2","C3","C4","C5","C6","C7","C9","C10")
Allsample<-c("A1-CBRAD5","G11-CBRAD5","G2-CBRAD5")

## loop for cbrad5 samples

# for (n in seq(Allsample)){
sample.df<-Allsample[1]
cv.intro.by.celltype<-lapply(seq(celltype),function(i){
  df.cell.type<-celltype[i]
  df.tree.infos<-filter(tip.long,celltype==df.cell.type&sample==sample.df) 
  df.gene.exp<-gene.exp[,colnames(gene.exp)%in%df.tree.infos$NodeName]
  ## filter gene based on cell number
  
  gene.filter<-mclapply(seq(nrow(df.gene.exp)),mc.cores=60,function(x){
    df.gene<-rownames(df.gene.exp)[x]
    if (sum(df.gene.exp[x,]>0,na.rm = TRUE)>1){     ## raw filter cutoff is ">3"
      gene.list<-df.gene
    }else{gene.list<-NA}
    gene.list<-data.frame(gene.list)
    rownames(gene.list)<-df.gene
    colnames(gene.list)<-"gene"
    return(gene.list)
  })%>%bind_rows()
  
  gene.filter<-na.omit(gene.filter)
  gene.exp.filter<-df.gene.exp[rownames(df.gene.exp)%in%rownames(gene.filter),]
  gene<-rownames(df.gene.exp)
  
  ### exp freq calculate
  expression.frequency<-mclapply(seq(nrow(gene.exp.filter)),mc.cores=60,function(x){
    this.gene.exp.freq<-sum(gene.exp.filter[x,1:ncol(gene.exp.filter)]>0)
    this.gene.exp.freq<-data.frame(this.gene.exp.freq)
    rownames(this.gene.exp.freq)<-rownames(gene.exp.filter)[x]
    colnames(this.gene.exp.freq)<-"expression.frequency"
    return(this.gene.exp.freq)
  })%>%bind_rows()
  
  path.filename<-paste0("~/fig2/",sample.df,".expression.frequency.filter3.",df.cell.type,".Rds")
  saveRDS(expression.frequency,path.filename)
  
  ## function
  
  ### combined internal node (common ancestor of sub_lineage leaves) with all gene's expression data of every leaf
  
  all.subtree.filter<-tip.long[tip.long$NodeName%in%colnames(df.gene.exp),]
  
  new.all.subtree.genes <- 
    mclapply(seq(nrow(all.subtree.filter)),mc.cores=60,function(x){                                            
      #this.gene.exp<-gene.exp[x,,drop=F]
      this.sub.exp <- gene.exp.filter[, match(all.subtree.filter$NodeName[x], colnames(gene.exp.filter))]%>%t()
      re.df <-
        data.frame(new.exp = this.sub.exp)
      colnames(re.df) <- rownames(gene.exp.filter)
      rownames(re.df)<-all.subtree.filter$NodeName[x]
      return(re.df)}) %>% bind_rows()
  
  all.subtree.gene.exp<-bind_cols(all.subtree.filter,new.all.subtree.genes)                                            
  
  #### calculate cv of per gene in subtree for real tree and simulated tree(average of 1000 tests)
  
  # for real tree
  
  # cv of intra lineage
  all.subtrees <- unique(all.subtree.gene.exp$subtree)%>%as.character()
  subtree.cv.real.tree <- 
    mclapply(seq(length(all.subtrees)),mc.cores=60,function(x){  
      subtree.this<-all.subtrees[x]
      subtree.gene.exp.df <-
        all.subtree.gene.exp %>% dplyr::filter(subtree == subtree.this) %>% `[`(,-c(1:5))
      cv.pergene<-lapply(seq(ncol(subtree.gene.exp.df)), function(i){
        cv<-sd(unlist(subtree.gene.exp.df[,i]))/mean(unlist(subtree.gene.exp.df[,i]))
        if(is.nan(cv)){cv<-0}
        cv.re <- data.frame(cv=cv)
        colnames(cv.re) <- colnames(subtree.gene.exp.df)[i]
        return(cv.re)
      }) %>% bind_cols()
      rownames(cv.pergene)<-subtree.this
      return(cv.pergene)
    }) %>% bind_rows()
  
  rownames(subtree.cv.real.tree)<-all.subtrees
  subtree.cv.real.tree<-na.omit(subtree.cv.real.tree)
  
  min.cv.real.tree<-lapply(seq(ncol(subtree.cv.real.tree)),function(i){
    min.cv.this.gene<-min(subtree.cv.real.tree[,i][subtree.cv.real.tree[,i]>0])%>%as.data.frame()
    colnames(min.cv.this.gene)<-colnames(subtree.cv.real.tree)[i]
    return(min.cv.this.gene)})%>%bind_cols()
  
  
  path.filename<-paste0("~/fig2/",sample.df,".subtree.cv.real.tree.intro.filter3.",df.cell.type,".Rds")
  saveRDS(subtree.cv.real.tree,path.filename)
  
  path.filename<-paste0("~/fig2/",sample.df,".min.cv.real.tree.intro.filter3.",df.cell.type,".Rds")
  saveRDS(min.cv.real.tree,path.filename)
  
  # for random tree
  # cv of intra lineage 
  min.cv.random.tree<-mclapply(seq(1:1000),mc.cores=60,function(f){
    random.gene.exp<-gene.exp.filter[,sample(1:ncol(gene.exp.filter),replace=T)]
    colnames(random.gene.exp)<-colnames(gene.exp.filter)
    sim.gene.exp.df<-random.gene.exp
    sim.subtree.gene.exp<-lapply(seq(nrow(all.subtree.filter)),function(x){                              
      this.sub.exp <- sim.gene.exp.df[, match(all.subtree.filter$NodeName[x], colnames(sim.gene.exp.df))]%>%t()
      re.df <-data.frame(new.exp = this.sub.exp)
      colnames(re.df) <-rownames(gene.exp.filter)
      rownames(re.df)<-all.subtree.filter$NodeName[x]
      return(re.df)}) %>% bind_rows()
    sim.subtree.gene.exp.cbind<-bind_cols(all.subtree.filter,sim.subtree.gene.exp)        
    all.subtrees <- unique(sim.subtree.gene.exp.cbind$subtree)                     
    subtree.cv.sim.tree<-lapply(seq(length(all.subtrees)),function(x){  
      subtree.this<-all.subtrees[x]
      sim.subtree.gene.exp.df <-sim.subtree.gene.exp.cbind %>% dplyr::filter(subtree == subtree.this) %>% `[`(,-c(1:5))  
      cv.pergene.sim<-lapply(seq(ncol(sim.subtree.gene.exp.df)),function(i){
        cv.sim<-sd(unlist(sim.subtree.gene.exp.df[,i]))/mean(unlist(sim.subtree.gene.exp.df[,i]))
        if(is.nan(cv.sim)){cv.sim<-0}
        cv.sim.re <- data.frame(cv.sim=cv.sim)
        colnames(cv.sim.re) <-colnames(sim.subtree.gene.exp.df)[i]
        # rownames(cv.re) <- 1
        return(cv.sim.re)
      }) %>% bind_cols()
      rownames(cv.pergene.sim)<-subtree.this
      return(cv.pergene.sim)
    }) %>% bind_rows()
    subtree.cv.sim.tree<-na.omit(subtree.cv.sim.tree)
    min.cv.sim.tree<-lapply(seq(ncol(subtree.cv.sim.tree)),function(i){
      sim.min.cv.this.gene<-min(subtree.cv.sim.tree[,i][subtree.cv.sim.tree[,i]>0])%>%as.data.frame()
      colnames(sim.min.cv.this.gene)<-colnames(subtree.cv.sim.tree)[i]
      return(sim.min.cv.this.gene)})%>%bind_cols()
    return(min.cv.sim.tree)
  })
  
  path.filename<-paste0("~/fig2/",sample.df,".min.cv.random.tree.intro.replace_sim.filter3.",df.cell.type,".Rds")
  saveRDS(min.cv.random.tree,path.filename)
  
  min.cv.random.tree.dataframe <-do.call("rbind",min.cv.random.tree)
  path.filename<-paste0("~/fig2/",sample.df,".min.cv.random.tree.intro.replace_sim.filter3.",df.cell.type,".dataframe.Rds")
  saveRDS(min.cv.random.tree.dataframe,path.filename)
  
  average.min.cv.random.tree<-colMeans(min.cv.random.tree.dataframe,na.rm = TRUE)%>%t()
  path.filename<-paste0("~/fig2/",sample.df,".min.cv.average.random.tree.intro.lineage.replace_sim.filter3.",df.cell.type,".Rds")
  saveRDS(average.min.cv.random.tree,path.filename)
  
  ### redefine for intra CV
  
  ### redefined the average.min.cv.random.tree, deal with the problem of "INF"
  
  min.cv.random.tree.dataframe.inf<-min.cv.random.tree.dataframe[,is.infinite(colMeans(min.cv.random.tree.dataframe))]
  min.cv.random.tree.dataframe.inf[mapply(is.infinite, min.cv.random.tree.dataframe.inf)] <- NA
  min.cv.random.tree.dataframe.inf.mean<-colMeans(min.cv.random.tree.dataframe.inf, na.rm = TRUE)%>%as.data.frame()
  
  average.min.cv.random.tree.normal<-average.min.cv.random.tree[,!is.infinite(colSums(average.min.cv.random.tree))]%>%as.data.frame()
  average.min.cv.random.tree<-rbind(average.min.cv.random.tree.normal,min.cv.random.tree.dataframe.inf.mean)
  
  gene.filter.n<-gene.filter
  
  # for lineage intra cv
  min.cv.real.random.combine<-rbind(min.cv.real.tree,average.min.cv.random.tree%>%t()%>%as.data.frame())
  rownames(min.cv.real.random.combine)<-c("mincv.real.tree","average.mincv.random.tree")
  min.cv.real.random.combine<-min.cv.real.random.combine%>%t()%>%as.data.frame()
  
  min.cv.real.random.combine$expression.frequency<-expression.frequency$expression.frequency[match(rownames(min.cv.real.random.combine),rownames(expression.frequency))]
  
  min.cv.real.random.combine$memory.index<-min.cv.real.random.combine$average.mincv.random.tree-min.cv.real.random.combine$mincv.real.tree ## intro lineage subtree CV of real tree is smaller than random tree
  
  min.cv.real.random.combine.filter<-min.cv.real.random.combine[rownames(min.cv.real.random.combine)%in%gene.filter.n$gene,]
  
  path.filename<-paste0("~/fig2/",sample.df,".min.cv.real.random.combine.filter3.",df.cell.type,".Rds")
  saveRDS(min.cv.real.random.combine.filter,path.filename)

 })

  # cell type count statistic of samples
  #
  celltype<-c("R1","R2","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")
  Allsample<-c("A1-CBRAD5","G11-CBRAD5","G2-CBRAD5","GS-HESC")
  
  cell.count.sample<-mclapply(seq(Allsample),mc.cores=50,function(f){
    sample.df<-Allsample[f]
    df<-lapply(seq(celltype),function(i){
      df.cell.type<-celltype[i]
      df.tree.infos<-tree.infos%>%filter(celltype==df.cell.type&sample==sample.df)%>%as.data.frame()
      celltype.df<-df.tree.infos%>%dplyr::select(label,celltype,sample)%>%as.data.frame()
      cell.count.df<-nrow(celltype.df)%>%as.data.frame()
      cell.count.df$celltype<-df.cell.type
      colnames(cell.count.df)<-c("cell.count","celltype")
      return(cell.count.df)
    })%>%bind_rows()
    df$sample<-sample.df
    return(df)
  })%>%bind_rows()
  
  ##### combine min cv of real and random tree group by celltype
  
  celltype<-c("R1","C6","C7","C9","C10")  ## from A1 G2 G11 CBRAD5 samples 
# or 
  celltype<-c("C1","C2","C3","C4") ## from GS HESC sample
# or
  celltype<-c("R2") ## from A1 G2 G11 CBRAD5 and GS HESC samples
  
  ### mean by sample and add the Ajusted memory index and filter by frequency
  
  for (i in seq(celltype)) {
    df.cell.type<-celltype[i]
    #sample<-c("A1-CBRAD5","G2-CBRAD5","G11-CBRAD5")
    #sample<-c("GS-HESC")
    sample<-c("A1-CBRAD5","G2-CBRAD5","G11-CBRAD5","GS-HESC")
    celltype.min.cv.name<-paste0(df.cell.type,".min.cv.real.random.combine.filter")
    min.cv.real.random.df<-mclapply(seq(sample),mc.cores=50,function(x){
      sample.df<-sample[x]
      path.filename<-paste0("~/fig2/",sample.df,".min.cv.real.random.combine.filter3.",df.cell.type,".Rds")
      min.cv.real.random.combine<-readRDS(path.filename)%>%dplyr::select("mincv.real.tree","average.mincv.random.tree","expression.frequency","memory.index")
      min.cv.real.random.combine$sample<-sample.df
      min.cv.real.random.combine$celltype<-df.cell.type
      min.cv.real.random.combine$gene<-rownames(min.cv.real.random.combine)
      cell.count.df<-cell.count.sample%>%filter(celltype==df.cell.type&sample==sample.df)
      min.cv.real.random.combine$expression.frequency<-min.cv.real.random.combine$expression.frequency/cell.count.df$cell.count
      return(min.cv.real.random.combine)})%>%bind_rows()
    min.cv.real.random.df<-aggregate(cbind(mincv.real.tree,average.mincv.random.tree,
                                           expression.frequency,memory.index) ~ gene,data = min.cv.real.random.df,
                                     FUN = mean, na.rm = TRUE)
    min.cv.real.random.df$celltype<-df.cell.type
    min.cv.real.random.df$Ajusted.memory.index<-min.cv.real.random.df$memory.index/min.cv.real.random.df$mincv.real.tree ## raw memory index/real min cv
    min.cv.real.random.df<-filter(min.cv.real.random.df,
                                  expression.frequency>quantile(min.cv.real.random.df$expression.frequency,.10)) ## filter gene by 10 percentile of gene expression frequency
    assign(celltype.min.cv.name,min.cv.real.random.df)
  }
  
  ## plot figure by celltype with Ajusted memory index
  ## add cutoff of high memory gene
  all.min.cv.real.random.combine<-rbind(C1.min.cv.real.random.combine.filter,C2.min.cv.real.random.combine.filter,
                                        C3.min.cv.real.random.combine.filter,C4.min.cv.real.random.combine.filter,
                                        C6.min.cv.real.random.combine.filter,C7.min.cv.real.random.combine.filter,
                                        C9.min.cv.real.random.combine.filter,C10.min.cv.real.random.combine.filter,
                                        R1.min.cv.real.random.combine.filter,R2.min.cv.real.random.combine.filter)
  
  all.min.cv.real.random.combine.mean<-aggregate(cbind(mincv.real.tree,average.mincv.random.tree,expression.frequency,memory.index,Ajusted.memory.index) ~ gene,
                                                 data = all.min.cv.real.random.combine, FUN = mean, na.rm = TRUE)
  
  cutoff.intra<-quantile(all.min.cv.real.random.combine$Ajusted.memory.index, probs = c(.90)) # cutoff of memory gene
  
  celltype.min.cv<-c("C1","C2","C3","C4","C6","C7","C9","C10","R1","R2")
  
  for (i in seq(celltype.min.cv)){
    file.name<-paste0(celltype.min.cv[i],".min.cv.real.random.combine.filter") 
    new.file.name<-paste0("new.",file.name)
    df<-get(file.name)
    df$fill.color<-"<=90 percentile"
    df$fill.color[which(df$Ajusted.memory.index>cutoff.intra)]<-">90 percentile"
    assign(new.file.name,df)
  }
  
  new.C4.min.cv.real.random.combine.filter$state<-"pluripotent" # define for C1,C2,C3,C4
  new.C10.min.cv.real.random.combine.filter$state<-"progenitor" # define for C6,C7,C9,C10
  new.R2.min.cv.real.random.combine.filter$state<-"spontaneous" # define for R1,R2
  
  new.all.min.cv.real.random.combine<-rbind(new.C1.min.cv.real.random.combine.filter,new.C2.min.cv.real.random.combine.filter,
                                            new.C3.min.cv.real.random.combine.filter,new.C4.min.cv.real.random.combine.filter,
                                            new.C6.min.cv.real.random.combine.filter,new.C7.min.cv.real.random.combine.filter,
                                            new.C9.min.cv.real.random.combine.filter,new.C10.min.cv.real.random.combine.filter,
                                            new.R1.min.cv.real.random.combine.filter,new.R2.min.cv.real.random.combine.filter)
  
  library(RColorBrewer)
  library("viridis")  
  display.brewer.all()
  brewer.pal.info
  color=c("#EDF8B1","#C7E9B4","#7FCDBB","#41B6C4","#1D91C0","#225EA8")
  
  pluri.stage.min.cv.real.random.combine.filter<-rbind(new.C1.min.cv.real.random.combine.filter,new.C2.min.cv.real.random.combine.filter,
                                                       new.C3.min.cv.real.random.combine.filter,new.C4.min.cv.real.random.combine.filter)
  
  pluri.stage.min.cv.real.random.combine.filter.mean<-aggregate(cbind(mincv.real.tree,average.mincv.random.tree,
                                                                      expression.frequency,memory.index,
                                                                      Ajusted.memory.index) ~ gene, 
                                                                data = pluri.stage.min.cv.real.random.combine.filter,
                                                                FUN = mean, na.rm = TRUE)
  pluri.stage.min.cv.real.random.combine.filter$fill.color<-"<=90 percentile"
  index<-which(pluri.stage.min.cv.real.random.combine.filter$Ajusted.memory.index>cutoff.intra)
  pluri.stage.min.cv.real.random.combine.filter$fill.color[index]<-">90 percentile"
  pluri.stage.min.cv.real.random.combine.filter
  
  prog.stage.min.cv.real.random.combine.filter<-rbind(new.C6.min.cv.real.random.combine.filter,new.C7.min.cv.real.random.combine.filter,
                                                      new.C9.min.cv.real.random.combine.filter,new.C10.min.cv.real.random.combine.filter)
  
  prog.stage.min.cv.real.random.combine.filter.mean<-aggregate(cbind(mincv.real.tree,average.mincv.random.tree,
                                                                     expression.frequency,memory.index,
                                                                     Ajusted.memory.index) ~ gene, 
                                                               data = prog.stage.min.cv.real.random.combine.filter,
                                                               FUN = mean, na.rm = TRUE)
  prog.stage.min.cv.real.random.combine.filter$fill.color<-"<=90 percentile"
  index<-which(prog.stage.min.cv.real.random.combine.filter$Ajusted.memory.index>cutoff.intra)
  prog.stage.min.cv.real.random.combine.filter$fill.color[index]<-">90 percentile"
  
  ## figure for point of per gene 
  # raw point
  pR2.raster<-R2.min.cv.real.random.combine.filter %>% # p1,p2,p3,p4,p6,p7,p9,p10,pR1,pR2 were reprensented with C1,C2,C3,C4,C6,C7,C9,C10,R1,R2 celltype figure
# or
  p1<-pluri.stage.min.cv.real.random.combine.filter.mean%>%
    arrange(expression.frequency,decreasing=T)%>%
    ggplot(aes(x = average.mincv.random.tree, y = mincv.real.tree,color=expression.frequency))+
    # ggrastr::rasterise(geom_point(size=12))+
    scale_colour_gradientn(colours=color)+
    geom_abline(color = "black",linetype=2,size=0.5)+
    theme_bw()+
    theme_classic()+
    theme(axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank())+
    theme(text = element_text(size = 24),axis.text.x=element_text(size=24),axis.text.y=element_text(size=24))+
    labs(title=paste0("pluripotent"," intra-subtree "),x="Expression variability(CV) from randomized lineage assignment", y = "Expression variability(CV) using barocdes")
  
  p2<-prog.stage.min.cv.real.random.combine.filter.mean%>%
    arrange(expression.frequency,decreasing=T)%>%
    ggplot(aes(x = average.mincv.random.tree, y = mincv.real.tree,color=expression.frequency))+
    ggrastr::rasterise(geom_point(size=12))+
    scale_colour_gradientn(colours=color)+
    geom_abline(color = "black",linetype=2,size=0.5)+
    theme_bw()+
    theme_classic()+
    theme(axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank())+
    theme(text = element_text(size = 24),axis.text.x=element_text(size=24),axis.text.y=element_text(size=24))+
    labs(title=paste0("progenitor"," intra-subtree "),x="Expression variability(CV) from randomized lineage assignment", y = "Expression variability(CV) using barocdes")
  
  pdf("~/fig2/Bystage.intra.cv.deviation.geneFreq.filter.pdf",width =60,height = 30)
  p1+p2
  dev.off()
  
  pdf("~/fig2/allcelltype.intra.cv.deviation.geneFreq.filter.pointbig.pdf",width = 40,height = 80)
  ggarrange(p1,p2,p3,p4,p6,p7,p9,p10,pR1,pR2,ncol=2,nrow=5)
  dev.off()
  
  pdf("~/fig2/allcelltype.intra.cv.deviation.geneFreq.filter.pointbig.raster.pdf",width = 40,height = 80)
  ggarrange(p1.raster,p2.raster,p3.raster,p4.raster,p6.raster,p7.raster,p9.raster,p10.raster,pR1.raster,pR2.raster,ncol=2,nrow=5)
  dev.off()
  
  #####  figure for distribution of Ajusted memory index for all celltype
  P1.intra.frequcy<-ggplot(new.C1.min.cv.real.random.combine.filter, aes(Ajusted.memory.index, color=as.factor(fill.color),fill=as.factor(fill.color))) +
    geom_histogram(binwidth=1/2.5,position = "identity",alpha=0.6)+ #bins=nlevels(as.factor(new.C1.min.cv.real.random.combine.filter$Ajusted.memory.index))
    scale_color_manual(values = c("#525252","#D73027"))+
    scale_fill_manual(values = c("#BDBDBD", "#F46D43"))+ # for pluripotent stage
    #scale_fill_manual(values = c("#FFFFFF","#FFFFFF"))+ # for progenitor stage
    scale_y_cut(breaks=2000, which=1, scales=0.5)+
    xlim(-2,12)+
    geom_vline(aes(xintercept = cutoff.intra), 
               linetype="dashed",
               size = 0.5, 
               color = "black")+
    theme_bw()+
    theme_classic()+
    theme(axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank())+
    theme(text = element_text(size = 20),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))+
    labs(title="C1 intra-subtree",x="Ajusted Memory index", y = "Frequency observed(Genes)")
  
  pdf("~/fig2/C1.intra.Ajusted.memory.index.distribution.pdf",width = 12,height = 12,onefile=F)
  P1.intra.frequcy
  dev.off()
  # for all stage
  pluri.stage<-ggplot(pluri.stage.min.cv.real.random.combine.filter, aes(Ajusted.memory.index,color=as.factor(fill.color),fill=as.factor(fill.color))) +
    geom_histogram(binwidth=1/5,position = "identity",alpha=0.6)+
    scale_color_manual(values = c("#525252","#D73027"))+
    scale_fill_manual(values = c("#BDBDBD","#F46D43"))+
    scale_y_cut(breaks=c(2000,10000), which=c(1, 3), scales=c(0.5, 3))+
    xlim(-2,12)+
    geom_vline(aes(xintercept = cutoff.intra), 
               linetype="dashed",
               size = 0.5, 
               color = "black")+
    theme_bw()+
    theme_classic()+
    theme(axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank())+
    theme(text = element_text(size = 20),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))+
    labs(title="pluri stage intra-subtree",x="Ajusted Memory index", y = "Frequency observed(Genes)")
  
  prog.stage<-ggplot(prog.stage.min.cv.real.random.combine.filter, aes(Ajusted.memory.index, color=as.factor(fill.color),fill=as.factor(fill.color))) +
    geom_histogram(binwidth=1/5,position = "identity",alpha=0.6,)+
    scale_color_manual(values = c("#525252","#D73027"))+
    scale_fill_manual(values = c("#FFFFFF","#FFFFFF"))+
    xlim(-2,12)+
    scale_y_cut(breaks=c(2000,10000), which=c(1, 3), scales=c(0.5, 3))+
    geom_vline(aes(xintercept = cutoff.intra), 
               linetype="dashed",
               size = 0.5, 
               color = "black")+
    theme_bw()+
    theme_classic()+
    theme(axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank())+
    theme(text = element_text(size = 20),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))+
    labs(title="prog stage intra-subtree",x="Ajusted Memory index", y = "Frequency observed(Genes)")
  pdf("~/fig2/Bystage.Ajusted.memory.index.distribution.pdf",width = 24,height = 8,onefile=F)
  pluri.stage+prog.stage
  dev.off()
  
  df<-rbind(prog.stage.min.cv.real.random.combine.filter,pluri.stage.min.cv.real.random.combine.filter)
  df$fill<-paste0(df$state,df$fill.color)
  
  pdf("~/fig2/pluri.vs.prog.stage.Ajusted.memory.index.distribution.2.pdf",width = 12,height = 8,onefile=F)
  df %>%
    ggplot(aes(x = Ajusted.memory.index, fill = as.factor(fill),color=as.factor(fill.color))) +
    geom_histogram(aes(y=..count..),binwidth=1/5, position = 'stack',alpha=0.6)+
    scale_color_manual(values = c("#525252","#D73027"))+
    scale_fill_manual(values = c("#BDBDBD","#F46D43","#FFFFFF","#FFFFFF" ))+
    xlim(-2,12)+
    geom_vline(aes(xintercept = cutoff.intra), 
               linetype="dashed",
               size = 0.5, 
               color = "black")+
    theme_bw()+
    theme_classic()+
    theme(axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank())+
    theme(text = element_text(size = 20),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))+
    labs(title="prog stage vs pluri stage",x="Ajusted Memory index", y = "Frequency observed(Genes)")
  dev.off()
  
  ## cutoff set figure of Ajusted memory index 
  pdf("~/fig2/allcelltype.intra.Ajusted.memory.index.cutoff.pdf",width = 12,height = 12,onefile=F)
  
  ggplot(new.all.min.cv.real.random.combine, aes(Ajusted.memory.index, color=as.factor(fill.color),fill=as.factor(fill.color))) +
    geom_histogram(binwidth=1/5,position = "identity",alpha=0.4)+
    scale_color_manual(values = c("#525252","#D73027"))+
    scale_fill_manual(values = c("#BDBDBD","#F46D43"))+
    xlim(-2,12)+
    scale_y_cut(breaks=5000, which=1, scales=0.5)+
    geom_vline(aes(xintercept = cutoff.intra), 
               linetype="dashed",
               size = 0.5, 
               color = "black")+
    theme_bw()+
    theme_classic()+
    theme(axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank())+
    theme(text = element_text(size = 20),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))+
    labs(title="all celltype intra-subtree",x="Ajusted Memory index", y = "Frequency observed(Genes)")
  
  dev.off()
  
  # memory gene fraction comparaness with Ajusted memory index
  library(scales)
  
  ## intra cv analysis data
  
  df.intra <- new.all.min.cv.real.random.combine %>%
    filter(celltype%in%c("C1","C2","C3","C4","C6","C7","C9","C10"))%>%
    group_by(celltype)%>% 
    dplyr::count(fill.color)
  
  df.intra <- df.intra%>%group_by(celltype)%>%mutate(prop = n / sum(n))
  
  df.intra$celltype.f = factor(df.intra$celltype, levels=c('C1','C2','C3',"C4","C6","C7","C9","C10"))
  
  df.intra.test<-df.intra%>%filter(.,fill.color==">90 percentile")
  index<-which(df.intra.test$celltype%in%c("C1","C2","C3","C4"))
  df.intra.pluri<-df.intra.test[index,]
  df.intra.prog<-df.intra.test[-index,]
  t.test(df.intra.pluri$prop,df.intra.prog$prop) # p-value = 0.003885  mean of x   0.17470061  mean of y 0.04707607
  
  df.intra$state<-"NA"
  df.intra$state[which(df.intra$celltype%in%c("C1","C2","C3","C4"))]<-"pluri"
  df.intra$state[which(df.intra$celltype%in%c("C6","C7","C9","C10"))]<-"prog"
  df.intra$state.new<-paste0(df.intra$fill.color,".",df.intra$state)
  
  pdf("~/fig2/high_memory_gene_fraction_comparation.intra.Ajusted.pdf",width = 12,height = 8,onefile=F)
  
  ggplot(df.intra,aes(x = celltype.f, y = prop, color=as.factor(fill.color),fill = as.factor(state.new))) +
    geom_col(size=1,alpha=0.8)+
    scale_y_break(c(0.4, 0.8))+
    scale_color_manual(values = c("#525252","#D73027"))+
    scale_fill_manual(values = c("#BDBDBD","#FFFFFF","#F46D43","#FFFFFF" ))+ 
    theme_classic()+
    theme(axis.text.y.right = element_blank(),
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank())+
    theme(text = element_text(size = 20),axis.text.x=element_text(size=20),axis.text.y=element_text(size=20))+
    labs(title="proportaion of memory gene,p=0.003885",x="Celltype of different stage", y = "Proportion")
  dev.off()  
  
  ############################ memory gene KEGG analysis#############################
  
  ### memory gene filter by intersect of intra and inter cv analysis with cutoff.intra of Ajusted memory index
  
  # intra
  ## select memory gene for all cell type of different state(pluripotent:C1,C2,C3,C4;progenitor:C6,C7,C9,C10)
	
  memory.gene.df.intra<-C10.min.cv.real.random.combine.filter%>% 
    filter(Ajusted.memory.index>cutoff.intra)%>%
    arrange(desc(Ajusted.memory.index))
	
  memory.gene.df<-memory.gene.df.intra[!grepl("\\.",memory.gene.df.intra$gene),]
  
  C10.memory.gene.intra<-memory.gene.df%>%select(gene)
  
  ## enrichR set
  
  dbs <- listEnrichrDbs()
  
  dbs <- dbs[order(dbs$libraryName),]
  
  # dbs[grep("Human",dbs$libraryName),]$libraryName
  
  websiteLive <- getOption("enrichR.live")
  
  listEnrichrSites()
  
  setEnrichrSite("Enrichr") # Human genes ?
  
  dbs_tf<-dbs[grep("TF",dbs$libraryName),]$libraryName[c(2,6)] #"ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X" "TF_Perturbations_Followed_by_Expression" 
  dbs_go <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023") ##  3
  dbs_pw <- c("KEGG_2021_Human", "WikiPathways_2021_Human", "BioPlanet_2019","MSigDB_Hallmark_2020",  
              "Elsevier_Pathway_Collection","Reactome_2022","NCI-Nature_2016")  ##  7
  
  ## enrichR enriched gene set
  
  #for (i in (seq(celltype.min.cv[1:8]))){
  
  celltype.df<-celltype.min.cv[1:8][i]
  
  # # use the memory gene calculated with Ajusted memory index
  enriched.result.name<-paste0(celltype.df,".enriched.result")
  
  # use the memory gene calculated with Ajusted memory index
  gene.intra<-get(paste0(celltype.df,".memory.gene.intra"))
  #######
  #encirhced_result<-mclapply(seq(gene.list),mc.cores=50,function(x){ # cycle for enrichr function probably was easily interrupted due to network connection prob
  
  gene<-gene.intra$gene
  
  enriched_tf_chip <-enrichr(gene, dbs_tf)
  
  enriched_tf_chip.df<-do.call("rbind",enriched_tf_chip)%>%
    filter(Adjusted.P.value<0.05)%>%filter(!grepl("mm",Term))%>%filter(!grepl("mouse",Term))%>%
    filter(!grepl("MOUSE",Term))%>%dplyr::select(Term,Adjusted.P.value,P.value,Overlap,Genes,Odds.Ratio,Combined.Score)%>%
    mutate(celltype=celltype.df)%>%arrange(desc(Combined.Score))%>%
    distinct(Term, .keep_all = TRUE)
  
  enriched_tf_chip.df$Term<-tolower(enriched_tf_chip.df$Term)
  
  enriched_go<-enrichr(gene, dbs_go)
  enriched_go.df<-do.call("rbind",enriched_go)%>%filter(Adjusted.P.value<0.05)%>%
    filter(!grepl("mm",Term))%>%filter(!grepl("mouse",Term))%>%
    filter(!grepl("MOUSE",Term))%>%dplyr::select(Term,Adjusted.P.value,P.value,Overlap,Genes,Odds.Ratio,Combined.Score)%>%
    mutate(celltype=celltype.df)%>%arrange(desc(Combined.Score))
  
  enriched_pw<-enrichr(gene, dbs_pw)
  enriched_pw.df<-do.call("rbind",enriched_pw)%>%
    filter(Adjusted.P.value<0.05)%>%filter(!grepl("mm",Term))%>%filter(!grepl("mouse",Term))%>%
    filter(!grepl("MOUSE",Term))%>%dplyr::select(Term,Adjusted.P.value,P.value,Overlap,Genes,Combined.Score)%>%
    mutate(celltype=celltype.df)%>%arrange(desc(Combined.Score))%>%
    distinct(Term, .keep_all = TRUE)
  
  enriched.result<-list(enriched_tf_chip.df,enriched_go.df,enriched_pw.df)
  
  names(enriched.result)<-c("tf_chip","go","pw")
  #return(enriched.result)
  ##})
  
  enriched.result.all<-enriched.result
  assign(enriched.result.name,enriched.result.all) 
  #}
  
  ## plot enrichr result of TF
  
  celltype<-c("C1","C2","C3","C4","C6","C7","C9","C10","R1","R2")
  
  ## pluritpotent cell types
  memory.info.list<-list(C1.min.cv.real.random.combine.filter,
                         C2.min.cv.real.random.combine.filter,
                         C3.min.cv.real.random.combine.filter,
                         C4.min.cv.real.random.combine.filter)
  
  celltype.df.hesc<-celltype[1:4]
  names(memory.info.list)<-celltype.df.hesc
  
  celltype.tf.info.pluri<-lapply(seq(celltype.df.hesc[1:4]),function(i){  
    celltype.df<-celltype.df.hesc[i]
    enriched.result<-get(paste0(celltype.df,".enriched.result"))
    table.intra<-enriched.result[[1]]  #
    table.intra<-table.intra%>%filter(grepl("ChEA",rownames(table.intra)))    ## Perturbations or ChEA
    table.intra<-table.intra%>%mutate(newTerm=str_extract(Term,"(\\w+)"))%>%
      distinct(newTerm, .keep_all = TRUE)
    memory.info<-memory.info.list[[celltype.df]]
    table.intra$expression.frequence<-memory.info$expression.frequency[match(toupper(table.intra$newTerm),
                                                                             memory.info$gene)]
    table.intra<-table.intra%>%arrange(desc(Combined.Score))
    return((table.intra))
  })
  
  names(celltype.tf.info.pluri)<-celltype.df.hesc
  
  tf.hesc.info<-mclapply(seq(celltype.df.hesc),mc.cores=50, function(i){
    celltype.df<-celltype.df.hesc[i]
    memory.info<-memory.info.list[[celltype.df]]
    df.info<-celltype.tf.info.pluri[[celltype.df]]
    df.info<-#df.info[-which(grepl("znf",df.info$Term)),]%>%
      #filter(.,expression.frequence>=perc50)%>% # TF could be filtered by exp freq percentile or not, and the results are simliar
      df.info%>%mutate(newTerm=paste0(str_extract(Term,"(\\w+)"),"_",i))%>% 
      distinct(newTerm, .keep_all = TRUE)%>%arrange(desc(Combined.Score))%>%'['(1:10,)
    return(df.info)
  })%>%bind_rows()
  
  ## progenitor cell types
  memory.info.list<-list(C6.min.cv.real.random.combine.filter,
                         C7.min.cv.real.random.combine.filter,
                         C9.min.cv.real.random.combine.filter,
                         C10.min.cv.real.random.combine.filter)
  
  celltype.df.cbra<-celltype[5:8]
  
  names(memory.info.list)<-celltype.df.cbra
  
  celltype.tf.info.prog<-lapply(seq(celltype.df.cbra[1:4]),function(i){  
    celltype.df<-celltype.df.cbra[i]
    enriched.result<-get(paste0(celltype.df,".enriched.result"))
    table.intra<-enriched.result[[1]]
    table.intra<-table.intra%>%filter(grepl("Perturbations",rownames(table.intra)))    ## Perturbations database  or ChEA database,the results are similar
    table.intra<-table.intra%>%mutate(newTerm=str_extract(Term,"(\\w+)"))%>%
      distinct(newTerm, .keep_all = TRUE)
    memory.info<-memory.info.list[[celltype.df]]
    table.intra$expression.frequence<-memory.info$expression.frequency[match(toupper(table.intra$newTerm),
                                                                             memory.info$gene)]
    table.intra<-table.intra[order(table.intra$Adjusted.P.value),]
    return((table.intra))
  })
  
  names(celltype.tf.info.prog)<-celltype.df.cbra
  
  tf.cbra.info<-mclapply(seq(celltype.df.cbra),mc.cores=50, function(i){
    celltype.df<-celltype.df.cbra[i]
    memory.info<-memory.info.list[[celltype.df]]
    #perc50<-quantile(memory.info$expression.frequency,.50) # # TF could be filtered by exp freq percentile or not, and the results are simliar
    df.info<-celltype.tf.info.prog[[celltype.df]]
    df.info<-#df.info[-which(grepl("znf",df.info$Term)),]%>%
      #filter(.,expression.frequence>=perc50)%>% # TF exp freq >50 percentile
      df.info%>%mutate(newTerm=paste0(str_extract(Term,"(\\w+)"),"_",i+4))%>% 
      distinct(newTerm, .keep_all = TRUE)%>%arrange(desc(Combined.Score))%>%'['(1:10,)
    return(df.info)
  })%>%bind_rows()#%>%na.omit()
  
  ############# plot for enrichment 
  
  #df<-tf.info.com.mean
  df<-tf.cbra.info
  index.1<-which(df$Adjusted.P.value<0.001)
  index.2<-which(df$Adjusted.P.value>0.001&df$Adjusted.P.value<0.01) 
  index.3<-which(df$Adjusted.P.value>0.01&df$Adjusted.P.value<0.05) 
  df$sig<-" "
  df$sig[index.1]<-"***"
  df$sig[index.2]<-"**"
  df$sig[index.3]<-"*"
  
  df.1<-df
  df.2<-df
  df<-rbind(df.1,df.2)
  df$newTerm <- factor(df$newTerm, levels = unique(df$newTerm))
  df$celltype<-factor(df$celltype,levels = c("C1","C2","C3","C4","C6","C7","C9","C10"))
  #df$state <- factor(df$newTerm, levels = unique(df$state))
  pdf("~/fig2/all.memory.gene.all.celltype.enriched.tf.top10.nofilterExpfreq.Perturbations.pdf",height=35,width = 28,onefile=F)
  df%>%ggplot(aes(newTerm,Combined.Score,color=as.factor(celltype),fill=as.factor(celltype)))+
    geom_col(alpha=0.6,position = "dodge")+
    geom_text(aes(label = sig),hjust = -0.5,vjust=1,size=6,position = position_dodge(width = 1),color="black")+
    # scale_color_manual(values = c("#41B6C4","#D73027"))+
    # scale_fill_manual(values = c("#74ADD1", "#F46D43"))+
    scale_color_manual(values = c("#FEE090","#FDAE61","#F46D43","#D73027",
                                  "#4EB3D3","#2B8CBE","#0868AC","#084081"))+
    
    scale_fill_manual(values = c("#FEE090","#FDAE61","#F46D43","#D73027",
                                 "#4EB3D3","#2B8CBE","#0868AC","#084081"))+
    coord_flip() +
    theme(plot.title = element_text(size=30, hjust = 0.5),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size=30),
          axis.title = element_text(size=30))+
    xlab("Enriched TFs")+
    ylab("Enriched combined.Score")+
    labs(title = "Genes in response to perturbation")+ # ENCODE and ChEA Consensus TFs for hesc celltype 
    facet_grid(as.factor(celltype)~.,scales="free_y")
  # Genes in response to perturbation for cbra celltype
  dev.off()

