library(ape)
library(ggtree)
library(tidyverse)
library(parallel)
library(magrittr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggraph)
library(RColorBrewer)
library(ggimage)
library(ggpubr)
library(ggsignif) 
library(plotrix)

# fig 4a,b,c

# modified depth of subtree
mdepth.df<-read.csv("~/fig2/AllNode.mdepth.bind.csv")%>%'['(-1)
colnames(mdepth.df)[1]<-"subtree"
mdepth.df$sample<-gsub("-",".",mdepth.df$sample)

all.subtree.depth<-mdepth.df%>%filter(grepl("A1",sample))  ## loop for sample names of cbrad5: A1 G11 G2 

colnames(all.subtree.depth)[1]<-"subtree"

# filter subtree by depth>0.3 
# subtree with depth<0.3 was not enough size of leaf for cell type composition test

df<-all.subtree.depth
df$mdepth.group<-df$mdepth/max(df$mdepth) # normalized from 0 to 1

# for A1 CBRAD5 and G11 CBRAD5 sample
zone<-seq(from=0,to=1.0,by=0.1)
# OR only for G2 sample
zone<-seq(from=0.1,to=1.0,by=0.1)

df.new.1<-lapply(seq((length(zone)-1)),function(x){
  min<-zone[x]
  max<-zone[x+1]
  index<-which(df$mdepth.group>min&df$mdepth.group<=max)
  if(length(index)<1){
    df.new<-as.data.frame(matrix(nrow=1,ncol=5))
    colnames(df.new)<-colnames(df)
  }else{
    df.new<-df[index,]
    df.new$mdepth.group<-max%>%as.numeric()}
  return(df.new)
})%>%bind_rows()%>%na.omit()

df.new.all<-df.new.1

all.subtree.depth<-df.new.all%>%dplyr::select("subtree","mdepth","mdepth.group","sample")

index<-which(all.subtree.depth$mdepth.group>=0.4) # subtree with depth>0.3
root<-min(all.subtree.depth$subtree)
subtree.filter<-all.subtree.depth$subtree[index]
subtree.filter<-subtree.filter[-which(subtree.filter==root)]
subtree.filter.info<-all.subtree.depth[which(all.subtree.depth$subtree%in%subtree.filter),]

#==================================================for real tree===========================================================
tree.infos <- readRDS("~/fig1/all_cbrad5_GS_hesc_tree_dataframe_modify_new.Rds")

tree.infos$celltype<-plyr::mapvalues(x= tree.infos$celltype,from = c("3","0","7","6","10", "5","2","11","8","4","9","1"),
                                     to = c(paste0("C",seq(1,10)), "R1", "R2"))

## calculate for all samples respectively
tree <- read.tree("~/fig1/A1-CBRAD5.nwk")  ## sample name
inode_filter<-subtree.filter
### celltype info
celltype<-filter(tree.infos,sample=="A1-CBRAD5"&celltype!="inter")
## next step
celltype<-celltype[,colnames(celltype)%in%c("label","BC","celltype")]
colnames(celltype)<-c("NodeName","BC","HESC.cbrad5.cluster")
celltype$leaftype<-gsub(".*_","",celltype$NodeName)
celltype.onecell<-celltype[which(celltype$leaftype=="1"),]%>%.[,-ncol(.)] 

N<-nrow(celltype) # celltype or celltype.onecell

celltype_info<-table(celltype$HESC.cbrad5.cluster)%>%as.data.frame() # celltype or celltype.onecell

celltype_counts<-celltype_info[,2]%>%as.numeric()
celltype_fraction<-celltype_counts/N
celltype_f<-data.frame(freq=celltype_fraction,
                       celltype=celltype_info$Var1)
					   
# subtree celltype fraction or count
inSubtrees <- subtrees(tree)
all.subtree <- lapply(seq(length(inSubtrees)), function(s){
  subtree <- inSubtrees[[s]]
  s.name <- subtree$name
  tip<-subtree$tip.label
  return(data.frame(subtree = s.name,
                    NodeName= tip,
                    stringsAsFactors = F))
}) %>% bind_rows()

all.subtree$leaftype<-gsub(".*_","",all.subtree$NodeName)
all.subtree.onecell<-all.subtree[which(all.subtree$leaftype=="1"),]
subtree_celltype<-merge(all.subtree,celltype,by="NodeName") # all.subtree.onecell or all.subtree /celltype.onecell or celltype
subtree_celltype<-subtree_celltype[,-which(grepl("leaftype",colnames(subtree_celltype)))]
total.subtree_celltype<-subtree_celltype
subtree_celltype<-filter(subtree_celltype,subtree%in%inode_filter)  
subtree<-split(subtree_celltype,subtree_celltype$subtree)

## table the freq of cell type for all subtree

celltype_counts_internode<-mclapply(seq(length(subtree)),mc.cores=50, function(s){
  sub_subtree <- subtree[[s]]
  s.name <- unique(sub_subtree$subtree)
  celltype.tab<-xtabs(~HESC.cbrad5.cluster,sub_subtree)%>%as.data.frame()
  celltype.data<-data.frame(matrix(NA, nrow = 1, ncol = nrow(celltype.tab)))
  colnames(celltype.data)<-celltype.tab[,1]
  celltype.data[1,]<-celltype.tab[,2]%>%as.numeric() # cell type counts
  if(nrow(celltype.data)<12){  # total 12 cell types for CBRAD5
    df<-setdiff(celltype$HESC.cbrad5.cluster,colnames(celltype.data))
    celltype.data[,df]<-0
    celltype.data<-celltype.data[ , order(names(celltype.data))]
  }else{celltype.data<-celltype.data[ , order(names(celltype.data))]}
  rownames(celltype.data)<-s.name
  return(celltype.data)
})

celltype_counts_internode<-rbindlist(celltype_counts_internode, fill = TRUE)%>%as.data.frame()
rownames(celltype_counts_internode)<-names(subtree)

#===================================simulated random tree====================================
ran.times=1000
    sim.subtree_celltype<-lapply(seq(ran.times),function(i){
      sim.celltype <- sample(total.subtree_celltype$HESC.cbrad5.cluster, size = nrow(total.subtree_celltype), replace = F)
      sim.pdata <- total.subtree_celltype %>% dplyr::select(-HESC.cbrad5.cluster)
      sim.pdata$HESC.cbrad5.cluster <- sim.celltype
      sim.pdata<-subset(sim.pdata,sim.pdata$subtree%in%inode_filter)
      return(sim.pdata)
    })

## put 1000 simulate tree to calculate all cell type counts for all interal node
sim.celltype_counts_internode<-mclapply(seq(length(sim.subtree_celltype)), mc.cores=50,function(i){
  sim.subtree.celltype<-sim.subtree_celltype[[i]]
  sim.subtree<-split(sim.subtree.celltype,sim.subtree.celltype$subtree)
  sim.celltype_counts_internode<-lapply(seq(length(sim.subtree)), function(s){
    sub_subtree <- sim.subtree[[s]]
    s.name <- unique(sub_subtree$subtree)
    celltype.tab<-xtabs(~HESC.cbrad5.cluster,sub_subtree)%>%as.data.frame()
    celltype.data<-data.frame(matrix(NA, nrow = 1, ncol = nrow(celltype.tab)))
    colnames(celltype.data)<-celltype.tab[,1]
    celltype.data[1,]<-celltype.tab[,2]%>%as.numeric() # celltype counts
    if(nrow(celltype.data)<12){  # total 12 cell types for CBRAD5 sample
      df<-setdiff(celltype$HESC.cbrad5.cluster,colnames(celltype.data))
      celltype.data[,df]<-0
      celltype.data<-celltype.data[ , order(names(celltype.data))]
    }else{celltype.data<-celltype.data[ , order(names(celltype.data))]}
    rownames(celltype.data)<-s.name
    return(celltype.data)
  })
  
  sim.celltype_counts_internode<-rbindlist(sim.celltype_counts_internode, fill = TRUE)%>%as.data.frame()
  rownames(sim.celltype_counts_internode)<-names(sim.subtree) 
  return(sim.celltype_counts_internode)
})

# put 1000 simulated tree to calculate expected all cell type counts for all internal node based on chi-square test logic

#=====================================chisq.test================================================

exp_f<-function(x){sum(x)*celltype_f$freq}  # celltype_f celltype_f.all
expect<-t(apply(celltype_counts_internode,1,exp_f))%>%as.data.frame()
colnames(expect)<-colnames(celltype_counts_internode)

##########   for real tree    ##########
chisq.value<-mclapply(seq(nrow(celltype_counts_internode)), mc.cores=50,function(i){
  subtree_df1<-celltype_counts_internode[i,]
  subtree_df2<-expect[i,]  # for celltype counts chisq
  subtree_df<-rbind(subtree_df1,subtree_df2)
  n_per_subtree_per_celltype<-lapply(seq(ncol(subtree_df)), function(j){
    n_per_celltype_1<-subtree_df[1,j]
    n_per_celltype_2<-subtree_df[2,j]
    df_chisq<-(n_per_celltype_1-n_per_celltype_2)^2/n_per_celltype_2
    return(data.frame(
      chisq.value=df_chisq,
      celltype=colnames(subtree_df)[j],
      internalnode=rownames(subtree_df)[1]))
  })%>% bind_rows()
})%>% bind_rows()

chisq.value$group<-"real"
colnames(chisq.value)<-c("value","celltype","internalnode","group")

real_value<-sum(chisq.value$value) 
real_value<-real_value%>%as.data.frame()%>%mutate(.,type="real")
colnames(real_value)<-c("Chisquare_value","type")

##########  for simulation   ##########
sim.chisq.value<-mclapply(seq(length(sim.celltype_counts_internode)),mc.cores = 50,function(x){
  sim.celltype_counts_internode_subset<-sim.celltype_counts_internode[[x]]
  chisq.value<-lapply(seq(nrow(sim.celltype_counts_internode_subset)),function(i){
    subtree_df1<-sim.celltype_counts_internode_subset[i,]
    subtree_df2<-expect[i,]  ## for celltype counts chisq
    subtree_df<-rbind(subtree_df1,subtree_df2)
    n_per_subtree_per_celltype<-lapply(seq(ncol(subtree_df)), function(j){
      n_per_celltype_1<-subtree_df[1,j]
      n_per_celltype_2<-subtree_df[2,j]
      df_chisq<-(n_per_celltype_1-n_per_celltype_2)^2/n_per_celltype_2
      return(data.frame(
        chisq.value<-df_chisq,
        celltype<-colnames(subtree_df)[j],
        internalnode<-rownames(sim.celltype_counts_internode_subset)[i]))
    })%>%bind_rows()
  })%>% bind_rows()
  chisq.value$group<-"sim"
  colnames(chisq.value)<-c("value","celltype","internalnode","group")
  return(chisq.value)
})

# sum the chisq value of all subtree
all.sim.value<-lapply(seq(length(sim.chisq.value)),function(i){
  if ("try-error" %in% class(sim.chisq.value[[i]]))
  {sum_chisq.value<-NA}
  else
  {chisq.value<-sim.chisq.value[[i]]
  sum_chisq.value<-sum(chisq.value[,1])
  return(sum_chisq.value)
  }})%>%unlist

all.sim.value<-sort(all.sim.value,decreasing = FALSE)
all.sim.value<-all.sim.value%>%as.data.frame()
colnames(all.sim.value)<-"Chisquare_value"
all.sim.value$type<-"sim"

#================================plot========================================== 
## for total celltype  
all.sim.real.value<-rbind(all.sim.value,real_value)
all.sim.real.value<-all.sim.real.value[order(all.sim.real.value$Chisquare_value),]
p.value<-which(grepl("real",all.sim.real.value$type))/1000  

df<-data.frame(x1= real_value$Chisquare_value,x2= real_value$Chisquare_value,y1=10,y2=0)

pdf("~/fig4/chisq_of_count.rootRemove.>0.3.A1.CBRAD5.allcellnode.pdf",height = 10,width =10)
# Histogram with density plot
ggplot(all.sim.value, aes(x=Chisquare_value)) +
  geom_histogram(aes(y = ..count..), binwidth=30, fill="#92C5DE") +
  geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2,colour="#D73027"),data=df,arrow = arrow(length = unit(0.4,"cm")))+guides(colour=FALSE)+
  theme_bw()+
  theme(plot.title = element_text(size=20, hjust = 0.5),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))+
  labs(title="A1 CBRAD5 p.value=0.002")+
  ylab("Count")+
  xlab("chisq for terminal cell type compositions summed among sub-CLTs")
dev.off()  

# fig 4d

all.subtree.depth<-mdepth.df%>%filter(grepl("G2",sample))  ## loop for samples: A1 G11 G2
colnames(all.subtree.depth)[1]<-"subtree"
#==================================================for real tree===========================================================
## calculate for all samples respectively
tree <- read.tree("~/fig1/A1-CBRAD5.nwk")  ## loop for samples: A1 G11 G2
## inode filter 
### remove root inode 
inode_filter<-all.subtree.depth$subtree[-which(all.subtree.depth$subtree==min(all.subtree.depth$subtree))]
### celltype info
celltype<-filter(tree.infos,sample=="A1-CBRAD5"&celltype!="inter")  # sample name
## next step
celltype<-celltype[,colnames(celltype)%in%c("label","BC","celltype")]
colnames(celltype)<-c("NodeName","BC","HESC.cbrad5.cluster")
celltype$leaftype<-gsub(".*_","",celltype$NodeName)
celltype.onecell<-celltype[which(celltype$leaftype=="1"),]%>%.[,-ncol(.)] 

N<-nrow(celltype) # celltype or celltype.onecell
celltype_info<-table(celltype$HESC.cbrad5.cluster)%>%as.data.frame() # celltype or celltype.onecell
celltype_counts<-celltype_info[,2]%>%as.numeric()
celltype_fraction<-celltype_counts/N
celltype_f<-data.frame(freq=celltype_fraction,
                       celltype=celltype_info$Var1)
inSubtrees <- subtrees(tree)

all.subtree <- lapply(seq(length(inSubtrees)), function(s){
  subtree <- inSubtrees[[s]]
  s.name <- subtree$name
  tip<-subtree$tip.label
  return(data.frame(subtree = s.name,
                    NodeName= tip,
                    stringsAsFactors = F))
}) %>% bind_rows()

all.subtree$leaftype<-gsub(".*_","",all.subtree$NodeName)
all.subtree.onecell<-all.subtree[which(all.subtree$leaftype=="1"),]
subtree_celltype<-merge(all.subtree,celltype,by="NodeName") # all.subtree.onecell or all.subtree /celltype.onecell or celltype
subtree_celltype<-subtree_celltype[,-which(grepl("leaftype",colnames(subtree_celltype)))]
total.subtree_celltype<-subtree_celltype
subtree_celltype<-filter(subtree_celltype,subtree%in%inode_filter)  
subtree<-split(subtree_celltype,subtree_celltype$subtree)

## table the freq of cell type for all subtree

celltype_counts_internode<-mclapply(seq(length(subtree)),mc.cores=50, function(s){
  sub_subtree <- subtree[[s]]
  s.name <- unique(sub_subtree$subtree)
  celltype.tab<-xtabs(~HESC.cbrad5.cluster,sub_subtree)%>%as.data.frame()
  celltype.data<-data.frame(matrix(NA, nrow = 1, ncol = nrow(celltype.tab)))
  colnames(celltype.data)<-celltype.tab[,1]
  celltype.data[1,]<-celltype.tab[,2]/nrow(sub_subtree)%>%as.numeric() ###  but cell type ratio not cell type counts
  if(nrow(celltype.data)<12){  # total 12 cell types for CBRAD5 sample
    df<-setdiff(celltype$HESC.cbrad5.cluster,colnames(celltype.data))
    celltype.data[,df]<-0
    celltype.data<-celltype.data[ , order(names(celltype.data))]
  }else{celltype.data<-celltype.data[ , order(names(celltype.data))]}
  rownames(celltype.data)<-s.name
  return(celltype.data)
})

celltype_counts_internode<-rbindlist(celltype_counts_internode, fill = TRUE)%>%as.data.frame()
rownames(celltype_counts_internode)<-names(subtree)

### statistics number of sub-CLTs by celltype which cell type ratio  give rise to +/-0.1 of expected ratio

A1.celltype_f<-celltype_f # loop for samples: A1 G11 G2
A1.celltype_counts<-celltype_counts_internode # loop for samples: A1 G11 G2
#===================================================================
# mean for all samples A1 G11 G2
mean.exp.f<-cbind(A1.celltype_f$freq,G2.celltype_f$freq,G11.celltype_f$freq)%>%rowMeans()%>%as.data.frame()
mean.exp.f$celltype<-A1.celltype_f$celltype
colnames(mean.exp.f)<-c("freq","celltype")

df.counts<-A1.celltype_counts[,c("C6","C7","C9","C10")]
celltype_f<-mean.exp.f%>%filter(celltype%in%c("C6","C7","C9","C10"))

subtree.f.stat<-lapply(seq(nrow(celltype_f)),function(i){
  celltype<-celltype_f$celltype[i]
  subtree.df<-df.counts[,which(colnames(df.counts)==celltype)]
  exp.f<-celltype_f$freq[i]
  exp.f.right<-exp.f+0.1  # *** exp.f+/-0.1
  exp.f.left<-exp.f-0.1
  sub.match<-rownames(df.counts)[which(subtree.df>=exp.f.left&subtree.df<=exp.f.right)]
  data<-data.frame(subtree=sub.match,
                   celltype=celltype)
  return(data)
})%>%bind_rows()

A1.subtree.f.stat.0.10<-subtree.f.stat # loop for samples: A1 G11 G2

df.1<-subtree.f.stat%>%filter(celltype=="C6")
df.2<-subtree.f.stat%>%filter(celltype=="C7")
df.3<-subtree.f.stat%>%filter(celltype=="C9")
df.4<-subtree.f.stat%>%filter(celltype=="C10")

N<-length(Reduce(intersect,list(df.1$subtree,df.2$subtree,df.3$subtree,df.4$subtree)))

A1.subtree.f.stat.0.10<-N
A1.subtree.f.stat.0.10.info<-Reduce(intersect,list(df.1$subtree,df.2$subtree,df.3$subtree,df.4$subtree))
# give rise to the averaged cell type composition of three differentiation samples
# Number of subtree: A1 9 G11 19 G2 7

# combine all samples
subtree.stereo<-list(A1.subtree.f.stat.0.10.info,G11.subtree.f.stat.0.10.info,G2.subtree.f.stat.0.10.info)
names(subtree.stereo)<-c("A1.CBRAD5","G11.CBRAD5","G2.CBRAD5")

dfTip.long <- readRDS("~/fig1/all_cbrad5_gs_hesc_tree_modify_tip_long_new.Rds")
all.subtree.depth<-mdepth.df%>%filter(grepl("A1",sample))  ## loop for samples: A1 G11 G2
colnames(all.subtree.depth)[1]<-"subtree"
df<-all.subtree.depth
df$mdepth.group<-df$mdepth/max(df$mdepth)
# For A1 CBRAD5 and G11 CBRAD5
zone<-seq(from=0,to=1.0,by=0.1)
# OR only for G2 CBRAD5
zone<-seq(from=0.1,to=1.0,by=0.1)
df.new.1<-lapply(seq((length(zone)-1)),function(x){
  min<-zone[x]
  max<-zone[x+1]
  index<-which(df$mdepth.group>min&df$mdepth.group<=max)
  if(length(index)<1){
    df.new<-as.data.frame(matrix(nrow=1,ncol=5))
    colnames(df.new)<-colnames(df)
  }else{
    df.new<-df[index,]
    df.new$mdepth.group<-max%>%as.numeric()}
  return(df.new)
})%>%bind_rows()%>%na.omit()
df.new.all<-df.new.1

all.subtree.depth<-df.new.all%>%dplyr::select("subtree","mdepth","mdepth.group","sample")
a1.all.subtree.depth<-all.subtree.depth
# combined all samples
mdepth.df.new<-rbind(a1.all.subtree.depth,g2.all.subtree.depth,g11.all.subtree.depth)

sample<-c("A1","G2","G11")

dfTip.long.depth.size<-mclapply(seq(sample),mc.cores=50,function(i){
  sample.df<-sample[i]
  subtree.info<-mdepth.df.new%>%filter(grepl(sample.df,sample))
  subtree.root<-subtree.info%>%select(subtree)
  size.info<-lapply(seq(subtree.root$subtree),function(f){
    subtree<-subtree.root[f,]
    depth<-subtree.info$mdepth.group[which(subtree.info$subtree==subtree)]
    subtree.df<-dfTip.long%>%filter(grepl(sample.df,intern))%>%filter(grepl(subtree,intern))
    subtree.df<-subtree.df[which(!grepl("add",subtree.df$intern)),]
    subtree.df$size<-length(subtree.df$tip)
    subtree.df$depth<-depth
    subtree.df$subtree<-subtree
    return(subtree.df)})%>%bind_rows()
  size.info<-size.info%>%dplyr::select("tip","celltype","size","depth","subtree")
  return(size.info)
})%>%bind_rows()

color.set<-readRDS("~/fig1/celltype_use_colors.Rds")

subtree.depth.size<-mclapply(seq(sample),mc.cores = 50,function(i){
  sample.df<-sample[i]
  subtree.this<-dfTip.long.depth.size%>%filter(grepl(sample.df,to))%>%dplyr::select(subtree,size,depth,celltype)
  subtree.info<-lapply(seq(unique(subtree.this$subtree)), function(f){
    df<-unique(subtree.this$subtree)[f]
    subtree.df<-subtree.this%>%filter(subtree==df)
    subtree.df.info<-data.frame(size=subtree.df$size%>%unique(),
                                depth=subtree.df$depth%>%unique(),
                                subtree=subtree.df$subtree%>%unique(),
                                sample=paste0(sample.df,".CBRAD5"))
    return(subtree.df.info)
    })%>%bind_rows()
  return(subtree.info)
  })%>%bind_rows()
 
subtree.celltype.pie<-mclapply(seq(sample),mc.cores = 50,function(i){
  sample.df<-sample[i]
  subtree.this<-dfTip.long.depth.size%>%filter(grepl(sample.df,to))%>%dplyr::select(subtree,size,depth,celltype)
  subtree.celltype.info<-lapply(seq(unique(subtree.this$subtree)), function(f){
    df<-unique(subtree.this$subtree)[f]
    subtree.df<-subtree.this%>%filter(subtree==df)
    celltype.info<-table(subtree.df$celltype)%>%as.data.frame()
    colnames(celltype.info)<-c("x","y")
    celltype.info$x<-as.character(celltype.info$x)
    color<-color.set[match(celltype.info$x,names(color.set))]
    pie <- ggplot(celltype.info, aes(x=1, y, fill=x)) + geom_bar(stat="identity", width=1) + coord_polar(theta="y") +
      theme_void() + theme(legend.position="none") + theme_transparent()+scale_fill_manual(values=color)
    return(pie)
    })
  names(subtree.celltype.info)<-unique(subtree.this$subtree)
  return(subtree.celltype.info)
  })
names(subtree.celltype.pie)<-paste0(sample,".CBRAD5")

# plot for A1.CBRAD5/G11.CBRAD5/G2.CBRAD5
subtree.info<-subtree.depth.size%>%filter(sample=="A1.CBRAD5")%>%filter(subtree%in%subtree.stereo$A1.CBRAD5)
pie.sub<-subtree.celltype.pie$A1.CBRAD5[names(subtree.celltype.pie$A1.CBRAD5)%in%subtree.stereo$A1.CBRAD5]
df<-tibble(x=subtree.info$size,
           y=subtree.info$depth,
           width=8,
           pie=pie.sub)
A1.df<-df # loop for samples: A1 G11 G2

# combine samples for plot
df.all<-rbind(A1.df,G11.df,G2.df)
df.all$x.log<-df.all$x%>%log2()
df.all$x<-df.all$x%>%log2()
df.all$width<-0.75
p <- ggplot(data=data.frame(x=c(0, 12), y=c(0, 1)), aes(x, y))+geom_blank()

pdf("~/fig4/pie_subtree_depth_size/cbrad5_stereo_subtree_pie.pdf",height=10,width=10)
p + geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=df.all)+
  theme_bw()+
  theme(plot.title = element_text(size=20, hjust = 0.5),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))+
  xlab("Sub-CLT size")+ylab("Normalized depth of its root")+
  scale_x_continuous(limits = c(0, 12), breaks = c(0, 2, 4, 6, 8, 10, 12), labels = c("0","4","16","64","256","1024","4096"))
dev.off()

#fig 4f
# simulate tree of 1024 tips with depth 10 and 1023 internal node
# three cell types :A B C D and cell type ratio is 1:2:2:4
# defined functional unit:1A-1B-2C-2D

# Three types of models for tree
# model 1: monoclonal origin of all cell types (as: A type ancestor produce all A type cells)
# model 2: randomized poly clonal origin of all cell types (as: all cell type was randomly distributed on the tree)
# model 3: stereotyped poly clonal origin of all cell types :each ancestor cell (total 128) could produce a sub clone with component of 2A-2B-4C  

# Test for robustness
# cell death rate: 1/1000-1/10, random ancestor cell ( 1023 internal node) dead according death rate
# calculate how many functional unit could be made up by the lived cells, and the number was divided by 128 to be defined as "functional capacity" 
# simulate 1000 times and mean the 1000 values of functional capacity to be as"development robustness"
# compare the robustness between 3 models, we expected model 3 > model 2 > model 1

##########  common step ###########
## input tree data
# tree file from the python script "Figure4_creat_full_balanced_bifurcated_phylogenetic_tree.py"
tree<-read.tree("~/fig4/full_balanced_bifurcated_phylogenetic_tree_1024.nwk")

#ggtree(tree)
death.rate<-c(0.001,0.005,0.01,0.05,0.1,0.5)

inSubtrees <- subtrees(tree)

all.subtree.real <- mclapply(seq(length(inSubtrees)), mc.cores=50,function(s){
  subtree <- inSubtrees[[s]]
  s.name <- subtree$name
  tip<-subtree$tip.label
  return(data.frame(subtree = s.name,
                    NodeName= tip,
                    stringsAsFactors = F))
}) %>% bind_rows()  

## Robustness calculate

DevRob.M3.random<-lapply(seq(length(death.rate)),function(i){
  treeTibble <- as_tibble(tree)
  inter.node<-treeTibble$node[treeTibble$label=="1"]
  deathRate<-death.rate[i]
  death.inter.counts<-length(inter.node)%>%`*`(deathRate)%>%round()
  treeTibble$celltype<-NA
  treeTibble$celltype[1:1024]<-Celltype.stero$celltype  # Celltype.mono;Celltype.random;Celltype.stero
  
  FunCap.Death<-mclapply(seq(1:1000),mc.cores=50,function(j){
    death.inter<-sample(x=inter.node,size=death.inter.counts,replace = FALSE)
    index.inter<-which(all.subtree.real$subtree%in%death.inter)
    d.tips<-all.subtree.real$NodeName[index.inter]%>%unique()
    index.tip<-which(treeTibble$label%in%d.tips)
    d.tree<-treeTibble[-index.tip,]
    
    live.tips<-d.tree%>%filter(.,celltype!="NA")
    fra.D<-length(which(live.tips$celltype=="D"))%>%`/`(4)
    fra.C<-length(which(live.tips$celltype=="C"))%>%`/`(2)
    fra.B<-length(which(live.tips$celltype=="B"))%>%`/`(1)
    fra.A<-length(which(live.tips$celltype=="A"))%>%`/`(1)
    FunCap<-min(fra.A,fra.B,fra.C,fra.D)%>%`/`(128)%>%as.data.frame()
    colnames(FunCap)<-"FunCap"
    return(FunCap)
  })%>%bind_rows()
  
  DevRob<-FunCap.Death$FunCap%>%as.vector()
  DevRob.data<-data.frame(Roubstness=DevRob,
                          death.rate=death.rate[i])
  return(DevRob.data)
})%>%bind_rows()

########### seperately step  ################

## tree mode 1 monoclonal tree

celltype=c(rep("A",128),rep("B",128),rep("C",256),rep("D",512))
# # or
# celltype=c(rep("A",64),rep("B",64),rep("C",128),rep("D",256),rep("E",512))

Celltype.mono <- as.character(celltype)%>%as.data.frame()
colnames(Celltype.mono)<-"celltype"

## tree mode 2 randomized poly clonal tree

celltype=c(rep("A",128),rep("B",128),rep("C",256),rep("D",512))
# # or
# celltype=c(rep("A",64),rep("B",64),rep("C",128),rep("D",256),rep("E",512))

Celltype.random <- sample(
  x = celltype,
  size = 1024,
  replace=FALSE)%>%as.data.frame()
colnames(Celltype.random)<-"celltype"
Celltype.random$celltype<-as.character(Celltype.random$celltype)

# tree mode 3 stereotype poly clonal tree
# step 1
# A : fixed stereotype
celltype=c("A","B","C","C","D","D","D","D")

Celltype<-sample(
  x = celltype,
  size = 8,
  replace=FALSE)%>%as.data.frame()

Celltype.stero<-rep(celltype,128)%>%as.data.frame()

# B: random stereotype
celltype=c("A","B","C","C","D","D","D","D")

Celltype.stero<-lapply(seq(c(1:128)), function(i){
  samplecell <- sample(
    x = celltype,
    size = 8,
    replace=FALSE)%>%as.data.frame()
})%>%bind_rows()

#step2
# stereotype assign to tips

colnames(Celltype.stero)<-"celltype"
Celltype.stero$celltype<-as.character(Celltype.stero$celltype)

## plot

colnames(DevRob.M1)<-c("Robustness","DeathRate")
DevRob.M1$type<-"Model_1"
colnames(DevRob.M2)<-c("Robustness","DeathRate")
DevRob.M2$type<-"Model_2"
colnames(DevRob.M3.fixed)<-c("Robustness","DeathRate")
DevRob.M3.fixed$type<-"Model_3.fixed"
colnames(DevRob.M3.random)<-c("Robustness","DeathRate")
DevRob.M3.random$type<-"Model_3.random"
DevRob.com<-rbind(DevRob.M1,DevRob.M2,DevRob.M3.fixed,DevRob.M3.random)
DevRob.com$DeathRate<-as.character(DevRob.com$DeathRate)

Data <-aggregate(Robustness ~ type + DeathRate, data = DevRob.com, FUN = function(x) c(mean = mean(x), sd=sd(x),se = std.error(x)))
group<-Data%>%select("type","DeathRate")
Mean<-Data$Robustness%>%as.data.frame()
Mean.data<-cbind(group,Mean)
Mean.data$log<-log(Mean.data$mean)

# try to remove model_1 when plot

DevRob.com.new<-DevRob.com%>%filter(!grepl("Model_1",type))

lvl = levels(as.factor(DevRob.com.new$type)) 

compare = expand.grid( a = lvl
                       , b = lvl
                       , stringsAsFactors = F) %>%
  filter( a != b) %>%
  mutate( comb = map2( a, b, function(x,y) sort( c(x,y) ) )
          , comb = map( comb, paste0, collapse = ',' )
  ) %>%
  unnest( comb ) #%>%
group_by( comb ) %>%
summarise() %>%
mutate( comb = stringr::str_split(comb, ',') ) %>%
.$comb

my_comparisons <- list( c("Model_2", "Model_3.fixed"), c("Model_2", "Model_3.random"), c("Model_3.fixed", "Model_3.random") )

colors <-c("#1E90FF","#FF8C00","#4DAF4A") #"#E41A1C",

pdf("~/fig4/development_robustness_model.signif.3.0.pdf",width =18,height = 10)
DevRob.com.new%>% 
  mutate(
    type = as.factor(type)
  ) %>%
  ggplot(aes(x = type, y = Robustness, color = type)) +
  facet_wrap(~DeathRate,scales = "free")+
  stat_boxplot(geom="boxplot",position = "dodge2",coef = 1.5)+ #outlier.shape = NA
  stat_summary(fun = "mean", geom = "point",size = 3, 
               alpha = .7, position = position_dodge(0.95)) +
  scale_color_manual(values = colors)+
  stat_compare_means(comparisons = my_comparisons, bracket.size = 0.3,label = "p.signif",method = "wilcox.test")+
  # or
  #geom_signif(comparisons = my_comparisons,
                 #map_signif_level=T, # T represent significance,F represent p value
                 #tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), # change the length of significance line
                 #textsize=10,# change the size of significant mark
                 #step_increase=0.2,
                 #y_position = c(1.2,1.6,2,2.4,2.8,3.2), # set the heigth of significant line
                 #size=1, # set the thickness of significant line
                 #test ="wilcox.test" )
  theme_bw()+
  theme(plot.title = element_text(size=20, hjust = 0.5),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        legend.text = element_text(size=20))+
  labs(x="Necrosis rate",y="Robustness",title = "wilcox test cell ratio:1:1:2:4")
dev.off()