library(SOT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scCustomize)
library(tidyverse)
library(purrr)
library(tidygraph)
library(ggraph)
library(cowplot)
library(grid)
library(igraph)
library(RColorBrewer)
library(stringi)
library(ggtree)
library(ape)
library(Seurat) # version 3.2.1
library(parallel)
library(cli)
library(GEOquery)
library(clusterProfiler)
library(org.Hs.eg.db)

# figure 1b,c

A1_CBRAD5_Seuratobj<-readRDS("~/fig1/Seuratobj/A1-CBRAD5_Seuratobj.Rds")
G11_CBRAD5_Seuratobj<-readRDS("~/fig1/Seuratobj/G11-CBRAD5_Seuratobj.Rds")
G2_CBRAD5_Seuratobj<-readRDS("~/fig1/Seuratobj/G2-CBRAD5_Seuratobj.Rds")
GS_HESC_Seuratobj<-readRDS("~/fig1/Seuratobj/GS-HESC_Seuratobj.Rds")

# allsample(diff+hesc) seurat object union
Seurat.list<-list(A1_CBRAD5_Seuratobj,G11_CBRAD5_Seuratobj,G2_CBRAD5_Seuratobj,GS_HESC_Seuratobj)
hesc.cbrad5_Seuratobj<-Reduce(function(x,y) merge(x,y) , Seurat.list)

### finding hvgs,otherwise,go to step 5

# 1. Normalize
hesc.cbrad5_Seuratobj <- NormalizeData(hesc.cbrad5_Seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)

## 2. The data manipulation of SOT
sce <- as.SingleCellExperiment(hesc.cbrad5_Seuratobj)
range(sce@assays@data$logcounts)

## Annotate Samples
sce$condition <- gsub(".*-", "", sce$orig.ident) 
metadata(sce)$`Sample color` = c("CBRAD5" = "#BEAED4","HESC" = "#386CB0")

## 3. Filter low expression gene
max(sce.raw@assays@data$logcounts["NKX2-1",])
sce = FilterLow(sce, minexp = 3, mincells = 10, datatype = "logcounts") ## minexp factor is key

## 4. Find high variable genes
sce = FindHVGs(sce, datatype = "logcounts", thr.bio = 0, thr.FDR = 0.1)
hvg.genes <- rowData(sce)$symbol[rowData(sce)$genes.use] %>% as.character()
#writeLines(hvg.genes,"~/fig1/hvg.genes.txt")
#hvg.genes<-readLines("~/fig1/hvg.genes.txt")

## 5. RUN PCA for seurat object by hvg.genes of SOT 
Seurat.wf <- function(Seuratobj,npcs,resolution){
  # 1. Normalize
  Seuratobj <- NormalizeData(Seuratobj, normalization.method = "LogNormalize", scale.factor = 10000)
  # 2. Find highly variable features
  Seuratobj <- FindVariableFeatures(
    object = Seuratobj,
    selection.method = "vst",
    num.bin = 20,
    mean.cutoff = c(0.1, 8),
    dispersion.cutoff = c(1, Inf)
  )
  # 3. cell cycle scoring
  Seuratobj <- CellCycleScoring(
    Seuratobj,
    s.features = cc.genes$s.genes,
    g2m.features = cc.genes$g2m.genes,
    set.ident = TRUE
  )
  # 4. scale data
  Seuratobj <- ScaleData(Seuratobj, features = rownames(Seuratobj))
  # 5. run PCA
  Seuratobj <-
    RunPCA(Seuratobj,
           features = hvg.genes,
           verbose = F,
           npcs = 100)
  # 6. cluster and run UMAP
  Seuratobj <- Seuratobj %>% 
    RunUMAP(reduction = "pca", dims = 1:npcs) %>% 
    FindNeighbors(dims = 1:npcs) %>% 
    FindClusters(resolution = resolution) %>% 
    identity()
  # 7. return results
  return(Seuratobj)
}


hesc.cbrad5_Seuratobj<-Seurat.wf(hesc.cbrad5_Seuratobj,npcs = 50, resolution = 0.6)

## 6. harmony
hesc.cbrad5_Seuratobj <- 
  RunHarmony(hesc.cbrad5_Seuratobj,
             group.by.vars = "orig.ident")
## 7. Umap by harmony
npcs <- 50
resolution <- 0.6
hesc.cbrad5_Seuratobj <- hesc.cbrad5_Seuratobj %>% 
  RunUMAP(reduction = "harmony", dims = 1:npcs) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:npcs) %>% 
  FindClusters(resolution = resolution) %>% 
  identity()

saveRDS(hesc.cbrad5_Seuratobj,"~/fig1/Seuratobj/hesc.cbrad5_Seuratobj.Rds")
HESC.CBRAD5_Seuratobj<-readRDS("~/fig1/Seuratobj/hesc.cbrad5_Seuratobj.Rds")
HESC.CBRAD5_Seuratobj@meta.data$seurat_clusters<-plyr::mapvalues(x = HESC.CBRAD5_Seuratobj@meta.data$seurat_clusters,
                                                                 from = c("3","0","7","6","10", "5","2","11","8","4","9","1"),
                                                                 to = c(paste0("C",seq(1,10)), "R1", "R2"))
Idents(HESC.CBRAD5_Seuratobj) <- HESC.CBRAD5_Seuratobj@meta.data$seurat_clusters
celltype.colors<-readRDS("~/fig1/celltype_use_colors.Rds")
names(celltype.colors)<-plyr::mapvalues(x = names(celltype.colors),
                                        from = c("3","0","7","6","10", "5","2","11","8","4","9","1"),
                                        to = c(paste0("C",seq(1,10)), "R1", "R2"))
DimPlot(HESC.CBRAD5_Seuratobj,group.by = "orig.ident",cols =c("#CAB2D6","#FDBF6F","#FF7F00","#B2DF8A"),pt.size = 1.5 ) # fig1b
FeaturePlot_scCustom(HESC.CBRAD5_Seuratobj, features = "NKX2-1",label = T,pt.size = 1) # fig1b

DimPlot(HESC.CBRAD5_Seuratobj,reduction = 'umap',label = T,cols = celltype.colors,pt.size = 1.5) #fig1c

# figure 1d

hesc_cbrad5_cluster_AVEexp<-read.csv("~/fig1/res0.6clusterAverages.csv")

##### GEO data analysis#####
GPL17930_ENTREZG<-read.table("~/fig1/GPL17930_hugene20st_Hs_ENTREZG_desc.annot.txt",
  header=T,sep="\t")

GSE83310_expr<-read.table("~/fig1/GSE83310_series_matrix_no_anno.txt",
  header = TRUE)
##### GEO data analysis#####

colnames(GPL17930_ENTREZG)[1]<-"ID_REF"
colnames(GPL17930_ENTREZG)[2]<-"probe_ID"
ENTREZID_SYMBLO<-read.table("~/fig1/mart_export.txt",header=T,sep="\t")
colnames(ENTREZID_SYMBLO)[3]<-"SYMBOL"
colnames(ENTREZID_SYMBLO)[2]<-"probe_ID"
colnames(ENTREZID_SYMBLO)[1]<-"ENSEMBL"
ID_REF_SYMBLO<-merge(GPL17930_ENTREZG,ENTREZID_SYMBLO,by="probe_ID")%>%`[`(-3)%>%`[`(-1)

sample_names<-colnames(GSE83310_expr)[-1]
probe_ids<-GSE83310_expr$ID_REF

GSE83310_expr<-merge(GSE83310_expr,ID_REF_SYMBLO,by="ID_REF")%>%`[`(-1)
GSE83310_expr <- GSE83310_expr[!duplicated(GSE83310_expr$SYMBOL), ]  ## GSE83310_expr <- GSE83310_expr %>% group_by(SYMBOL) %>% filter (! duplicated(SYMBOL))
rownames_GSE83310<-GSE83310_expr$SYMBOL
row.names(GSE83310_expr) <- as.character(rownames_GSE83310)
GSE83310_expr<-GSE83310_expr[, -which(colnames(GSE83310_expr)=="SYMBOL")]
GSE83310_expr<-GSE83310_expr[nchar(rownames(GSE83310_expr)) > 0, ]
class(rownames(GSE83310_expr))

# GSM2199282-84 neural GFP+ 
# GSM2199285-87 day 0 undifferentiated iPSC17 
# GSM2199288-90 day 3 definitive endoderm
# GSM2199291-93 day 6 anterior foregut endoderm 
# GSM2199294-96 day 15 Nkx2-1-GFP+
# GSM2199297-99 day 15 Nkx2-1-GFP-
# GSM2199300-302 day 28 Nkx2-1-GFP+
# GSM2199303-305 day 28 Nkx2-1-GFP-
# GSM2199306 differentiated AT2 cells
# GSM2199307 uncultured naïve lung epithelium 1
# GSM2199308 uncultured naïve lung epithelium 2
GSE83310_expr<-GSE83310_expr[,c(1:(which(colnames(GSE83310_expr)=="GSM2199300")-1))]
df.3<-GSE83310_expr[,which(colnames(GSE83310_expr)%in%c("GSM2199294","GSM2199295","GSM2199296"))]
df.2<-GSE83310_expr[,which(colnames(GSE83310_expr)%in%c("GSM2199297","GSM2199298","GSM2199299"))]
df.1<-GSE83310_expr[,setdiff(colnames(GSE83310_expr),c(colnames(df.2),colnames(df.3)))]
GSE83310_expr<-cbind(df.1,df.2,df.3)

GSE83310_SeuratObj <- CreateSeuratObject(
  counts = GSE83310_expr) 
GSE83310_exp<-as.matrix(GSE83310_SeuratObj@assays$RNA@data)
# Find highly variable features
GSE83310_SeuratObj <- FindVariableFeatures(
  object = GSE83310_SeuratObj,
  selection.method = "vst",
  num.bin = 20,
  mean.cutoff = c(0.1, 8),
  dispersion.cutoff = c(1, Inf)
)
# scale data
GSE83310_SeuratObj <- ScaleData(GSE83310_SeuratObj, features = rownames(GSE83310_SeuratObj))

GSE83310_SeuratObj@meta.data$times <-c("neural GFP+","neural GFP+","neural GFP+",
                                       "day 0 undifferentiated iPSC17","day 0 undifferentiated iPSC17","day 0 undifferentiated iPSC17",
                                       "day 3 definitive endoderm","day 3 definitive endoderm","day 3 definitive endoderm",
                                       "day 6 anterior foregut endoderm","day 6 anterior foregut endoderm","day 6 anterior foregut endoderm",
                                       "day 15 Nkx2-1-GFP-","day 15 Nkx2-1-GFP-","day 15 Nkx2-1-GFP-",
                                       "day 15 Nkx2-1-GFP+","day 15 Nkx2-1-GFP+","day 15 Nkx2-1-GFP+")
Idents(GSE83310_SeuratObj) <- GSE83310_SeuratObj@meta.data$times

GSE83310_SeuratObj.markers <- FindAllMarkers(object = GSE83310_SeuratObj, only.pos = TRUE, min.pct = 0.1, 
                                             thresh.use = 0.25,test.use = "bimod")

GSE83310_SeuratObj.markers.list<-lapply(seq(unique(GSE83310_SeuratObj.markers$cluster)),function(i){
  df<-unique(GSE83310_SeuratObj.markers$cluster)[i]
  df.data<-GSE83310_SeuratObj.markers%>%filter(cluster==df)%>%dplyr::select(gene)%>%as.data.frame()
  return(df.data)
})
names(GSE83310_SeuratObj.markers.list)<-unique(GSE83310_SeuratObj.markers$cluster)

colnames(sce@meta.data)
sce@meta.data<-sce@meta.data[,-c(18:93)]

sce_mutate_score<-mclapply(seq(GSE83310_SeuratObj.markers.list),mc.cores = 50,function(x){
  dataset<-GSE83310_SeuratObj.markers.list[[x]]
  gene<-dataset$gene
  score.name<-names(GSE83310_SeuratObj.markers.list)[x]%>%as.character()
  sce.new<- Seurat::AddModuleScore(object = sce,features = list(gene),name = paste(score.name,"1",sep=""))
  score<-sce.new@meta.data[,length(sce.new@meta.data)]%>%as.numeric()
  return(score)})%>%bind_cols()
colnames(sce_mutate_score)<-names(GSE83310_SeuratObj.markers.list)
rownames(sce_mutate_score)<-rownames(sce@meta.data)

sce@meta.data<-cbind(sce@meta.data,sce_mutate_score)

Features<-names(GSE83310_SeuratObj.markers.list)%>%as.character()
                              
Feature.dotplot<-DotPlot(sce,features = Features,cols = "RdYlBu",dot.scale=6)+RotatedAxis()+ # scaled the averaged expression of signal by dotplot
  ggplot2::coord_flip()+theme(axis.text=element_text(size=13))+scale_size_manual()

# figure 1e

cell.counts <- table(HESC.CBRAD5_Seuratobj$orig.ident, HESC.CBRAD5_Seuratobj$new.ident) %>% as.data.frame()

colnames(cell.counts) <- c("sample", "cluster", "cell.num")

cell.counts <- cell.counts %>% dplyr::group_by(sample) %>% do({
  df <- .
  df$ratio <- df$cell.num / sum(df$cell.num) * 100
  df
})

cell.counts$new.ident <- 
  plyr::mapvalues(cell.counts$cluster,
                  from=c("3","0","7","6","10", "5","2","11","8","4","9","1"),
                  to=c(paste0("C",seq(1,10)), "R1", "R2"))


new.ident.color <- colors.use
names(new.ident.color) <- plyr::mapvalues(names(new.ident.color),
                                          from=c("3","0","7","6","10", "5","2","11","8","4","9","1"),
                                          to=c(paste0("C",seq(1,10)), "R1", "R2"))

new.ident.color <- new.ident.color[c(paste0("C",seq(1,10)), "R1", "R2")]

celltype.ratio.plot <- 
  cell.counts %>% 
    ggplot(aes(x=sample, y=ratio, fill=cluster)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values=new.ident.color) +
    scale_y_continuous(limits = c(0,100), expand = c(0.01,0.01), 
                       breaks = seq(0,100, 25), labels = paste0(seq(0,100,25), "%"))+
    theme_bw() + labs(x="", y="", fill="celltype") +
    theme(axis.text.x = element_text(size=16, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size=16))
                              
# figure 1f

## ===================== define and source functions =========================
# require subfig.R and subplot.R from ggsubplot2 on github webset https://github.com/mkoohafkan/ggsubplot2/tree/master
# get tree tibble
get.node.depths <- function(info.df){
  all.inodes <- unique(info.df$from)
  all.inodes <- all.inodes[all.inodes != "0"]
  
  node.depth.df <- 
    lapply(all.inodes, function(n){
      d <- 1
      use.node <- n
      repeat{
        source.df <- filter(info.df, from %in% use.node)
        if (nrow(source.df) == 0) {
          break
        } 
        d <- d + 1
        use.node <- source.df$to %>% unlist() 
      }
      return(data.frame(to = n,
                        depth = d))
    }) %>% bind_rows()
  node.depth.df$to <- as.character(node.depth.df$to)
  return(node.depth.df)
}

count.ct <- function(ref.ct, test.ct){
  ref.ct <- ref.ct[!(ref.ct %in% c("inter", "root"))] %>% as.character()
  test.ct <- as.character(test.ct)
  ct.counts <- c()
  for (ct in ref.ct){
    ct.counts <- c(ct.counts, mean(test.ct == ct))
  }
  names(ct.counts) <- ref.ct
  ct.counts.df <- data.frame(as.list(ct.counts))
  colnames(ct.counts.df) <- ref.ct
  return(ct.counts.df)
}      

cp.plot.materials <- function(all.tree.info, sample.name = "all", max.depth = NULL, tip.tables, ct.colors = NULL, subfig.x.adj=-0.15, subfig.y.adj=-0.14){
  blank_theme <- theme(axis.line=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank(),
                       panel.background=element_blank(),
                       panel.border=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       plot.background=element_blank())
  
  ## extract sample tree infos
  if (sample.name == "all") {
    tree.info <- all.tree.info
  } else if (!(sample.name %in% all.tree.info$sample)) {
    stop("please use sample name within tree info dataframe.")
  } else {
    tree.info <- tree.info %>% filter(sample %in% sample.name)
  } 
  ## filter by max.depth 
  if (is.null(max.depth)) {
    tree.info <- tree.info
  } else {
    tree.info <- 
      tree.info %>% dplyr::group_by(sample) %>% do({
        s.df <- .
        s.max.depth <- max(s.df$depth, na.rm = T) - max.depth
        s.infos <- s.df %>% filter(., from.depth >= s.max.depth)
        s.infos
      })
  }
  ## extract edges and node infos to plot
  node.info <- list()
  i <- 1
  ref.cts <- unique(all.tree.info$celltype)
  newEdges <- 
    tree.info %>% 
    group_by(., from) %>% do({
      mydf <- .
      # check is there have terminal node 
      terminal.df <- mydf %>% filter(., !(to %in% tree.info$from))
      ## @condition 1: if no terminal node, then do not count celltype 
      if (nrow(terminal.df) == 0) {
        # add node infos to list
        info.df <- data.frame(name=mydf$from[1], plot.pie=FALSE, leaves=0)
        # return edges
        re.df <- mydf %>% select(from, to)
      } else {
        ## @condition 2: if have terminal node
        ## @ condition 2.1: terminal node all is leaves
        if (all(terminal.df$type == "leaf")) {
          # add node infos to list
          count.ct.df <- count.ct(ref.ct = ref.cts,
                                  test.ct = unlist(terminal.df$celltype))
          if (nrow(terminal.df) == 1){
            # return edges
            re.df <- mydf %>% select(from, to)
            # 
            info.df <- data.frame(name=c(terminal.df$to[1], mydf$from[1]),
                                  leaves = c(1, 0),
                                  plot.pie = c(TRUE, FALSE))
            info.df <- bind_cols(info.df, count.ct.df)
          } else {
            info.df <- count.ct.df %>% 
              mutate(
                plot.pie=TRUE,
                name=mydf$from[1],
                leaves= nrow(terminal.df)
              )
            # return edges
            re.df <- data.frame(NULL)
          }
          
        } else {
          ### @ condition 2.2: terminal node not all leaves
          # 1. get the from node infos
          # from's info
          from.info <- count.ct(ref.ct = ref.cts,
                                test.ct = "") %>% 
            mutate(name=mydf$from[1],
                   leaves=0,
                   plot.pie = FALSE)
          ## 2. count celltypes, 
          # when terminal is leaf, count own self
          # when terminal is internal node, count celltype ratio for each 
          # to's info
          info.df <- 
            lapply(seq(nrow(terminal.df)), function(t){
              
              sub.tm <- terminal.df[t,]
              #print(sub.tm)
              if (sub.tm$type == "leaf") {
                # add node infos to list
                count.ct.df <- count.ct(ref.ct = ref.cts,
                                        test.ct = unlist(sub.tm$celltype))
                sub.info.df <- count.ct.df %>% 
                  mutate(
                    plot.pie=TRUE,
                    name=sub.tm$to[1],
                    leaves = 1
                  )
              } else {
                node.leaves <- filter(tip.tables, intern %in% sub.tm$to[1]) %>% `[`("celltype") %>% unlist() %>% unique()
                count.ct.df <- count.ct(ref.ct = ref.cts,
                                        test.ct = node.leaves)
                sub.info.df <- count.ct.df %>% 
                  mutate(plot.pie=TRUE,
                         name=sub.tm$to[1],
                         leaves=length(node.leaves))
              }
              return(sub.info.df)
            }) %>% bind_rows()
          ## bind_all infos df
          info.df <- bind_rows(from.info, info.df)
          # return edges
          re.df <- mydf %>% select(from, to)
        }
      }
      node.info[[i]] <<- info.df
      i <<- i +1
      re.df
    })
  # bind all node infos to dataframe
  newNode <- data.table::rbindlist(node.info, fill = TRUE) %>% as.data.frame()
  newNode[is.na(newNode)] <- 0
  
  newEdges <- newEdges %>% filter(., from != 0)
  newNode <- newNode %>% filter(., name != 0)
  
  # set colors
  if (is.null(ct.colors)) {
    all.ct <- colnames(newNode)[!(colnames(newNode) %in% c("name", "plot.pie", "leaves"))] %>% as.character()
    all.ct <- all.ct[!is.na(all.ct) & (all.ct != "inter")] %>% as.character()
    ct.colors <- structure(colorRampPalette(brewer.pal(12, "Paired"))(length(all.ct)), names=all.ct)
  } else {
    ct.colors <- ct.colors
  }
  
  for (c in all.ct) {
    newNode[,c] <- as.character(newNode[,c]) %>% as.numeric()
  }
  newNode$name <- as.character(newNode$name)
  ### plot pie chart for each node
  ## plot pie 
  plot.list <- lapply(seq(nrow(newNode)), function(n){
    pie.df <- newNode %>% filter(name == newNode$name[n])
    if (pie.df$plot.pie == FALSE){
      return(ggplot() + 
               blank_theme + 
               theme(axis.text.x=element_blank()) + NoLegend() +
               theme(
                 rect = element_rect(fill = "transparent")
               ))
    } else {
      pie.df <- pie.df[, all.ct] ### TODO: change select celltype base on situation
      pie.df[is.na(pie.df)] <- 0
      pie.df <- pie.df %>% t() %>% as.data.frame() %>% set_names("ratio") %>%
        mutate(celltype = rownames(.))
      pie.df$ratio <- as.character(pie.df$ratio) %>% as.numeric()
      pie <- ggplot(pie.df)+
        geom_bar(aes(x="", y=ratio, fill=celltype), color=NA, width = 2, stat = "identity") + coord_polar("y", start=0) +
        blank_theme
      pie <- pie + scale_fill_manual(values = ct.colors) +
        theme(axis.text.x=element_blank(),
              legend.position = "none",
              validate = TRUE) +
        theme(
          rect = element_rect(fill = "transparent"),
          plot.margin = unit(c(0,0,0,0), "lines")
        )
      
      return(pie)
    }
  })
  
  ### create graph object
  cp_graph <- tbl_graph(edges = newEdges, nodes = newNode)
  cp_graph <- cp_graph %N>% 
    mutate(plot = plot.list)
  ## plot circle packing 
  cp_graph.cp.p <- 
    ggraph(cp_graph, layout = 'circlepack', direction="out", weight = leaves) + 
    geom_node_circle(aes(filter=!leaf)) +
    geom_node_circle(aes(filter=leaf)) +
    #geom_node_label(aes(label=name)) +
    geom_subfig(aes(x = x + subfig.x.adj, y = y + subfig.y.adj, plot = plot, width = 2.5*r, height = 2.5*r)) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "lines"))
  
  ### plot a test to get legend
  select.name <- sample(filter(newNode, plot.pie == TRUE) %>% select(name) %>% unlist(), 1) %>% as.character()
  test.df <- newNode %>% filter(name %in% select.name)
  test.df <- test.df[, all.ct] # TODO: change select celltype base on situation
  test.df[is.na(test.df)] <- 0
  test.df <- test.df %>% t() %>% as.data.frame() %>% set_names("ratio") %>%
    mutate(celltype = rownames(.))
  test.df$ratio <- as.character(test.df$ratio) %>% as.numeric()
  
  pie <- ggplot(test.df)+
    geom_bar(aes(x="", y=ratio, fill=celltype), width = 2, stat = "identity") + coord_polar("y", start=0) +
    blank_theme
  
  pie <- 
    pie + scale_fill_manual(values = ct.colors) +
    theme(axis.text.x=element_blank()) +
    theme(
      rect = element_rect(fill = "transparent")
    )
  
  legends = plot_grid(
    as_ggplot(get_legend(pie)))
  #### return results
  return(list(node.df = newNode,
              edge.df = newEdges,
              color.vector = ct.colors,
              pie.list = plot.list,
              graph.obj = cp_graph,
              circle.pack.p = cp_graph.cp.p,
              legend.p = legends
  ))
}


fine.adjust.cp.plots <- function(graph.obj, legend.plot, subfig.x.adj=-0.15, subfig.y.adj=-0.14, pie.r.scale = 2.5, rel.w=c(20,3)){
  # fine adjustment of circle packing plot
  cp_graph.cp.p <- 
    ggraph(graph.obj, layout = 'circlepack', direction="out", weight = leaves) + 
    geom_node_circle(aes(filter=!leaf)) +
    geom_node_circle(aes(filter=leaf)) +
    #geom_node_label(aes(label=name)) +
    geom_subfig(aes(x = x + subfig.x.adj, y = y + subfig.y.adj, plot = plot, width = pie.r.scale*r, height = pie.r.scale*r)) +
    theme_void() +
    theme(plot.margin = unit(c(0,0,0,0), "lines"))
  # re-combine circle packing plot and legend
  combine_plots <- 
    plot_grid(
      cp_graph.cp.p + theme(plot.margin = unit(c(0,0,0,0), "lines"),
                            rect = element_rect(fill = "transparent")), 
      legend.plot + theme(plot.margin = unit(c(0,0,0,0), "lines"),
                          rect = element_rect(fill = "transparent")),
      rel_widths = rel.w, nrow = 1, ncol = 2, align = "v")
  return(combine_plots)
}

## ================================================ load and save vars ====================

tree.info <- readRDS("~/fig1/all_cbrad5_GS_hesc_tree_dataframe_modify_new.Rds")
dfTip.long <- readRDS("~/fig1/all_cbrad5_gs_hesc_tree_modify_tip_long_new.Rds")

all.ct <- tree.info$celltype %>% unique() %>% as.character()
all.ct <- all.ct[!is.na(all.ct) & (all.ct != "inter")] %>% as.character()
ct.colors <- structure(colorRampPalette(brewer.pal(12, "Paired"))(length(all.ct)), names=all.ct)
saveRDS(ct.colors,file = "~/fig1/celltype_use_colors.Rds")
## =============================== plot =======================================================
#### combine legend and graph
all.sample.plot <- 
  cp.plot.materials(all.tree.info = tree.info,
                    sample.name = "all", 
                    max.depth = 12, 
                    tip.tables = dfTip.long)

### G11-CBRAD5
g11.sample.plot <- 
  cp.plot.materials(all.tree.info = tree.info,
                    sample.name = "G11-CBRAD5", 
                    max.depth = 15, 
                    tip.tables = dfTip.long,
                    ct.colors = ct.colors,
                    subfig.x.adj=-0.15, 
                    subfig.y.adj=-0.14)


g11cbrad5.combine_plots <- 
  fine.adjust.cp.plots(graph.obj = g11.sample.plot$graph.obj,
                       legend.plot = g11.sample.plot$legend.p,
                       subfig.x.adj=-0.15, 
                       subfig.y.adj=-0.14, 
                       pie.r.scale = 2.5,
                       rel.w=c(20,3))

g11cbrad5.combine_plots

## A1-CBRAD5
a1.sample.plot <- 
  cp.plot.materials(all.tree.info = tree.info,
                    sample.name = "A1-CBRAD5", 
                    max.depth = 15, 
                    tip.tables = dfTip.long,
                    ct.colors = ct.colors,
                    subfig.x.adj=-0.15, 
                    subfig.y.adj=-0.14)

a1cbrad5.combine_plots <- 
  fine.adjust.cp.plots(graph.obj = a1.sample.plot$graph.obj,
                       legend.plot = a1.sample.plot$legend.p,
                       subfig.x.adj=-0.15, 
                       subfig.y.adj=-0.14, 
                       pie.r.scale = 2.5,
                       rel.w=c(20,3))

a1cbrad5.combine_plots

## G2-CBRAD5
g2.sample.plot <- 
  cp.plot.materials(all.tree.info = tree.info,
                    sample.name = "G2-CBRAD5", 
                    max.depth = 15, 
                    tip.tables = dfTip.long,
                    ct.colors = ct.colors,
                    subfig.x.adj=-0.14, 
                    subfig.y.adj=-0.13)


g2cbrad5.combine_plots <- 
  fine.adjust.cp.plots(graph.obj = g2.sample.plot$graph.obj,
                       legend.plot = g2.sample.plot$legend.p,
                       subfig.x.adj=-0.15, 
                       subfig.y.adj=-0.14, 
                       pie.r.scale = 2.5,
                       rel.w=c(20,3))

g2cbrad5.combine_plots

## GS-HESC
gs.sample.plot <- 
  cp.plot.materials(all.tree.info = tree.info,
                    sample.name = "GS-HESC", 
                    max.depth = 15, 
                    tip.tables = dfTip.long,
                    ct.colors = ct.colors,
                    subfig.x.adj=-0.14, 
                    subfig.y.adj=-0.13)

gshesc.combine_plots <- 
  fine.adjust.cp.plots(graph.obj = gs.sample.plot$graph.obj,
                       legend.plot = gs.sample.plot$legend.p,
                       subfig.x.adj=-0.15, 
                       subfig.y.adj=-0.14, 
                       pie.r.scale = 2.5,
                       rel.w=c(20,3))

gshesc.combine_plots

# figure 1g

subtree.depths <- function(phylo.object){
  tree.tibble <- as_tibble(phylo.object) %>% as.data.frame()
  all.subtrees <- subtrees(phylo.object)
  distMatrix <- dist.nodes(phylo.object)
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

get.adistVStdist.dist <- function(allele.info_path, tree_path){
  allele.info <- read.csv(allele.info_path, sep = '\t', stringsAsFactors = F)
  node.edits.df <- allele.info %>% select(NodeName, !!!syms(c(paste0("target", seq(13)))))
  node.edits.df <- node.edits.df[!duplicated(node.edits.df),]
  rownames(node.edits.df) <- node.edits.df$NodeName
  node.edits.df <- node.edits.df %>% dplyr::select(-NodeName)
  node.edits <- node.edits.df %>% as.matrix()
  # get node compare df
  node.comp.df <- t(combn(rownames(node.edits), 2)) %>% as.data.frame(., stringsAsFactors =F) 
  colnames(node.comp.df) <- c("node1", "node2")
  # get tree depth as distance
  tree.obj <- read.tree(tree_path)
  distMatrix <- dist.nodes(tree.obj)
  treeTibble <- as_tibble(tree.obj) %>% as.data.frame()
  all.subtree.depth <- subtree.depths(tree.obj)
  dist.comp <- 
    mclapply(seq(nrow(node.comp.df)), 
             mc.cores = 60,
             function(n){
               n1 <- node.comp.df$node1[n]
               n2 <- node.comp.df$node2[n]
               # extract node1 and node2 edits
               n1.edit <- node.edits[n1,] %>% as.vector()
               n2.edit <- node.edits[n2,] %>% as.vector()
               # compare allelic distance between these two nodes
               ad <- 0
               for (ie in seq(length(n1.edit))) {
                 if ((n1.edit[ie] != n2.edit[ie]) & (n1.edit[ie] != "NONE") & (n2.edit[ie] != "NONE")){
                   ad <- ad + 2
                 } else if ((n1.edit[ie] == "NONE" | n2.edit[ie] == "NONE") & (n1.edit[ie] != n2.edit[ie])) {
                   ad <- ad + 1
                 } else {
                   ad <- ad + 0
                 }
               }
               mrca <- getMRCA(tree.obj, c(n1, n2))
               # compare tree dist (use the MRCA's depth as distance for two node)
               mrca.depth <- filter(all.subtree.depth, subtree.root == mrca) %>% `[`("depth") %>% unlist()
               # or tree dist used the depth of max distrance for pairs to most recent common ancestors 
               n1.seqName <- filter(treeTibble, label == n1) %>% `[`("node") %>% unlist()
               n2.seqName <- filter(treeTibble, label == n2) %>% `[`("node") %>% unlist()
               treedist <- max(distMatrix[n1.seqName, mrca], distMatrix[n2.seqName, mrca])
               # return results
               return(data.frame(node1 = n1,
                                 node2 = n2,
                                 allelic.dist = ad,
                                 tree.dist = treedist,
                                 stringsAsFactors = F))
             }) %>% bind_rows()
  
  dist.comp$nor.allelic.dist <- dist.comp$allelic.dist / 26
  dist.comp$nor.tree.dist <- dist.comp$tree.dist / max(dist.comp$tree.dist)
  dist.comp$node1 <- as.character(dist.comp$node1)
  dist.comp$node2 <- as.character(dist.comp$node2)
  return(dist.comp)
}

### ================================= figure 1e tree distance VS allelic barcode distance =========================
####### G11-CBRAD5
g11.cbrad5.atdist.df <- 
  get.adistVStdist.dist(allele.info_path = "~/fig1/G11-CBRAD5.AllelesInfo.txt",
                        tree_path = "~/fig1/G11-CBRAD5.nwk")

cor.test(g11.cbrad5.atdist.df$allelic.dist, g11.cbrad5.atdist.df$tree.dist, method = "spearman")
cor.test(g11.cbrad5.atdist.df$nor.allelic.dist, g11.cbrad5.atdist.df$nor.tree.dist, method = "pearson")

g11.cbrad5.alDvsTrD.denbin.p <- 
  g11.cbrad5.atdist.df %>% 
  ggplot(aes(x=nor.allelic.dist, y=nor.tree.dist)) +
  #geom_tile() +
  geom_bin2d(aes(fill=..density..), bins = 10) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x="Normalized allelic distance", y="Normalized tree distance", title = "G11-CBRAD5\npearson cor:0.5067688 , pvalue < 2.2e-16")

####### A1-CBRAD5
a1.cbrad5.atdist.df <- 
  get.adistVStdist.dist(allele.info_path = "~/fig1/A1-CBRAD5.AllelesInfo.txt",
                        tree_path = "~/fig1/A1-CBRAD5.nwk")

cor.test(a1.cbrad5.atdist.df$allelic.dist, a1.cbrad5.atdist.df$tree.dist, method = "spearman")
cor.test(a1.cbrad5.atdist.df$nor.allelic.dist, a1.cbrad5.atdist.df$nor.tree.dist, method = "pearson")

a1.cbrad5.alDvsTrD.denbin.p <- 
  a1.cbrad5.atdist.df %>% 
  ggplot(aes(x=nor.allelic.dist, y=nor.tree.dist)) +
  #geom_tile() +
  geom_bin2d(aes(fill=..density..), bins = 10) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x="Normalized allelic distance", y="Normalized tree distance", title = "A1-CBRAD5\npearson cor:0.4357036 , pvalue < 2.2e-16")

###### G2-CBRAD5
g2.cbrad5.atdist.df <- 
  get.adistVStdist.dist(allele.info_path = "~/fig1/G2-CBRAD5.AllelesInfo.txt",
                        tree_path = "~/fig1/G2-CBRAD5.nwk")

cor.test(g2.cbrad5.atdist.df$allelic.dist, g2.cbrad5.atdist.df$tree.dist, method = "spearman")
cor.test(g2.cbrad5.atdist.df$nor.allelic.dist, g2.cbrad5.atdist.df$nor.tree.dist, method = "pearson")

g2.cbrad5.alDvsTrD.denbin.p <- 
  g2.cbrad5.atdist.df %>% 
  ggplot(aes(x=nor.allelic.dist, y=nor.tree.dist)) +
  #geom_tile() +
  geom_bin2d(aes(fill=..density..), bins = 10) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  labs(x="Normalized allelic distance", y="Normalized tree distance", title = "G2-CBRAD5\npearson cor:0.4196391, pvalue < 2.2e-16")

###### GS-HESC
gs.hesc.atdist.df <- 
  get.adistVStdist.dist(allele.info_path = "~/fig1/GS-HESC.AllelesInfo.txt",
                        tree_path = "~/fig1/GS-HESC.nwk")

cor.test(gs.hesc.atdist.df$allelic.dist, gs.hesc.atdist.df$tree.dist, method = "spearman")
cor.test(gs.hesc.atdist.df$nor.allelic.dist, gs.hesc.atdist.df$nor.tree.dist, method = "pearson")

### all 
all.atdist.df <- bind_rows(a1.cbrad5.atdist.df, g11.cbrad5.atdist.df, g2.cbrad5.atdist.df, gs.hesc.atdist.df)

cor.test(all.atdist.df$nor.allelic.dist, all.atdist.df$nor.tree.dist, method = "spearman") # 0.5511992 p-value < 2.2e-16
cor.test(all.atdist.df$nor.allelic.dist, all.atdist.df$nor.tree.dist, method = "pearson") # 0.6057589 p-value < 2.2e-16

all.atdist.df.0.2<-all.atdist.df[all.atdist.df$nor.tree.dist<=0.2,]
all.atdist.df.0.4<-all.atdist.df[all.atdist.df$nor.tree.dist<=0.4&all.atdist.df$nor.tree.dist>0.2,]
all.atdist.df.0.6<-all.atdist.df[all.atdist.df$nor.tree.dist<=0.6&all.atdist.df$nor.tree.dist>0.4,]
all.atdist.df.0.8<-all.atdist.df[all.atdist.df$nor.tree.dist<=0.8&all.atdist.df$nor.tree.dist>0.6,]
all.atdist.df.1.0<-all.atdist.df[all.atdist.df$nor.tree.dist<=1&all.atdist.df$nor.tree.dist>0.8,]

all.atdist.df.0.2$nor.tree.dist.group<-"1"
all.atdist.df.0.4$nor.tree.dist.group<-"2"
all.atdist.df.0.6$nor.tree.dist.group<-"3"
all.atdist.df.0.8$nor.tree.dist.group<-"4"
all.atdist.df.1.0$nor.tree.dist.group<-"5"

new.all.atdist.df<-rbind(all.atdist.df.0.2,all.atdist.df.0.4,all.atdist.df.0.6,all.atdist.df.0.8,all.atdist.df.1.0)

#install.packages('ggrastr')
new.all.atdist.df$nor.tree.dist.group<-new.all.atdist.df$nor.tree.dist.group%>%as.factor()

new.all.atdist.df %>%
  ggplot(aes(x=nor.tree.dist.group, y=nor.allelic.dist),fill=nor.tree.dist.group) +
  geom_boxplot(fill=c("#A6BDDB","#74A9CF","#3690C0","#0570B0","#045A8D"), #ggrastr::rasterise raster.dpi = 2400
                                  color="#D0D1E6")+
  #geom_boxplot()+
  stat_summary(fun.y=mean, geom="point", shape=23, size=8,fill="white",color="#D0D1E6")+
  theme_bw()+
  labs(x="normalized topological distance",
       y="normalized allelic distance",
       title="All CBRAD5 and GS-HESC samples,Pearson's cor=0.6057589 p-value < 2.2e-16")+
  theme(plot.title = element_text(size=20, hjust = 0.5),
        panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20))

dev.off()

# figure 1h

####               some functions               ####
#####################################################
checkCon2 <- function(Node1, Node2, alleleInfo, rawEvents){
  con1 <- filter(alleleInfo, NodeName == Node1)
  con2 <- filter(alleleInfo, NodeName == Node2)
  # get cell BC
  BC1 <- con1$BC[1]
  BC2 <- con2$BC[1]
  # get raw events
  raw1 <- filter(rawEvents, BC == BC1)[,c(paste0("target", seq(1,13)))]
  raw2 <- filter(rawEvents, BC == BC2)[,c(paste0("target", seq(1,13)))]
  # overlap edits
  overlap <- dplyr::intersect(raw1, raw2)
  return(data.frame(node1 = Node1,
                    node2 = Node2,
                    overlapE = nrow(overlap), stringsAsFactors = F))
}

getOverlapInfo.dist <- function(Pair.df, alleleInfo, rawEvents){
  OverlapInfo <- mclapply(seq(nrow(Pair.df)), function(i){
    n1 <- Pair.df$node1[i]
    n2 <- Pair.df$node2[i]
    td <- Pair.df$treeDist[i]
    subCheck <- checkCon2(n1, n2, alleleInfo, rawEvents)
    subCheck$treeDist <- td
    return(subCheck)
  }, mc.cores = 50) %>% bind_rows(.)
  return(OverlapInfo)
}

filter.ANC.ratio <- function(moreInfo, alleleInfo, ancInfo, cutoff){
  # 1. add event tag in moreInfo
  moreInfo$eventTag <- with(moreInfo, paste(target1, target2, target3, target4, target5, target6, 
                                            target7, target8, target9, target10, target11, target12, 
                                            target13, sep = "_"))
  # 2. calculate the anc event ratio in cells
  anc.moreInfo <- filter(moreInfo, checkInfo == "anc")
  anc.moreInfo <- anc.moreInfo[, c("BC", "eventTag")] %>% unique()
  anc.edit.freq <- anc.moreInfo %>% dplyr::group_by(eventTag) %>% 
    dplyr::summarise(cell.num=n())
  anc.edit.freq$freq <- anc.edit.freq$cell.num / sum(anc.edit.freq$cell.num)
  # 3. filter eventTag freq bigger than cutoff
  filter.event.df <- filter(anc.edit.freq, freq >= cutoff)
  if (nrow(filter.event.df) == 0) {
    message("no anc event frequency >= ", cutoff)
    fix.df <- filter(moreInfo, checkInfo != "dis")
  } else {
    filter.event <- filter.event.df$eventTag %>% unlist()
    print(filter.event)
    fix.df <- filter(moreInfo, checkInfo != "dis" & !(eventTag %in% filter.event))
  }
  # 4. calculate the share ANC ratio when filter out frequency more than cutoff
  node.pair.df <- filter(ancInfo$overlapInfo_dist, node1 %in% ancInfo$haveANC_oneCellNode & node2 %in% ancInfo$haveANC_oneCellNode)
  fix.overlapInfo.df <- getOverlapInfo.dist(node.pair.df, alleleInfo, fix.df)
  fix.haveANC.share_ratio <- fix.overlapInfo.df %>% dplyr::group_by(treeDist) %>% dplyr::summarise(share_ratio=sum(overlapE > 0) / n())
  return(fix.haveANC.share_ratio)
}

### filter out the most frequency 

###################
### calculate the exp level of the have anc cells
exp.level.compare <- function(fix.moreInfo, alleleInfo_path, anc.nodes){
  alleleInfo <- read.csv(alleleInfo_path, sep = '\t', stringsAsFactors = F)
  ## 1. use all oneCell node 
  oneCellNode.BC <- filter(alleleInfo, cellNum == 1) %>% `[`("BC") %>% unlist()
  anc.BC <- filter(alleleInfo, NodeName %in% anc.nodes) %>% `[`("BC") %>% unlist()
  ## 2. extract counts matrix
  BC.umiCounts <- fix.moreInfo %>% filter(checkInfo != "dis") %>% group_by(BC) %>% summarise(umi.num=n())
  oneCellNode.umis <- BC.umiCounts %>% filter(BC %in% oneCellNode.BC)
  oneCellNode.umis$haveANC <- sapply(oneCellNode.umis$BC, function(x) ifelse(x %in% anc.BC, "Y", "N"))
  return(oneCellNode.umis)
}

######################
### plot functions
plot.ShareANC <- function(shareANC.info, plot.title){
  ggplot(shareANC.info, aes(x=treeDist, y = share_ratio)) +
    geom_bar(stat = "identity") +
    scale_x_continuous(breaks = unique(shareANC.info)) +
    labs(x="Tree distance", y= "Share ancestral barcode ratio", title = plot.title) +
    theme(plot.title = element_text(size=20, hjust = 0.5),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size=20),
          axis.title = element_text(size=20))
}

plot.exp.level <- function(level.df, sample.name){
  level.df %>% ggplot(aes(fill=haveANC)) +
    geom_histogram(aes(x=newCount, y=..count../sum(..count..)), alpha=0.5, binwidth = 1, position="identity") + 
    scale_x_continuous(breaks = seq(0,50, 10), labels = c(seq(0, 40, 10), ">50")) +
    scale_fill_discrete(labels = c("NO", "YES")) +
    labs(x= paste0("Expression level of ", sample.name, " in each cells"), y="Frequency", fill="Cells with ancestral\ntranscripts") +
    theme(plot.title = element_text(size=20, hjust = 0.5),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size=20),
          axis.title = element_text(size=20),
          legend.position = c(0.7, 0.7))
  
}

filter.ANC.most <- function(moreInfo_path, alleleInfo_path, ancInfo_path){
  moreInfo <- read.csv(moreInfo_path, sep = '\t', stringsAsFactors = F)
  alleleInfo <- read.csv(alleleInfo_path, sep = '\t', stringsAsFactors = F)
  ancInfo <- readRDS(ancInfo_path)
  # 1. add event tag in moreInfo
  moreInfo$eventTag <- with(moreInfo, paste(target1, target2, target3, target4, target5, target6, 
                                            target7, target8, target9, target10, target11, target12, 
                                            target13, sep = "_"))
  # 2. calculate the anc event ratio in cells
  anc.moreInfo <- filter(moreInfo, checkInfo == "anc")
  anc.moreInfo <- anc.moreInfo[, c("BC", "eventTag")] %>% unique()
  anc.edit.freq <- anc.moreInfo %>% dplyr::group_by(eventTag) %>% 
    dplyr::summarise(cell.num=n())
  anc.edit.freq$freq <- anc.edit.freq$cell.num / sum(anc.edit.freq$cell.num)
  cutoff <- max(anc.edit.freq$freq)
  # 3. filter eventTag freq bigger than cutoff
  filter.event.df <- filter(anc.edit.freq, freq >= cutoff)
  if (nrow(filter.event.df) == 0) {
    message("no anc event frequency >= ", cutoff)
    fix.df <- filter(moreInfo, checkInfo != "dis")
  } else {
    filter.event <- filter.event.df$eventTag %>% unlist()
    print(filter.event)
    fix.df <- filter(moreInfo, checkInfo != "dis" & !(eventTag %in% filter.event))
  }
  # 4. calculate the share ANC ratio when filter out frequency more than cutoff
  fix.ANC_BC <- filter(fix.df, checkInfo == "anc") %>% `[`("BC") %>% unlist()
  filter.node <- filter(alleleInfo, BC %in% fix.ANC_BC) %>% `[`("NodeName") %>% unique() %>% unlist()
  node.pair.df <- filter(ancInfo$overlapInfo_dist, node1 %in% filter.node & node2 %in% filter.node)
  ## 4.1 bootstrap
  boot.times <- 100
  boot.ratio <- 0.9
  
  # ==== return node and BC
  final_nodes <- unique(c(node.pair.df$node1, node.pair.df$node2))
  fix.overlapInfo.df <- getOverlapInfo.dist(node.pair.df, alleleInfo, fix.df)
  ## bootstrap 
  boot.re.replace <- lapply(seq(boot.times), function(i){
    boot.fix.overlapInfo.df <- sample_n(fix.overlapInfo.df, size = nrow(fix.overlapInfo.df), replace = T)
    boot.fix.haveANC.share_ratio <- boot.fix.overlapInfo.df %>% dplyr::group_by(treeDist) %>% dplyr::summarise(share_ratio=sum(overlapE > 0) / n())
    boot.fix.haveANC.share_ratio$boot.time <- i
    #print(i)
    return(boot.fix.haveANC.share_ratio)
  }) %>% bind_rows()
  
  boot.re.replace.summary <- boot.re.replace %>% group_by(treeDist) %>% summarise(boot.share_ratio.mean = mean(share_ratio),
                                                                                  boot.share_ratio.sd = sd(share_ratio))
  ### output results
  fix.haveANC.share_ratio <- fix.overlapInfo.df %>% dplyr::group_by(treeDist) %>% dplyr::summarise(share_ratio=sum(overlapE > 0) / n())
  fix.haveANC.share_ratio <- left_join(fix.haveANC.share_ratio, boot.re.replace.summary)
  return(list(ANC.ratio.df = fix.haveANC.share_ratio,
              boot.ANC.ratio.df = boot.re.replace,
              haveANC.nodes = final_nodes,
              fix.moreInfo = fix.df,
              most.ratio.ANC.event = filter.event))
}

plot.ShareANC.errorbar <- function(shareANC.info, plot.title){
  ggplot(shareANC.info, aes(x=depth.group, y = mean_share_ratio)) + # treeDist
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin=mean_share_ratio-mean_boot.share_ratio.sd, ymax=mean_share_ratio+mean_boot.share_ratio.sd), width=0.5) +
    #scale_x_continuous(breaks = unique(shareANC.info$depth.group)) +
    labs(x="Tree distance", y= "Share ancestral barcode ratio", title = plot.title) +
    theme(plot.title = element_text(size=20, hjust = 0.5),
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(size=20),
          axis.title = element_text(size=20))
}

# A1-CBRAD5
A1.CBRAD5.ANC.most <- 
  filter.ANC.most(moreInfo_path = "~/fig1/A1-CBRAD5_moreInfos.txt",
                  alleleInfo_path = "~/fig1/A1-CBRAD5.AllelesInfo.txt",
                  ancInfo_path = "~/fig1/A1-CBRAD5_ANC_share_analysis_results.Rds")

A1.CBRAD5.exp.level <- 
  exp.level.compare(fix.moreInfo = A1.CBRAD5.ANC.most$fix.moreInfo,
                    alleleInfo_path = "~/fig1/A1-CBRAD5.AllelesInfo.txt",
                    anc.nodes = A1.CBRAD5.ANC.most$haveANC.nodes)

plot.ShareANC.errorbar(A1.CBRAD5.ANC.most$ANC.ratio.df, "A1-CBRAD5")

A1.CBRAD5.exp.level$newCount <- sapply(A1.CBRAD5.exp.level$umi.num, function(x) ifelse(x>=50, 50, x))
plot.exp.level(A1.CBRAD5.exp.level, "A1-CBRAD5")

# G11-CBRAD5
G11.CBRAD5.ANC.most <- 
  filter.ANC.most(moreInfo_path = "~/fig1/G11-CBRAD5_moreInfos.txt",
                  alleleInfo_path = "~/fig1/G11-CBRAD5.AllelesInfo.txt",
                  ancInfo_path = "~/fig1/G11-CBRAD5_ANC_share_analysis_results.Rds")

G11.CBRAD5.exp.level <- 
  exp.level.compare(fix.moreInfo = G11.CBRAD5.ANC.most$fix.moreInfo,
                    alleleInfo_path = "~/fig1/G11-CBRAD5.AllelesInfo.txt",
                    anc.nodes = G11.CBRAD5.ANC.most$haveANC.nodes)

plot.ShareANC.errorbar(G11.CBRAD5.ANC.most$ANC.ratio.df, "G11-CBRAD5")
G11.CBRAD5.exp.level$newCount <- sapply(G11.CBRAD5.exp.level$umi.num, function(x) ifelse(x>=50, 50, x))
plot.exp.level(G11.CBRAD5.exp.level, "G11-CBRAD5")
# G2-CBRAD5
G2.CBRAD5.ANC.most <- 
  filter.ANC.most(moreInfo_path = "~/fig1/G2-CBRAD5_moreInfos.txt",
                  alleleInfo_path = "~/fig1/G2-CBRAD5.AllelesInfo.txt",
                  ancInfo_path = "~/fig1/G2-CBRAD5_ANC_share_analysis_results.Rds")

G2.CBRAD5.exp.level <- 
  exp.level.compare(fix.moreInfo = G2.CBRAD5.ANC.most$fix.moreInfo,
                    alleleInfo_path = "~/fig1/G2-CBRAD5.AllelesInfo.txt",
                    anc.nodes = G2.CBRAD5.ANC.most$haveANC.nodes)

plot.ShareANC.errorbar(G2.CBRAD5.ANC.most$ANC.ratio.df, "G2-CBRAD5")
G2.CBRAD5.exp.level$newCount <- sapply(G2.CBRAD5.exp.level$umi.num, function(x) ifelse(x>=50, 50, x))
plot.exp.level(G2.CBRAD5.exp.level, "G2-CBRAD5")

# sample as A1 
group.A1<-seq(from=1,to=11,by=1)%>%as.data.frame()
group.A1$depth.rel<-group.A1[,1]/11
group.A1.1<-group.A1%>%filter(depth.rel<=0.2)%>%mutate(group.depth=1)
group.A1.2<-group.A1%>%filter(depth.rel>0.2&depth.rel<=0.4)%>%mutate(group.depth=2)
group.A1.3<-group.A1%>%filter(depth.rel>0.4&depth.rel<=0.6)%>%mutate(group.depth=3)
group.A1.4<-group.A1%>%filter(depth.rel>0.6&depth.rel<=0.8)%>%mutate(group.depth=4)
group.A1.5<-group.A1%>%filter(depth.rel>0.8&depth.rel<=1.0)%>%mutate(group.depth=5)

group.A1<-rbind(group.A1.1,group.A1.2,group.A1.3,group.A1.4,group.A1.5)
colnames(group.A1)<-c("depth","depth.rel","group.depth")

A1.CBRAD5.ANC.most$ANC.ratio.df$depth.group<-group.A1$group.depth[match(A1.CBRAD5.ANC.most$ANC.ratio.df$treeDist,group.A1$depth)]

# combine all samples
cbrad5.ANC.most.new<-rbind(A1.CBRAD5.ANC.most$ANC.ratio.df,G11.CBRAD5.ANC.most$ANC.ratio.df,G2.CBRAD5.ANC.most$ANC.ratio.df)

cbrad5.ANC.most.new <- cbrad5.ANC.most.new %>% group_by(depth.group) %>% 
  summarise(mean_share_ratio=mean(share_ratio),
            mean_boot.share_ratio.mean= mean(boot.share_ratio.mean),
            mean_boot.share_ratio.sd=mean(boot.share_ratio.sd)) %>%
  as.data.frame()

plot.ShareANC.errorbar(cbrad5.ANC.most.new, "All CBRAD5")

#========== HESC ====

# GS
GS.HESC.ANC.most <- 
  filter.ANC.most(moreInfo_path = "~/fig1/GS-HESC_moreInfos.txt",
                  alleleInfo_path = "~/fig1/GS-HESC.AllelesInfo.txt",
                  ancInfo_path = "~/fig1/GS-HESC_ANC_share_analysis_results.Rds")

GS.HESC.exp.level <- 
  exp.level.compare(fix.moreInfo = GS.HESC.ANC.most$fix.moreInfo,
                    alleleInfo_path = "~/fig1/GS-HESC.AllelesInfo.txt",
                    anc.nodes = GS.HESC.ANC.most$haveANC.nodes)

group.GS<-seq(from=1,to=12,by=1)%>%as.data.frame()
group.GS$depth.rel<-group.GS[,1]/12
group.GS.1<-group.GS%>%filter(depth.rel<=0.2)%>%mutate(group.depth=1)
group.GS.2<-group.GS%>%filter(depth.rel>0.2&depth.rel<=0.4)%>%mutate(group.depth=2)
group.GS.3<-group.GS%>%filter(depth.rel>0.4&depth.rel<=0.6)%>%mutate(group.depth=3)
group.GS.4<-group.GS%>%filter(depth.rel>0.6&depth.rel<=0.8)%>%mutate(group.depth=4)
group.GS.5<-group.GS%>%filter(depth.rel>0.8&depth.rel<=1.0)%>%mutate(group.depth=5)

group.GS<-rbind(group.GS.1,group.GS.2,group.GS.3,group.GS.4,group.GS.5)
colnames(group.GS)<-c("depth","depth.rel","group.depth")

GS.HESC.ANC.most$ANC.ratio.df$depth.group<-group.GS$group.depth[match(GS.HESC.ANC.most$ANC.ratio.df$treeDist,group.GS$depth)]

GS.HESC.ANC.most.new <-GS.HESC.ANC.most$ANC.ratio.df%>%na.omit()
GS.HESC.ANC.most.new<-GS.HESC.ANC.most.new %>% group_by(depth.group) %>% 
  summarise(mean_share_ratio=mean(share_ratio),
            mean_boot.share_ratio.mean= mean(boot.share_ratio.mean),
            mean_boot.share_ratio.sd=mean(boot.share_ratio.sd)) %>%
  as.data.frame()

plot.ShareANC.errorbar(GS.HESC.ANC.most.new, "GS-HESC")
=======
# FigureS_2C

   ###################### add some gene of lung progenitor or day3/day6 marker ####################

day0.gene<-c("SIX3", "OTX2", "PAX6","NODAL","SOX2","ETV5","PDPN","LPCAT1","SCGB3A2","NANOG","OCT4")

day3.gene<-c("OTX2","PAX6","NODAL","GATA4","GATA6","SOX17","LAMP3","ETV5")

day6.gene<-c("OTX2","SALL1","ZFP42","FOXC1","HHEX","DMBX1","GATA3","OTX1","ALX1",
             "MAPK10","GATA4","SOX15","NR6A1","TAF4B","ISL1","EOMES","PITX2") # day6 vs day15+

day15.positive.gene.vs.neural<-c("GATA6","FOXP1","FOXP2","GRHL2","AHR","HOXB2","ELF3","HNF1B","KLF5","HOXA1","FOSL2","FOXA1") #(vs neural nkx2.1)

day15.positive.gene.vs.neg<-c("SFTA3","SPOCK3","PALMD","FAM189A2","CRH","NKX2-1","NKX2-1-AS1","GRM8","BMP3","WDR49")#(day15 postive vs day15 nkx2.1 negative)

neural.nkx2.1.gene.vs.day15.pos<-c("SIX3","OTX1","OTX2","PAX6","SIX3","DBX1","OTX2","OTX1","FOXD1","PAX6","HESX1","LHX5","ZEB2","ZEB1") # neural nkx2-1+ vs day15 nkx2.1+

HESC.CBRAD5_Seuratobj <- AddModuleScore(HESC.CBRAD5_Seuratobj, features =list(day15.positive.gene.vs.neural), name ="day15.positive.gene.vs.neural")
day15.pos.vs.neural<-FeaturePlot_scCustom(seurat_object = HESC.CBRAD5_Seuratobj, features = "day15.positive.gene.vs.neural1",label = T,pt.size = 1)
#day0,day3,day6,day15.pos.vs.neural,neural.vs.day15.pos,day15.pos.vs.neg
