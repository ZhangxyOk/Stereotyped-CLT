library(rjson)
library(ggplot2)
library(ggraph)
library(igraph)
library(scales)
library(tidyverse)
library(RColorBrewer)
library(parallel)
library(dplyr)
library(ggrepel)

# input: json file from mdelta alignment between each sample pair
# within sample
df<-"g2_cbrad5_onecell_g2_cbrad5_onecell_top100_diff90_pv-0.2_mg100.0.json" # loop for CBRAD5 sample pair
path<-"~/fig5/"
g2.top100<-fromJSON(file=paste0(path,df)) # loop for CBRAD5 sample pair
top100<-g2.top100 #
match.tree<-mclapply(seq(length(top100)),mc.cores=50,function(i){
  group<-"g2_cbrad5" # loop for CBRAD5 sample pair
  df.match<-top100[[i]]
  node1<-df.match$Root1_label
  node1<-paste0(group,"_",node1) 
  node2<-df.match$Root2_label
  node2<-paste0(group,"_",node2)  
  size1<-df.match$Root1_leaves_count
  size2<-df.match$Root2_leaves_count
  data1<-data.frame(from=group,
                    to=node1,  
                    size=size1)
  data2<-data.frame(from=group,
                    to=node2,
                    size=size2)
  data3<-data.frame(from=node1,
                    to=node2,
                    size=NA) # size of "to" node
  data<-rbind(data1,data2,data3)
  
  return(data)
})%>%bind_rows%>%unique()

match.tree.g2<-match.tree # loop for CBRAD5 sample pair

edges.df<-mclapply(seq(length(top100)),mc.cores=50,function(i){
  group<-"g2_cbrad5" # loop for CBRAD5 sample pair
  df.match<-top100[[i]]
  node1<-df.match$Root1_label
  node1<-paste0(group,"_",node1) 
  node2<-df.match$Root2_label
  node2<-paste0(group,"_",node2)
  data1<-data.frame(from="origin",
                    to=group)
  data2<-data.frame(from=group,
                    to=node1)
  data3<-data.frame(from=group,
                    to=node2)
  data<-rbind(data1,data2,data3)
  return(data)
})%>%bind_rows%>%unique()
edges.g2<-edges.df  #loop for CBRAD5 samples

connect<-match.tree%>%filter(.,from!="g2_cbrad5")%>%dplyr::select("from","to") # loop for CBRAD5 sample pair
value<-table(connect$from)%>%as.data.frame()
colnames(value)<-c("from","freq")
connect$value<-value$freq[match(connect$from,value$from)]
connect.g2<-connect # loop for CBRAD5 sample pair

# between sample

df<-"g2_cbrad5_onecell_g11_cbrad5_onecell_top100_diff90_pv-0.2_mg100.0.json" # loop for CBRAD5 sample pair
path<-"~/fig5/" 
g11.g2.top100<-fromJSON(file=paste0(path,df)) # loop for CBRAD5 sample pair
top100<-g11.g2.top100 # loop for CBRAD5 sample pair
match.tree<-mclapply(seq(length(top100)),mc.cores=50,function(i){
  group.1<-"g11_cbrad5" # loop for CBRAD5 sample pair
  group.2<-"g2_cbrad5" # loop for CBRAD5 sample pair
  df.match<-top100[[i]]
  node1<-df.match$Root1_label
  node1<-paste0(group.1,"_",node1) 
  node2<-df.match$Root2_label
  node2<-paste0(group.2,"_",node2) 
  size1<-df.match$Root1_leaves_count
  size2<-df.match$Root2_leaves_count
  data1<-data.frame(from=group.1,
                    to=node1,  
                    size=size1) # size of verticles (group to)
  data2<-data.frame(from=group.2,
                    to=node2,
                    size=size2)  # size of verticles (group to)
  data3<-data.frame(from=node1,
                    to=node2,
                    size=NA) 
  data=rbind(data1,data2,data3)
  return(data)
})%>%bind_rows%>%unique()

index<-c(which(grepl("root",match.tree$from)),which(grepl("root",match.tree$to)))%>%unique()
match.tree.g11.g2<-match.tree  # loop for CBRAD5 sample pair

edges.df<-mclapply(seq(length(top100)),mc.cores=50,function(i){
  group.1<-"g11_cbrad5" # loop for CBRAD5 sample pair
  group.2<-"g2_cbrad5" # loop for CBRAD5 sample pair
  df.match<-top100[[i]]
  node1<-df.match$Root1_label
  node1<-paste0(group.1,"_",node1) 
  node2<-df.match$Root2_label
  node2<-paste0(group.2,"_",node2)
  data1<-data.frame(from="origin",
                    to=group.1)
  data2<-data.frame(from="origin",
                    to=group.2)
  data3<-data.frame(from=group.1,
                    to=node1)
  data4<-data.frame(from=group.2,
                    to=node2)
  data<-rbind(data1,data2,data3,data4)
  return(data)
})%>%bind_rows%>%unique()
edges.g11.g2<-edges.df  # loop for CBRAD5 sample pair

connect<-match.tree%>%filter(.,from!="g11_cbrad5"&from!="g2_cbrad5")%>%dplyr::select("from","to") # loop for CBRAD5 sample pair
value<-table(connect$from)%>%as.data.frame()
colnames(value)<-c("from","freq")
connect$value<-value$freq[match(connect$from,value$from)]
connect.g11.g2<-connect # loop for CBRAD5 sample pair

## connection plot
match.df.all<-rbind(match.tree.a1,match.tree.g2,match.tree.g11,match.tree.a1.g2,match.tree.a1.g11,match.tree.g11.g2)
connect.all<-rbind(connect.a1,connect.g2,connect.g11,connect.a1.g2,connect.a1.g11,connect.g11.g2)
edges.all<-rbind(edges.a1,edges.g2,edges.g11,edges.a1.g2,edges.a1.g11,edges.g11.g2)

df<-connect.all%>%dplyr::select("from","value")%>%unique
df<-aggregate(df$value, by=list(from=df$from), FUN=sum)
colnames(df)<-c("from","value")
connect.all$value<-df$value[match(connect.all$from,df$from)]

vertices<- data.frame(name = unique(c(as.character(connect.all$from), as.character(connect.all$to))))
vertices.size<-vertices
vertices.size$size<-match.df.all$size[match(vertices$name,match.df.all$to)]
vertices.value<-vertices
vertices.value$value<-connect.all$value[match(vertices$name,connect.all$from)]
vertices.value$value[which(is.na(vertices.value$value))]<-1

vertices<-vertices.value
edges<-edges.all%>%unique
a<-data.frame(name="a1_cbrad5",value=sum(vertices$value[which(grepl("a1_cbrad5",vertices$name))],na.rm=TRUE))
b<-data.frame(name="g2_cbrad5",value=sum(vertices$value[which(grepl("g2_cbrad5",vertices$name))],na.rm=TRUE))
c<-data.frame(name="g11_cbrad5",value=sum(vertices$value[which(grepl("g11_cbrad5",vertices$name))],na.rm=TRUE))
d<-data.frame(name="origin",value=sum(a$value,b$value,c$value))
vertices<-rbind(vertices,a,b,c,d)

vertices$group  <-edges$from[ match( vertices$name, edges$to ) ]
vertices$id <- NA
myleaves <- which(is.na( match(vertices$name, edges$from) ))
nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)
vertices$angle <- 90 - 360 * vertices$id / nleaves

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
vertices$hjust <- ifelse( vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
vertices$angle <- ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

# Basic usual argument
# Create a graph object
mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices ) 
# The connection object must refer to the ids of the leaves:
connect<-connect.all
from  <-  match( connect$from, vertices$name)
to  <-  match( connect$to, vertices$name)
pdf("~/fig5/top100_network_Quantitive_alignment_with_and_between_samples_prune_perc20.text.pdf",height=10,width=10) 
# Basic usual argument
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=value, alpha=0.2))+
  geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, colour="skyblue", width=0.9) +
  geom_node_text(aes(x = x*1.1, y=y*1.1, filter = leaf, label=name, angle = angle, hjust=hjust), size=1.5, alpha=1) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.3, 1.3), y = c(-1.3, 1.3))
dev.off()
