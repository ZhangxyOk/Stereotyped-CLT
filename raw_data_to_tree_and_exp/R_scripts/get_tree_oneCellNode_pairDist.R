suppressMessages({
  library(dplyr)
  library(parallel)
  library(ape)
  library(ggtree)
  library(argparse)
})

## accept argument from command line
parser <- ArgumentParser(description = "get oneCelll node pair distance in tree")
parser$add_argument("-i", "--intree", type="character", help="PATH of tree, nwk")
parser$add_argument("-o", "--outDist", type="character", help="PATH of output tree dist")

args <- parser$parse_args()

## run
tree <- read.tree(args$intree)

# get oneCell node pair
oneCellNode <- tree$tip.label[grepl("_1$", tree$tip.label)]
oneCellNode.pair <- combn(oneCellNode, 2)
oneCellNode.pair <- data.frame(node1=oneCellNode.pair[1,], node2=oneCellNode.pair[2,], stringsAsFactors = F)

# get tree tibble
inSubtrees <- subtrees(tree)
distMatrix <- dist.nodes(tree)
treeTibble <- as_tibble(tree)

getDist <- function(nodeName, node){
  n2 <- dplyr::filter(treeTibble, label == nodeName) %>% `[`("node") %>% unlist()
  dis <- distMatrix[node, n2]
  return(dis)
}

all.subtree.depth <- lapply(seq(length(inSubtrees)), function(s){
  subtree <- inSubtrees[[s]]
  s.name <- subtree$name
  treeDeep <- Map(getDist, subtree$tip.label, s.name) %>% unlist() %>% as.vector() %>% max()
  return(data.frame(subtree = s.name,
                    subtree.depth = treeDeep,
                    stringsAsFactors = F))
}) %>% bind_rows()

oneCellNode.pair.dist <- mclapply(seq(nrow(oneCellNode.pair)), function(i){
  node1 <- oneCellNode.pair$node1[i]
  node2 <- oneCellNode.pair$node2[i]
  n1 <- filter(treeTibble, label == node1) %>% `[`("node") %>% unlist()
  n2 <- filter(treeTibble, label == node2) %>% `[`("node") %>% unlist()
  # dist1
  n1n2.mrca <- getMRCA(tree, c(node1, node2))
  treedist <- max(distMatrix[n1, n1n2.mrca], distMatrix[n2, n1n2.mrca])
  # dist normalize as subtree depth
  treedist2 <- filter(all.subtree.depth, subtree == n1n2.mrca) %>% `[`("subtree.depth") %>% unlist()
  return(data.frame(node1=node1,
                    node2=node2,
                    mrca = n1n2.mrca,
                    treeDist=treedist,
                    treeDist.asDepth = treedist2,
                    stringsAsFactors = F))
}, mc.cores = 40) %>% bind_rows(.)



# write out results
write.table(oneCellNode.pair.dist, args$outDist, sep = '\t', quote = F, row.names = F)




