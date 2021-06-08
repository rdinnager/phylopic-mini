## function to convert igraph tree to phylo object
as_phylo.igraph <- function(ig) {
  out <- degree(ig, mode = "out")
  tips <- names(out)[out == 0]
  n_tips <- length(tips)
  
  reorder_num <- numeric(length(V(ig)))
  
  reorder_num[out == 0] <- 1:sum(out == 0)
  reorder_num[reorder_num == 0] <- (sum(out == 0) + 1L):(sum(out == 0) + sum(reorder_num == 0))
  
  ig <- permute(ig, reorder_num)
  
  out2 <- degree(ig, mode = "out")
  
  phylo_edge <- as_edgelist(ig, names = FALSE)
  
  phylo <- list(edge = phylo_edge,
                tip.label = names(V(ig)[1:n_tips]),
                node.label = names(V(ig)[(n_tips + 1L):length(V(ig))]),
                Nnode = length(V(ig)) - n_tips)
  
  if(!is.null(igraph::vertex_attr(ig, "brlen"))) {
    phylo$edge.length <- igraph::vertex_attr(ig, "brlen")[phylo_edge[ , 2]]
  }
  
  class(phylo) <- "phylo"
  
  ## trick to fix formatting of phylo object, save to newick, then read it back in
  ## works most times
  tt <- tempfile(fileext = ".nw")
  write.tree(phylo, tt)
  phylo <- read.tree(tt)
  
  phylo <- ladderize(phylo)
  
  phylo
}

## function that takes a phylo object, and a vector of nodes (numbers or labels)
## and outputs a phylo with only the requested nodes, with nodes as tips
## if they are no longer internal
keep_tips_and_nodes <- function(phy, nodes, collapse_singles = FALSE) {
  
  temp_phy <- phy
  temp_phy$tip.label <- as.character(seq_along(temp_phy$tip.label))
  temp_phy$node.label <- as.character(length(temp_phy$tip.label) +
                                        seq_len(temp_phy$Nnode))
  
  
  ig <- as.igraph(temp_phy, directed = TRUE)
  
  if(inherits(nodes, "character")) {
    nodes <- as.character(c(which(phy$tip.label %in% nodes),
                            length(phy$tip.label) + which(phy$node.label %in% nodes)))
  } else {
    nodes <- as.character(nodes)
  }
  
  if(!is.null(temp_phy$edge.length)) {
    edges_subtending <- match(as.numeric(names(V(ig))), temp_phy$edge[ , 2])
    igraph::vertex_attr(ig, "brlen") <- temp_phy$edge.length[edges_subtending]
    igraph::vertex_attr(ig, "brlen")[is.na(igraph::vertex_attr(ig, "brlen"))] <- 0
  }
  
  ## find root
  degs <- igraph::degree(ig, mode = "in")
  root <- names(degs)[degs == 0]
  
  ## This is the bit that filters the nodes, it finds shortest paths from the
  ## root to each of the nodes (which for trees should be the one and only
  ## path between them) and keeps all nodes and edges it finds along the way
  paths_to_keep <- igraph::shortest_paths(ig, from = root, to = nodes, mode = "out")
  nodes_to_keep <- unique(unlist(paths_to_keep$vpath))
  ig <- igraph::induced_subgraph(ig, nodes_to_keep)
  
  new_tree <- as_phylo.igraph(ig)
  
  old_labels <- c(phy$tip.label, phy$node.label)
  new_tip_labels <- old_labels[as.numeric(new_tree$tip.label)]
  
  new_tree$tip.label <- new_tip_labels
  
  if(!is.null(phy$node.label)) {
    new_node_labels <- old_labels[as.numeric(new_tree$node.label)]
    new_tree$node.label <- new_node_labels
  }
  
  if(collapse_singles) {
    new_tree <- ape::collapse.singles(new_tree)
  }
  
  new_tree
  
}


set.seed(2)

library(ape)
library(tidygraph)
library(igraph)

## test on some random tree and node lists
phy <- ape::rcoal(100)
phy$node.label <- paste0("n", seq_len(phy$Nnode))

## sample 5 tips and 5 nodes randomly
nodes <- c(sample(phy$tip.label, 5), sample(phy$node.label, 5))

## show sampled nodes as red points on tree
plot(phy, type = "f", cex = 0.55)
nodelabels(node = c(which(phy$tip.label %in% nodes),
                    length(phy$tip.label) + which(phy$node.label %in% nodes)), 
           pch = 19,
           col = "red")

## test the function!
test <- keep_tips_and_nodes(phy, nodes)

## plot looks good:
plot(test, type = "f")
## rectangular plot shows that tips that were originally internal nodes
## have are not aligned to nodes that were originally tips, so branchlengths
## are handled correctly.

plot(test)

## tip labels match original nodes requested except for one, which stayed
## as an internal node in the filtered tree (n52)
nodes
test$tip.label

## has many extra nodes
test
## collapse.singles removes internal nodes along branches
collapse.singles(test)
