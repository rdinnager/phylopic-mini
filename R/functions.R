##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param png_files
##' @param local_files
download_phylopics <- function(png_files, local_files) {
  
  Sys.sleep(2)
  
  download.file(png_files, local_files, quiet = TRUE, mode = "wb")
  
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param pngs
make_mini_unedited <- function(pngs, files_unedited) {
  
  pic <- imager::load.image(pngs) %>%
    imager::flatten.alpha() %>%
    imager::grayscale()
  
  h <- imager::height(pic)
  w <- imager::width(pic)
  
  dim_diff <- h - w
  
  if(dim_diff >= 0) {
    xy = "x"
  } else {
    xy = "y"
    dim_diff <- -dim_diff
  }

  pic <- pic %>%
    imager::pad(dim_diff, xy, val = 1) %>%
    imager::resize(size_x = 28L,
                   size_y = 28L,
                   interpolation_type = 6) 
  
  pic <- 1 - pic
  
  imager::save.image(pic, files_unedited, quality = 1)
  
  pic
  
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param pic_data
get_taxonomy <- function(pic_data) {
  
  Sys.sleep(2)
  
  purrr::map(pic_data[[1]]$directNames,
             ~rphylopic::name_taxonomy(.x$uid,
                                       supertaxa = "all",
                                       options = c("string",
                                                   "citationStart")))
  
  
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param phylopic_taxo
##' @param pic_data
make_phylopic_tree <- function(phylopic_taxo, pic_data) {
  
  ## get longest lineage for each tip
  
  lens <- purrr::map(phylopic_taxo,
                     ~sapply(.x, function(x) nrow(x$taxa)))
  
 
  lineages <- purrr::map2(phylopic_taxo,
                          lens,
                         ~if(length(.y) > 0) {.x[[which.max(.y)]]}
                           else {NULL})
  
  tip_taxa <- purrr::map2(pic_data,
                           lens,
                           ~if(length(.y) > 0) {.x$directNames[[which.max(.y)]]}
                           else {NULL}) 
  
  names_uncited <- map(lineages,
                       ~list(string = stringr::str_sub(.x$taxa$canonicalName$string,
                                         1L,
                                         ifelse(is.na(.x$taxa$canonicalName$citationStart),
                                                -1L,
                                                .x$taxa$canonicalName$citationStart - 1L)),
                             inclusions = .x$inclusions))
  
  tip_names <- map(tip_taxa,
                   ~stringr::str_sub(.x$string,
                                     1L,
                                     ifelse(is.null(.x$citationStart),
                                            -1L,
                                            .x$citationStart - 1L)))
  
  ## add edge for missing tips
  missing_tips <- map2(names_uncited,
                       tip_names,
                       possibly(~!.y %in% .x$string,
                                otherwise = FALSE)) %>%
    modify_if(~length(.x) == 0,
              ~FALSE) %>%
    unlist()
  
  #x <- names_uncited[missing_tips][[1]]
  add_tip_edge <- function(x, tip_name) {
    
    if(length(x$inclusions) > 0) {
      graph_df <- x$inclusions %>%
        as_tibble()
      
      graph_ig <- graph_from_data_frame(graph_df, directed = TRUE)
      
      ## find end vertice
      v_degs <- degree(graph_ig, mode = "out")
      end <- names(v_degs[v_degs == 0])[1] %>%
        as.numeric()
      
      x$string <- c(x$string, tip_name)
      x$inclusions <- rbind(x$inclusions,
                            c(end, max(x$inclusions) + 1L))
    
    }
    
    x
    
  }
  
  new_names <- map_if(seq_along(names_uncited),
                          ~missing_tips[.x],
                          ~add_tip_edge(names_uncited[[.x]],
                                        tip_names[[.x]]),
                      .else = ~names_uncited[[.x]]) 
  
  
  lineage_dfs <- imap_dfr(new_names,
                          possibly(~.x$string[.x$inclusions + 1L] %>%
                                     matrix(nrow = nrow(.x$inclusions)) %>%
                                     as_tibble() %>%
                                     mutate(tip_num = .y),
                                   otherwise = tibble())) %>%
    drop_na()
     
  
  tree_ig <- graph_from_data_frame(lineage_dfs)
  
  tree_ig <- igraph::simplify(tree_ig, edge.attr.comb = toString)

  reachable <- subcomponent(tree_ig, "Pan-Biota",
                            "out")
  
  tree_ig <- induced_subgraph(tree_ig, reachable)
  
  tree <- dominator_tree(tree_ig, "Pan-Biota")$domtree
  
  #sum(tip_names %in% names(V(tree)))
  
  
  
  return(tree)
  
  out <- degree(tree, mode = "out")
  tips <- names(out)[out == 0]
  n_tips <- length(tips)
  
  reorder_num <- numeric(length(V(tree)))
  
  reorder_num[out == 0] <- 1:sum(out == 0)
  reorder_num[reorder_num == 0] <- (sum(out == 0) + 1L):(sum(out == 0) + sum(reorder_num == 0))
  
  tree <- permute(tree, reorder_num)
  
  out2 <- degree(tree, mode = "out")
  
  phylo_edge <- as_edgelist(tree, names = FALSE)
  
  phylo <- list(edge = phylo_edge,
                tip.label = names(V(tree)[1:n_tips]),
                node.label = names(V(tree)[(n_tips + 1L):length(V(tree))]),
                Nnode = length(V(tree)) - n_tips,
                edge.length = rep(1, nrow(phylo_edge)))
  
  class(phylo) <- "phylo"
  phylo <- ladderize(phylo)
  
  tt <- tempfile(fileext = ".nw")
  write.tree(phylo, tt)
  phylo <- read.tree(tt)
  
  sum(gsub(" ", "_", unlist(tip_names)) %in% phylo$tip.label)
  sum(gsub(" ", "_", unlist(tip_names)) %in% phylo$node.label)
  
  sum(gsub(" ", "_", unlist(tip_names)) %in% phylo$tip.label) +
  sum(gsub(" ", "_", unlist(tip_names)) %in% phylo$node.label)
  
  phylo <- collapse.singles(phylo)
  
  not_in_tree <- unlist(tip_names[!tip_names %in% names(V(tree))])
  
  where_are_they <- purrr::map(not_in_tree,
                               ~grep(.x, names(V(tree_ig)), value = TRUE))
  
  which(tip_names == "Erlikosaurus")
  
  
}





