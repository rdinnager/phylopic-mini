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
##' @param otol
##' @param pic_data
map_phylopic_to_otol <- function(otol, pic_data) {
  
  
  tax <- rphylopic::
  
  spec_names_list <- purrr::map(pic_data,
                             ~sapply(.x$directNames, function(x) x$string)) %>%
    purrr::compact()
  
  spec_names_o <- otol$tip.label %>%
    strsplit("_", fixed = TRUE) %>%
    sapply(function(x) paste(x[1:2], collapse = "_"))
  
  spec_names_p <- spec_names_list %>%
    unlist() %>%
    stringr::str_remove_all("[[:punct:]]") %>%
    strsplit(" ", fixed = TRUE) %>%
    sapply(function(x) paste(x[1:2], collapse = "_")) %>%
    unique()
  
  test <- taxizedb::classification(spec_names_p)
  
  test_tree <- taxize::class2tree(test)
  
  classified <- purrr::map_int(test,
                               ~sum(!is.na(.x)))
 
  gen_names_o <- otol$node.label
 sum(spec_names_p %in% spec_names_o)
  
}





