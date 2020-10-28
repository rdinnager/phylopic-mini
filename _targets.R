library(targets)
# This is an example target script.
# Read the tar_script() help file for details.

# Define custom functions and other global objects.
# This is where you write source(\"R/functions.R\")
# if you keep your functions in external scripts.
source("R/functions.R")

# Set target-specific options such as packages.
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("dplyr", "rphylopic", "tidyr", "purrr", "ape",
                            "ape", "igraph"))

# Define targets
targets <- list(
  tar_target(phylopic_count, image_count(),
             cue = tar_cue("always")),
  
  tar_target(pic_data, 
             rphylopic::image_list(length = phylopic_count, 
                                   options=c('credit','licenseURL','pngFiles','submitted','submitter',
                                             'svgFile','taxa','canonicalName','string','uri','type',
                                             'citationStart', 'directNames'))),
  
  tar_target(max_res, map_int(pic_data, ~length(.x$pngFiles))),
  
  tar_target(png_urls, map2_chr(pic_data, max_res,
                                ~.x$pngFiles[[.y]]$url)),
  
  tar_target(png_files, paste0("http://phylopic.org", png_urls), 
             pattern = map(png_urls)),
  
  tar_target(local_files, file.path("data/phylopic/pngs", basename(png_urls)),
             pattern = map(png_urls)),
  
  tar_target(pngs, download_phylopics(png_files, local_files),
             pattern = map(png_files, local_files)),
  
  tar_target(files_unedited, file.path("data/phylopic/mini_unedited", basename(local_files)),
             pattern = map(local_files)),
  
  tar_target(mini_unedited, make_mini_unedited(local_files, files_unedited),
             pattern = map(local_files, files_unedited)),
  
  tar_target(phylopic_taxo, get_taxonomy(pic_data),
             pattern = map(pic_data),
             iteration = "list",
             error = "continue"),
  
  tar_target(phylopic_tree, make_phylopic_tree(phylopic_taxo, pic_data)),
  
  NULL
  
)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
tar_pipeline(targets)
