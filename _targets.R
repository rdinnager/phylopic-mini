library(targets)
# This is an example target script.
# Read the tar_script() help file for details.

# Define custom functions and other global objects.
# This is where you write source(\"R/functions.R\")
# if you keep your functions in external scripts.
source("R/functions.R")

# Set target-specific options such as packages.
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("dplyr", "rphylopic", "tidyr", "purrr"))

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
  
  tar_target(png_files, paste0("http://phylopic.org", png_urls), format = "url", 
             pattern = map(png_urls)),
  
  tar_target(local_files, file.path("data/phylopic/pngs", basename(png_urls)), format = "file",
             pattern = map(png_files)),
  
  tar_target(pngs, download_phylopics(png_files, local_files), format = "file",
             pattern = map(png_files, local_files)),
  
  NULL
  
)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
tar_pipeline(targets)
