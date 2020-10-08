library(targets)
# This is an example target script.
# Read the tar_script() help file for details.

# Define custom functions and other global objects.
# This is where you write source(\"R/functions.R\")
# if you keep your functions in external scripts.
source("R/functions.R")

# Set target-specific options such as packages.
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("dplyr", "rphylopic", "tidyr"))

# Define targets
targets <- list(
  tar_target(phylopic_count, image_count(),
             cue = tar_cue("always")),
  
  tar_target(pic_data, 
             rphylopic::image_list(length = phylopic_count, 
                                   options=c('credit','licenseURL','pngFiles','submitted','submitter',
                                             'svgFile','taxa','canonicalName','string','uri','type',
                                             'citationStart', 'directNames')))
)

# End with a call to tar_pipeline() to wrangle the targets together.
# This target script must return a pipeline object.
tar_pipeline(targets)
