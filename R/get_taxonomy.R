library(rphylopic)
library(purrr)
library(dplyr)
library(glue)
library(readr)
library(conflicted)

conflict_prefer("filter", "dplyr")

## collect full phylopic taxonomy, this might take awhile
root_uid <- "7FF2E439-19FE-4ADB-AC8B-0C8052AE313F"
# uids <- get_names(root_uid, subtaxa = "immediate", options = "string")$uid
get_descendents <- function(uuid, sup_name) {
  subtaxa <- get_names(uuid, subtaxa = "immediate", options = "string")
  subtaxa <- subtaxa %>%
    filter(uid != tolower(uuid)) 
  
  if(nrow(subtaxa) > 0) {
    subtaxa <- subtaxa %>%
      mutate(supertaxa = sup_name)
  }
  subtaxa
}


subtaxa <- map2_dfr(root_uid, "Panbiota",
                    ~get_descendents(.x, .y))

full_df <- subtaxa
i <- 1

saver <- file.path("data/taxonomy", glue("taxa_{i}.csv"))
if(!file.exists(saver)) {
  write_csv(subtaxa, saver)
} else {
  subtaxa <- read_csv(saver)
}

while(nrow(subtaxa) > 0) {
  i <- i + 1
  
  saver <- file.path("data/taxonomy", glue("taxa_{i}.csv"))
  
  if(!file.exists(saver)) {
  
    subtaxa <- map2_dfr(subtaxa$uid, subtaxa$name,
                     ~get_descendents(.x, .y))
    
  
    write_csv(subtaxa, saver)
    
  } else {
    
    subtaxa <- read_csv(saver)
    
  }
  
  full_df <- bind_rows(full_df,
                       subtaxa)
  
  cat("Added ", nrow(subtaxa), " taxa; Total: ", nrow(full_df), "; We have done ", i, " levels.\n")

}

write_csv(full_df, "data/phylopic_taxonomy.csv")
