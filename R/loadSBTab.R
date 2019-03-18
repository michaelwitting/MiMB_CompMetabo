# load required libraries ======================================================
library(tidyverse)

# helper function ==============================================================
read_sbtab <- function(folderPath) {
  
  # get all files
  sbtab_files <- list.files(folderPath,
                            pattern = "-SBTab.tsv$",
                            full.names = TRUE)
  
  # make new list
  sbtab_list <- list()
  
  # iterate over all files
  for(i in 1:length(sbtab_files)) {
    
    #get current file and add to list
    sbtab_file <- sbtab_files[i]
    sbtab_list[[i]] <- read_tsv(sbtab_file, comment = "!!")
    
  }
  
  # correct names
  sbtab_names <- str_replace_all(basename(sbtab_files), "-SBTab.tsv", "")
  names(sbtab_list) <- paste0(sbtab_names, "_table")
  
  # make tibble for eacht table
  for(i in 1:length(sbtab_list)) {
    
    assign(names(sbtab_list)[i], sbtab_list[[i]], envir = parent.frame())
    
  }
  
}
