# load required libraries ======================================================
library(tidyverse)

# ==============================================================================
# Peform only if data is not loaded
# ==============================================================================
# source loading function ------------------------------------------------------
source("R\\loadSBTab.R")

# load tables ------------------------------------------------------------------
read_sbtab("SBTab")

# ==============================================================================
# Step 1: Retrieve genes related to pyruvate snythase reaction
# ==============================================================================
# get genes involed in pyruvate synthase reaction ------------------------------
gene_list <- reactions_table$`!GeneAssociation`[20] %>% 
  map(.f=~str_extract_all(.x, "WBGene\\d+")) %>% 
  unlist()

# filter gene table according to entries in gene_list --------------------------
genes_table %>% 
  filter(.$`!ID` %in% gene_list)

# single command to obtain a filtered gene table -------------------------------
reactions_table$`!GeneAssociation`[20] %>% 
  map(.f=~str_extract_all(.x, "WBGene\\d+")) %>% 
  unlist() %>% 
  tibble(`!ID` = .) %>% 
  left_join(., genes_table, by = c("!ID"))

# ==============================================================================
# Step 2: Retrieve chemical information for pyruvate
# ==============================================================================
# get chemical information for pyruvate ----------------------------------------
compounds_table %>% 
  filter(str_detect(.$`!ID`, "M_pyr"))

# get chemical information for pyruvate in the cytosol -------------------------
compounds_table %>% 
  filter(.$`!ID` == "M_pyr_c")

# ==============================================================================
# Step 3: Retrieve chemical information on multipe compounds
# ==============================================================================
# get all metabolites related to pyruvate synthase reaction --------------------
metabolite_list <- reactions_table$`!ReactionFormula`[20] %>% 
  map(.f=~str_extract_all(.x, "M_\\w+_(c|m|e|n)")) %>% 
  unlist()

# filter gene table according to entries in metabolite list --------------------
compounds_table %>% 
  filter(.$`!ID` %in% metabolite_list)

# single command to obtain a filtered compound table ---------------------------
reactions_table$`!ReactionFormula`[20] %>% 
  map(.f=~str_extract_all(.x, "M_\\w+_(c|m|e|n)")) %>% 
  unlist() %>% 
  tibble(`!ID` = .) %>% 
  left_join(., compounds_table, by = c("!ID"))

# ==============================================================================
# Step 4: Retrieve chemical information from neutral metabolitse
# ==============================================================================
# single command to obtain a filtered compound table and select only neutral ---
reactions_table$`!ReactionFormula`[20] %>% 
  map(.f=~str_extract_all(.x, "M_\\w+_(c|m|e|n)")) %>% 
  unlist() %>% 
  tibble(`!ID` = .) %>% 
  left_join(., compounds_table, by = c("!ID")) %>% 
  select(contains("neutral"))
