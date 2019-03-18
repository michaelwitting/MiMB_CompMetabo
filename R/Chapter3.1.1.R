# load required libraries ======================================================
library(tidyverse)

# ==============================================================================
# Step 1: loading data from SBTab files
# ==============================================================================
# source loading function ------------------------------------------------------
source("R\\loadSBTab.R")

# load tables ------------------------------------------------------------------
read_sbtab("SBTab")

# ==============================================================================
# Step 2: Working with the reaction table
# ==============================================================================
# print complete reaction table to console -------------------------------------
reactions_table

# print only reaction formulas to console --------------------------------------
reactions_table$`!ReactionFormula`

# print pyruvate snythase reaction ---------------------------------------------
reactions_table$`!ReactionFormula`[20]

# ==============================================================================
# Step 3: Retrieval of metabolites from reactions
# ==============================================================================
# isolate all compounds and count occurence ------------------------------------
metabolite_counts <- reactions_table$`!ReactionFormula` %>%
  map(.f=~str_extract_all(.x, "M_\\w+_(c|m|e|n)")) %>%
  unlist() %>% 
  tibble(`!ID` = .) %>% 
  add_count(`!ID`) %>% 
  distinct(`!ID`, .keep_all = TRUE) %>% 
  arrange(desc(n))

# filter to remove hub metabolites =============================================
hub_metabolites <- c("M_h_c", "M_h_m", "M_h2o_c", "M_h2o_m")

metabolite_counts %>% 
  filter(!`!ID` %in% hub_metabolites)

# ==============================================================================
# Step 4: Retrieval of gene associations from reactions
# ==============================================================================
# get genes involed in pyruvate synthase reaction ------------------------------
reactions_table$`!GeneAssociation`[20]

reactions_table$`!GeneAssociation`[20] %>% 
  map(.f=~str_extract_all(.x, "WBGene\\d+")) %>% 
  unlist()

# ==============================================================================
# Step 5: Retrieval of reactions containing a specific metabolite
# ==============================================================================
# get all reactions with a pyruvate involved -----------------------------------
filtered_reactions_table <- reactions_table %>%
  filter(str_detect(.$`!ReactionFormula`, "M_pyr"))

filtered_reactions_table
