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
# Step 1: Retrieve all metabolites and their neutral information
# ==============================================================================
# metabolites not of interest
metabolite_exclusion <- c("M_h_c", "M_h2o_c", "M_h_m", "M_h2o_c", "M_q_m", "M_qh2_m")

# isolate all metabolites from all reactions and leave unique ------------------
neutral_metabolites <- reactions_table$`!ReactionFormula` %>% 
  map(.f=~str_extract_all(.x, "M_\\w+_(c|m|e|n)")) %>% 
  unlist() %>% 
  tibble(`!ID` = .) %>% 
  left_join(., compounds_table, by = c("!ID")) %>% 
  filter(!.$`!ID` %in% metabolite_exclusion) %>% 
  select(contains("neutral")) %>% 
  distinct()

# create new tibble that can be used with masstrixR ----------------------------
suspect_List <- tibble(
  id = neutral_metabolites$`!Notes:ChEBI_neutral`,
  smiles = NA,
  inchi = neutral_metabolites$`!Notes:InChI_neutral`,
  inchikey = neutral_metabolites$`!Notes:InChIKey_neutral`,
  formula = neutral_metabolites$`!Notes:FORMULA_Neutral`,
  name = neutral_metabolites$`!Notes:ChEBI_Name_neutral`,
  exactmass = NA 
)

# ==============================================================================
# Step 2: generate list with m/z values for metabolites
# ==============================================================================
# load masstrixR library -------------------------------------------------------
library(masstrixR)

# create compound list with adduct masses --------------------------------------
neg_adducts <- c("[M-H]-")
neg_suspect_List <- as_tibble(prepareCompoundList(suspect_List, adductList = neg_adducts))
