# load required libraries ======================================================
library(tidyverse)

# ==============================================================================
# Step 1: loading of metabolomics data and mapping file
# ==============================================================================
# read files -------------------------------------------------------------------
metabolomics_data <- read_tsv("MetabolomicsData\\metabolomics_table.tsv")
model_mapping <- read_tsv("MetabolomicsData\\metabo-model-mapping.tsv")

# ==============================================================================
# Step 2: filter metabolites on pathways of interest
# ==============================================================================
# filter metabolites from glycolysis, TCA and PPP ------------------------------
metabolomics_data_filtered <- metabolomics_data %>% 
  filter(str_detect(`METABOLIC PATHWAY`,
                    "Glycolysis|TCA|Pentose phosphate pathway"))

# ==============================================================================
# Step 3: reformat data table
# ==============================================================================
# generate long format from data table -----------------------------------------
metabolomics_long <- gather(metabolomics_data_filtered, key, value,
                            -`Compound (Metabolite)`, -`METABOLIC PATHWAY`,
                            convert = TRUE, na.rm = TRUE)

# generate long format and split sample column to sub groups -------------------
metabolomics_long <- gather(metabolomics_data_filtered, key, value,
                            -`Compound (Metabolite)`, -`METABOLIC PATHWAY`,
                            convert = TRUE, na.rm = TRUE) %>% 
  separate(into = c("time", "strain", "Day", "sampleID"),
           col = key,
           sep = "-")

# ==============================================================================
# Step 4: calculate group mean
# ==============================================================================
# group data and summarize -----------------------------------------------------
metabolomics_summary <- metabolomics_long %>% 
  group_by(`Compound (Metabolite)`, time, strain) %>% 
  summarize(mean = mean(value),
            sd = sd(value),
            n = n())

# ==============================================================================
# Step 5: add model metabolite names
# ==============================================================================
# use join to combine table ----------------------------------------------------
metabolomics_summary <- metabolomics_summary %>% left_join(.,
                                   model_mapping,
                                   by = c("Compound (Metabolite)" = "METABOLITE")) %>%
  filter(!is.na(MODEL))

# ==============================================================================
# Step 6: filter groups and time points of interest
# ==============================================================================
# use filter function to get strain and time points of interest ----------------
metabolomics_summary %>% filter(strain == "glp",
                                time %in% c("65", "73")) %>% 
  ungroup() %>% 
  unite(sampleID, strain, time) %>% 
  select(sampleID, MODEL, mean) %>% 
  spread(key = sampleID, value =mean) %>% 
  rename(ID = MODEL) %>% 
  write_csv("Escher\\escher_metabolites.csv")
