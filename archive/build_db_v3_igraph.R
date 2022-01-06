
### create igraph object from db
##for fendr

this.file <- "https://raw.githubusercontent.com/Sage-Bionetworks/polypharmacology-db/develop/Database/build_db_v3_igraph.R"

library(tidyverse)
library(synapser)
synLogin()
library(igraph)

db <- readRDS(synGet("syn17091507")$path) %>%
  filter(!is.na(hugo_gene)) %>%
  select(internal_id, hugo_gene, mean_pchembl, n_quantitative, n_qualitative) %>%
  group_by(internal_id, hugo_gene) %>%
  mutate(total_n = sum(n_quantitative, n_qualitative, na.rm = T)) %>%
  ungroup()

db.names <- readRDS(synGet("syn17091507")$path) %>%
  filter(!is.na(hugo_gene)) %>%
  select(internal_id, std_name) %>%
  distinct()

db.igraph <- graph.data.frame(db)

saveRDS(db.igraph, "drug-target_explorer_igraph.rds")
saveRDS(db.names, "drug-target_explorer_igraph_name_map.rds")

synStore(File("drug-target_explorer_igraph.rds", parentId = "syn11802193"),
         used = "syn11712148",
         executed = this.file)

synStore(File("drug-target_explorer_igraph_name_map.rds", parentId = "syn11802193"),
         used = "syn11712148",
         executed = this.file)


