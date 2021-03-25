library(synapser)
library(tidyverse)

synLogin()

input_df <- read.csv('remaining_unmatched_nf_drugs_11_26_2018.csv') %>% distinct()

mint_new_ids <- function(input_df){
  
  new_minimum_id <- synTableQuery('SELECT max(internal_id) FROM syn17015110')$asDataFrame()[1,1] + 1 
  
  df <- input_df %>% 
    distinct() %>% 
    group_by(inchikey) %>% 
    add_column('internal_id' = group_indices(.)+new_minimum_id) %>% 
    gather(foo, names, c("std_name", "synonym")) %>% 
    select(-foo) %>% 
    ungroup()
  
  syn17015110 <- df %>% 
    group_by(inchikey) %>% 
    slice(1) %>% 
    mutate(std_name = names) %>% 
    select(internal_id, std_name, inchikey, inchi, std_smiles) %>% 
    ungroup() %>% 
    distinct() %>% 
    mutate_all(as.character)
  
  syn17015111 <- df %>% 
    mutate(common_name = names) %>% 
    select(internal_id, common_name) %>% 
    distinct() %>% 
    mutate_all(as.character)
  
  syn17015113 <- df %>% 
    mutate(original_smiles = smiles) %>% 
    mutate(external_id = NA) %>% 
    mutate(database = "nf_preclinical") %>% 
    mutate(v2_smiles = NA) %>% 
    select(internal_id, v2_smiles, original_smiles, external_id, database) %>% 
    distinct() %>% 
    mutate_all(as.character)
  
  
  syn17015110_new <- synTableQuery("select * from syn17015110", includeRowIdAndRowVersion = F)$asDataFrame() %>%
    mutate_all(as.character) %>% 
    bind_rows(syn17015110)
  
  syn17015111_new <- synTableQuery("select * from syn17015111", includeRowIdAndRowVersion = F)$asDataFrame() %>% 
    mutate_all(as.character) %>% 
    bind_rows(syn17015111)
  
  syn17015113_new <- synTableQuery("select * from syn17015113", includeRowIdAndRowVersion = F)$asDataFrame() %>% 
    mutate_all(as.character) %>% 
    bind_rows(syn17015113) 
  
  synBuildTable("v3.0.1 IDs, Structures, Standard Name", "syn11672851", syn17015110_new) %>% synStore()
  synBuildTable("v3.0.1 IDs, Synonyms", "syn11672851", syn17015111_new) %>% synStore()
  
  list(
    Column(name = "internal_id", columnType = "INTEGER"),
    Column(name = "v2_smiles", columnType = "LARGETEXT"),
    Column(name = "internal_id_v2", columnType = "STRING", maximumSize = 50),
    Column(name = "original_smiles", columnType = "LARGETEXT"),
    Column(name = "external_id", columnType = "LARGETEXT"),
    Column(name = "database", columnType = "STRING", maximumSize = 50)) %>% 
  Schema(name = "v3.0.1 IDs, v2 IDs, and Source Structures", columns = ., parent = "syn11672851") %>% 
  Table(., syn17015113_new) %>%
  synStore()
  
}

mint_new_ids(input_df)

##also add all synonyms for previously mapped in drugs (separate script)
syn17090820_addition <- synTableQuery("select * from syn17090820", includeRowIdAndRowVersion = F)$asDataFrame() %>% 
  mutate(name_id = paste0(common_name, "_", internal_id))
  
input_df_2 <- read.csv('nf_dtex_mapped_drugs.csv') %>% 
  select(common_name, alt_common_name, internal_id) %>% 
  gather(foo, common_name, c("common_name", "alt_common_name")) %>% 
  select(-foo) %>% 
  ungroup() %>% 
  distinct() %>% 
  mutate(name_id = paste0(common_name, "_", internal_id)) %>% 
  filter(!name_id %in% syn17090820_addition$name_id) %>% 
  select(-name_id)

synStore(Table("syn17090820", input_df_2))
