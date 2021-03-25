library(synapser)
library(tidyverse)
synLogin()
###DGIDB contains drug-target relationships beyond binding/inhibiting, so need to eliminate those for the purposes of this db
###Then, annotate remaining molecules through pubchem id exchange etc

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/polypharmacology-db/master/curateDGIDBstructures_db_v2.R"

types <- c("inhibitor",  "antagonist", "agonist", "binder",  "modulator", "blocker", "channel blocker", 
           "positive allosteric modulator", "allosteric modulator", "activator", "inverse agonist", "partial agonist",
           "activator,channel blocker", "gating inhibitor", "agonist,antagonist", "agonist,allosteric modulator", "activator,antagonist",
           "stimulator", "negative modulator", "allosteric modulator,antagonist",  "channel blocker,gating inhibitor antagonist,inhibitor",            
           "inhibitory allosteric modulator")  

dgidb.targets <- read.table(synGet("syn12684108")$path, sep = "\t", quote = "",header = T, strip.white = TRUE) %>%
  filter(interaction_claim_source != "ChemblInteractions") %>% 
  filter(drug_claim_primary_name != "", gene_name != "") %>% 
  filter(interaction_types %in% types)

chembl_struct <- read.table(synGet("syn12973248")$path, sep = "\t", header = T, comment.char = "") %>% 
  select(1,4) %>% 
  distinct()

chembl_names <- read.table(synGet("syn12972665")$path, sep = "\t", header = T, quote = "", comment.char = "") %>% 
  select(1,4) %>% 
  distinct() %>% 
  inner_join(chembl_struct) %>% 
  select(-molregno)

colnames(chembl_names) <- c("drug_chembl_id", "original_smiles")

dgidb.targets.chemblids <- filter(dgidb.targets, drug_chembl_id != "") %>% 
  left_join(chembl_names) 

dgidb.targets.noid <- filter(dgidb.targets, drug_chembl_id == "")
  
unmapped <- dgidb.targets.chemblids %>% 
  filter(is.na(original_smiles)) %>% 
  select(-original_smiles) %>% 
  bind_rows(dgidb.targets.noid)

unmapped$drug_name <- as.character(unmapped$drug_name)
unmapped$drug_claim_primary_name<-as.character.Date(unmapped$drug_claim_primary_name)
unmapped$drug_name[unmapped$drug_name==""] <- unmapped$drug_claim_primary_name[unmapped$drug_name==""]

old_dgidb_struct <- read.table(synGet("syn11832827", version = 5)$path, sep = "\t", header = T, comment.char = "") %>% 
  distinct() %>% 
  set_names(c("drug_name", "original_smiles")) %>% 
  right_join(unmapped) %>% 
  filter(!is.na(original_smiles))

dgidb_structures <- dgidb.targets.chemblids %>% 
  filter(!is.na(original_smiles)) %>% 
  bind_rows(old_dgidb_struct) %>% 
  select(drug_name, original_smiles)

still_unmapped <- read.table(synGet("syn11832827", version = 5)$path, sep = "\t", header = T, comment.char = "") %>% 
  distinct() %>% 
  set_names(c("drug_name", "original_smiles")) %>% 
  right_join(unmapped) %>% 
  filter(is.na(original_smiles)) %>% 
  select(-original_smiles)

names <- as.data.frame(still_unmapped$drug_name) %>% distinct()
write.table(names, "NoGit/dgidb_names_v2.txt", sep = "\t", quote = F, row.names = F)


# ###export, run through pubchem id exhange, manually curate remaining with chemspider if possible
struct_cur <- read.table("Data/dgidb_3_0_2_idex_structures.txt", sep = "\t", header = F, comment.char = "", quote = "\"") %>% 
  set_names("drug_name", "original_smiles") %>%
  filter(original_smiles != "") %>%
  group_by(drug_name) %>%
  top_n(1) %>%
  ungroup() %>%
  distinct() %>% 
  bind_rows(dgidb_structures)

write.table(struct_cur, "Data/dgidb_structures_curated_3_0_2_v1.txt", sep = "\t", row.names = F, quote = T)
synStore(File("Data/dgidb_structures_curated_3_0_2_v1.txt", parentId = "syn12683808"), executed = this.file, 
         used = c("syn11832827","syn12972665","syn12684108","syn12684108"))
