library(synapser)
library(tidyverse)

###DGIDB contains drug-target relationships beyond binding/inhibiting, so need to eliminate those for the purposes of this db
###Then, annotate remaining molecules through pubchem id exchange etc

types <- c("inhibitor",  "antagonist", "agonist", "binder",  "modulator", "blocker", "channel blocker", 
           "positive allosteric modulator", "allosteric modulator", "activator", "inverse agonist", "partial agonist",
           "activator,channel blocker", "gating inhibitor", "agonist,antagonist", "agonist,allosteric modulator", "activator,antagonist",
           "stimulator", "negative modulator", "allosteric modulator,antagonist",  "channel blocker,gating inhibitor antagonist,inhibitor",            
           "inhibitory allosteric modulator")  

dgidb.targets <- read.table(synGet("syn11672978")$path, sep = "\t", quote = "",header = T, strip.white = TRUE) %>%
  filter(interaction_claim_source != "ChemblInteractions") %>% 
  filter(drug_claim_primary_name != "", gene_name != "") %>% 
  filter(interaction_types %in% types)

names <- as.data.frame(dgidb.targets$drug_claim_primary_name) %>% distinct()
write.table(names, "NoGit/dgidb_names.txt", sep = "\t", quote = F, row.names = F)

###export, run through pubchem id exhange, manually curate remaining with chemspider if possible

names_cur <- read.table("dgidb_curated.txt", sep = "\t", header = F, comment.char = "", quote = "\"") %>% 
  set_names("dgidb_id", "smiles") %>% 
  filter(!is.na(smiles)) %>% 
  group_by(dgidb_id) %>% 
  top_n(1) %>% 
  ungroup() %>% 
  distinct()

names_cur$dgidb_id <- gsub("\x9a", "ö", names_cur$dgidb_id)
names_cur$dgidb_id <- gsub("17_-ESTRADIOL", "17β-ESTRADIOL", names_cur$dgidb_id)
names_cur$dgidb_id <- gsub("19-Apr", "APR19", names_cur$dgidb_id)

write.table(names_cur, "dgidb_structures_curated.txt", sep = "\t", row.names = F, quote = T)
synStore(File("dgidb_structures_curated.txt", parentId = "syn11672881"))
