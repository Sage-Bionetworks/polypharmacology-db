library(dplyr)
source('helpers.R')

temp <- read.table("Data/interactions.tsv", sep = "\t", quote = "",header = T) ##from DGIDB downloads section
temp2 <- temp %>% 
  select(entrez_gene_symbol, drug_primary_name) %>% 
  group_by(entrez_gene_symbol,drug_primary_name) %>% 
  summarize(count = n()) %>% 
  ungroup()

#clean up weird name formatting for pubchem
temp2$drug_primary_name <- gsub("<.*?>", "", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&ALPHA\\;", "ALPHA", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&BETA\\;", "BETA", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&DELTA\\;", "DELTA", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&GAMMA\\;", "GAMMA", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&OMEGA\\;", "OMEGA", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&TAU\\;", "TAU", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&KAPPA\\;", "KAPPA", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&LAMBDA\\;", "LAMBDA", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&MU\\;", "MU", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&RHO\\;", "RHO", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("\\&PSI\\;", "PSI", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("(\\&PLUSMN\\;)", "+/-", temp2$drug_primary_name)
temp2$drug_primary_name <- gsub("(\\&OUML\\;)", "o", temp2$drug_primary_name)

temp3 <- as.data.frame(unique(temp2$drug_primary_name))
write.table(temp3, "dgidb_mols.txt", sep = "\t", quote = F, col.names = F, row.names = F)
##upload this table to pubchem ID exhange to get smiles

smiles <- read.table("Data/dgidb_smiles_from_pubchem.txt", 
                     sep = "\t",
                     quote = "",
                     header = F, 
                     comment.char = "",
                     na.strings = NA)
smiles$count <- c(1:nrow(smiles))
smiles <- smiles %>% 
  group_by(V1) %>% 
  top_n(1, -count) %>% 
  select(V1,V2) %>% 
  filter(V2 != "") %>% 
  ungroup()

colnames(smiles) <- c("drug_primary_name", "smiles")

dgidb.int <- inner_join(temp2, smiles)

colnames(dgidb.int) <- c("Hugo_Gene", "Common_Name", "N_DGIDB", "Original_molecule_SMILES")


evo2 <- bind_rows(evo, dgidb.int)
saveRDS(evo2, "Data/evotec_dgidb.RDS")

##map DGIDB with evotec where chemical overlap exists
library(parallel)
evo<-readRDS("Data/evotec.RDS")
fp.dgi <- parseInputFingerprint(as.character(unique(dgidb.int$smiles)))
fp.evo <- parseInputFingerprint(as.character(unique(evo$Original_molecule_SMILES)))

sims <- mclapply(fp.dgi, function(i) { ##this will take a really long time and yield a multi gb file
  sim <- sapply(fp.evo, function(j) {
    distance(i, j)
  })
  bar <- as.data.frame(sim)
  bar$match <- rownames(bar)
  bar
}, mc.cores = detectCores())

sims2 <- mclapply(sims, function(i) {
  print(i)
  x <- top_n(i, 1, sim)
  x
}, mc.cores = detectCores())

colnames(dgidb.int) <- c("Hugo_Gene", "Common_Name2", "N_DGIDB", ".id")

#join perfect matches to evo db
sims3 <- ldply(sims2) %>% filter(sim == 1) %>% left_join(dgidb.int) %>% select(match, Hugo_Gene, Common_Name2, N_DGIDB)
colnames(sims3) <- c("Original_molecule_SMILES", "Hugo_Gene", "Common_Name_DGIDB", "N_DGIDB")
evo <- left_join(evo, sims3, by = c("Original_molecule_SMILES","Hugo_Gene"))

#add non perfect matches to bottom as additional rows
sims4 <- ldply(sims2) %>% filter(sim < 1) %>% left_join(dgidb.int)
sims5 <- filter(dgidb.int, .id %in% sims4$.id)
colnames(sims5) <- c("Hugo_Gene", "Common_Name", "N_DGIDB", "Original_molecule_SMILES")
evo <- bind_rows(evo, sims5)
saveRDS(evo, "Data/evotec_dgidb.rds")

###refresh fingerprint
evo2<-readRDS("Data/evotec_dgidb.RDS")
fp.evo <- parseInputFingerprint(unique(evo2$Original_molecule_SMILES))
saveRDS(fp.evo, "Data/fpevo.rds")
