source('helpers.R')
library(tibble)

evo <- readRDS("Data/evotec.RDS")

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

##this long string of dplyr pipes is to combine all compounds with an identical SMILES, of which there are several
dgidb.int <- inner_join(temp2, smiles)
colnames(dgidb.int) <- c("Hugo_Gene", "Common_Name", "N_DGIDB", "Original_molecule_SMILES")

dgidb.cons <- dgidb.int %>% 
  group_by(Hugo_Gene, Original_molecule_SMILES) %>% 
  mutate(sum(N_DGIDB)) %>% 
  ungroup()

dgidb.topname <- dgidb.cons %>% 
  group_by(Original_molecule_SMILES, Common_Name) %>% 
  summarize(n()) %>% 
  ungroup() %>% 
  group_by(Original_molecule_SMILES) %>% 
  top_n(1) %>% 
  ungroup()

dgidb.topname <- dgidb.topname %>% 
  add_column("count" = c(1:nrow(dgidb.topname))) %>% 
  group_by(Original_molecule_SMILES) %>% 
  top_n(1, -count) %>% 
  ungroup() %>% 
  select(Original_molecule_SMILES, Common_Name)

dgidb.cons <- dgidb.cons %>% 
  rename("nonstandard_name" = Common_Name) %>% 
  left_join(dgidb.topname) %>%
  select(Hugo_Gene, Common_Name, `sum(N_DGIDB)`, Original_molecule_SMILES) %>% 
  distinct()

colnames(dgidb.cons) <- c("Hugo_Gene", "Common_Name", "N_DGIDB", "Original_molecule_SMILES")


##straight adding - no chemical matching
#evo2 <- bind_rows(evo, dgidb.int)
#saveRDS(evo2, "Data/evotec_dgidb.RDS")

##map DGIDB with evotec where chemical overlap exists
library(parallel)
parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
}

evo<-readRDS("Data/evotec.RDS")
fp.dgi <- parseInputFingerprint(as.character(unique(dgidb.cons$Original_molecule_SMILES)))
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

colnames(dgidb.cons) <- c("Hugo_Gene", "Common_Name2", "N_DGIDB", ".id")

#join perfect matches to evo db
sims3 <- ldply(sims2) %>% filter(sim == 1) %>% left_join(dgidb.cons) %>% select(match, Hugo_Gene, Common_Name2, N_DGIDB)
colnames(sims3) <- c("Original_molecule_SMILES", "Hugo_Gene", "Common_Name_DGIDB", "N_DGIDB")
evo <- left_join(evo, sims3, by = c("Original_molecule_SMILES","Hugo_Gene"))

#add non perfect matches to bottom as additional rows
sims4 <- ldply(sims2) %>% filter(sim < 1) %>% left_join(dgidb.cons)
sims5 <- filter(dgidb.cons, .id %in% sims4$.id)
colnames(sims5) <- c("Hugo_Gene", "Common_Name", "N_DGIDB", "Original_molecule_SMILES")
evo <- bind_rows(evo, sims5)
saveRDS(evo, "Data/evotec_dgidb.rds")

###refresh fingerprint
evo2<-readRDS("Data/evotec_dgidb.RDS")
fp.evo <- parseInputFingerprint(unique(evo2$Original_molecule_SMILES))
saveRDS(fp.evo, "Data/fpevo.rds")

###refresh confidence score
evo2<-readRDS("Data/evotec_dgidb.RDS")
evo2 <- evo2 %>% rowwise() %>% mutate("neglogActivity" = -log10(MedianActivity_nM))
evo2 <- evo2 %>% rowwise() %>% mutate("neglogActivityWeighted" = -MedianActivity_nM^10)
evo2 <- evo2 %>% rowwise() %>% mutate("Confidence_Score" = sum(N_quantitative,N_qualitative,N_DGIDB,-N_inactive,neglogActivity, na.rm=T))
evo2 <- evo2 %>% rowwise() %>% mutate("Activity_Weighted_Confidence_Score" = sum(N_quantitative/116859,N_qualitative/4912,N_DGIDB,-N_inactive,neglogActivityWeighted, na.rm=T))
saveRDS(evo2, "Data/evotec_dgidb.RDS")

###refresh common names
evo2<-readRDS("Data/evotec_dgidb.RDS")
allsmiles <- unique(evo2$Original_molecule_SMILES)
write.table(allsmiles, "Data/smiles.csv", sep = ",", quote = F, col.names = F, row.names = F)
common.names<-read.table("Data/PubChem.txt", sep = "\t", comment.char = "", quote = "")
    ##get all common names from pubchem chemical id resolver
colnames(common.names) <- c("Original_molecule_SMILES", "Common_Name")
    ##clean up to enrich for human readable names and eliminate uncommon/vendor ids (not CHEMBL, HMS, ZINC etc)
    ##find common prefixes for IDs

common.names <- common.names %>%  mutate(short = substr(Common_Name, 0, 5))
tab <- table(common.names$short) %>% as.data.frame() %>% filter(Freq>500)

common.names <- common.names %>%  mutate(short = substr(Common_Name, 0, 3))
tab <- table(common.names$short) %>% as.data.frame() %>% filter(Freq>500)


common.names.filt <- common.names %>% 
  filter(!str_detect(Common_Name, "NA")) %>% 
  filter(!str_detect(Common_Name, "\\_")) %>% 
  filter(!str_detect(Common_Name, "\\'")) %>% 
  filter(!str_detect(Common_Name, "\\,")) %>% 
  filter(!str_detect(Common_Name, "\\[")) %>% 
  filter(!str_detect(Common_Name, "\\]")) %>% 
  filter(!str_detect(Common_Name, "\\(")) %>% 
  filter(!str_detect(Common_Name, "\\)")) %>% 
  filter(!str_detect(Common_Name, "^\\d+-\\d+-\\d+$")) %>%
  filter(!str_detect(Common_Name, "&")) %>%
  filter(!str_detect(Common_Name, "#")) %>%
  filter(!str_detect(Common_Name, ":")) %>%
  filter(!str_detect(Common_Name, "^MolPort.+")) %>%
  filter(!str_detect(Common_Name, "^AKOS.+")) %>%
  filter(!str_detect(Common_Name, "^KB-.+")) %>%
  filter(!str_detect(Common_Name, "^I02-.+")) %>%
  filter(!str_detect(Common_Name, "^SPECTRUM.+")) %>%
  filter(!str_detect(Common_Name, "^BDBM\\d+")) %>%
  filter(!str_detect(Common_Name, "^CTK.+")) %>%
  filter(!str_detect(Common_Name, "^EU-\\d+")) %>%
  filter(!str_detect(Common_Name, "^NU\\d+")) %>%
  filter(!str_detect(Common_Name, "^SBI-\\d+\\..+")) %>%
  filter(!str_detect(Common_Name, "^\\d+-EP\\d+A\\d")) %>%
  filter(!str_detect(Common_Name, "^BRN \\d+")) %>%
  filter(!str_detect(Common_Name, "^AB00\\d+")) %>%
  filter(!str_detect(Common_Name, "^AJ-\\d+")) %>%
  filter(!str_detect(Common_Name, "^AK-\\d+")) %>%
  filter(!str_detect(Common_Name, "^AKOS\\d+")) %>%
  filter(!str_detect(Common_Name, "^ANW-\\d+")) %>%
  filter(!str_detect(Common_Name, "^MFCD\\d+")) %>%
  filter(!str_detect(Common_Name, "^SMR000\\d+")) %>%
  filter(!str_detect(Common_Name, "^UNII-.+")) %>%
  filter(!str_detect(Common_Name, "^EN300.+")) %>%
  filter(!str_detect(Common_Name, "^UPCMLD.+")) %>%
  filter(!str_detect(Common_Name, "^ZX-\\d+")) %>%
  filter(!str_detect(Common_Name, "^DTXSID\\d+")) %>%
  filter(!str_detect(Common_Name, "/")) %>%
  filter(!str_detect(Common_Name, "^BB 0.+")) %>%
  filter(!str_detect(Common_Name, "^BBL0.+")) %>%
  filter(!str_detect(Common_Name, "^CAS-\\d+-\\d+-\\d+")) %>%
  filter(!str_detect(Common_Name, "^CCG-\\d+")) %>%
  filter(!str_detect(Common_Name, "^LP\\d+")) %>%
  filter(!str_detect(Common_Name, "^MLS\\d+")) %>%
  filter(!str_detect(Common_Name, "^NE\\d.+")) %>%
  filter(!str_detect(Common_Name, "^NSC \\d+")) %>%
  filter(!str_detect(Common_Name, "^NSC-\\d+")) %>%
  filter(!str_detect(Common_Name, "^NSC\\d+")) %>%
  filter(!str_detect(Common_Name, "^ZINC.+")) %>%
  filter(!str_detect(Common_Name, "^Spect.+")) %>%
  filter(!str_detect(Common_Name, "^Tox21.+")) %>%
  filter(!str_detect(Common_Name, "^KBio.+")) %>%
  filter(!str_detect(Common_Name, "^Prest.+")) %>%
  filter(!str_detect(Common_Name, "^BDBM.+")) %>%
  filter(!str_detect(Common_Name, "^FT-0.+")) %>%
  filter(!str_detect(Common_Name, "^MFCD0.+")) %>%
  filter(!str_detect(Common_Name, "^EINEC.+")) %>%
  filter(!str_detect(Common_Name, "^BSPBi.+")) %>%
  filter(!str_detect(Common_Name, "^KS-00.+")) %>%
  filter(!str_detect(Common_Name, "^SR-01.+")) %>%
  filter(!str_detect(Common_Name, "^SPBio.+")) %>%
  filter(!str_detect(Common_Name, "^BSPBi.+")) %>%
  mutate(Common_Name = gsub("-Supplied by Selleck Chemicals", "", Common_Name)) %>% 
  filter(Common_Name != "") %>%  
  filter(!is.numeric(Common_Name)) %>% 
  distinct() %>% 
  select(-short) %>% 
  group_by(Common_Name) %>% 
  filter(row_number() == 1) %>% 
  ungroup()

common.names.filt$Original_molecule_SMILES <- as.character(common.names.filt$Original_molecule_SMILES)

saveRDS(common.names.filt, "Data/commname.rds")

###DrugBank database
library(XML)
drugbank<-xmlTreeParse("NoGit/full database.xml")
