##options(java.parameters = "-Xmx36g" ) 
##set above parameter before loading rJavaif generating fingerprints
##should be on c4.8xlarge or equiv

library(plyr)
library(tidyverse)
library(rJava)
library(rcdk)
library(fingerprint)
library(parallel)
library(synapser)
library(rlist)
library(pbapply)
synLogin()

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/polypharmacology-db/master/build_db_v1.R"

fp.to.simple.matrix <- function( fplist ) {
  size <- fplist[[1]]@nbit
  m <- matrix(0, nrow=length(fplist), ncol=2)
  cnt <- 1
  for ( i in fplist ) {
    foo <- c(rep(0, size))
    foo[i@bits] <- 1
    m[cnt,1] <- paste0(foo, collapse = "")
    m[cnt,2] <- names(fplist)[cnt]
    cnt <- cnt + 1
  }
  m
}

is.smiles <- function(x, verbose = TRUE) { ##corrected version from webchem
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("rcdk needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # x <- 'Clc(c(Cl)c(Cl)c1C(=O)O)c(Cl)c1Cl'
  if (length(x) > 1) {
    stop('Cannot handle multiple input strings.')
  }
  out <- try(rcdk::parse.smiles(x), silent = TRUE)
  if (inherits(out[[1]], "try-error") | is.null(out[[1]])) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

###CODE for DB build using ChEMBLv23, Drugbank 5.0.11, DGIDB 3.0

##import drugbank 5.0.11 structure data 
x <- read.csv(synGet("syn11673040")$path, header = T, comment.char = "", sep = ",", allowEscapes = FALSE)
db_struct <- x %>% dplyr::select(DrugBank.ID, SMILES) %>% 
  set_names(c("external_id", "smiles")) %>% 
  filter(smiles != "")

db_struct$db <- "drugbank"

##import CHeMBL v23 structure data downloaded from ftp site
chembl_struct <- read.csv(synGet("syn11672894")$path, header = T, comment.char = "", allowEscapes = FALSE) %>% distinct()
colnames(chembl_struct) <- c("external_id", "smiles")
chembl_struct$external_id <- as.factor(chembl_struct$external_id)

chembl_struct$db <- "chembl"

##import DGIDB 3.0
dgidb_struct <- read.table(synGet("syn11673095")$path, 
                           sep = "\t",
                           quote = "",
                           header = F, 
                           comment.char = "",
                           na.strings = NA,
                           allowEscapes = FALSE)
dgidb_struct$count <- c(1:nrow(dgidb_struct))
dgidb_struct <- dgidb_struct %>% 
  group_by(V1) %>% 
  top_n(1, -count) %>% 
  select(1:2) %>% 
  filter(V2 != "") %>% 
  ungroup() %>% 
  set_names(c("external_id", "smiles"))

dgidb_struct$db <- "dgidb"
structures <- bind_rows(chembl_struct, db_struct, dgidb_struct)

valid.smiles<-pbsapply(structures$smiles, is.smiles)

smiles <- names(valid.smiles)
valid.smiles <- as.data.frame(valid.smiles)
valid.smiles$smiles <- smiles
valid <- valid.smiles$smiles[valid.smiles==TRUE]

parseInputFingerprint <- function(input) {
  print("parsing smiles")
  input.mol <- parse.smiles(as.character(input))
  print("typing smiles")
  pblapply(input.mol, do.typing)
  print("imprinting aromaticity")
  pblapply(input.mol, do.aromaticity)
  print("identifying isotopes")
  pblapply(input.mol, do.isotopes)
  print("generating fingerprints")
  pblapply(input.mol, get.fingerprint, type = "circular")
}

##takes a while
mat <- matrix(ncol = 2)
ct <- 1

for(i in 1:ceiling(length(valid)/5000)){
  if((length(valid)-(i*5000))>=0){
    print(ct)
    print(i*5000)
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- parseInputFingerprint(valid[ct:(i*5000)])
    mat <- rbind(mat, fp.to.simple.matrix(foo))
    ct<-ct+5000
  }else{
    print(ct)
    print(length(valid))
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- parseInputFingerprint(valid[ct:length(valid)])
    mat <- rbind(mat, fp.to.simple.matrix(foo))
  }
}

mat.df <- as.data.frame(mat)

mat.df2 <- mat.df %>% 
  set_names(c("fp", "smiles")) %>% 
  distinct() %>% ##critical for mapping because there are repeated smiles
  inner_join(structures)

mat.df2$internal_id <- group_indices(mat.df2, fp)  

write.table(mat.df2, "grouped_structures.txt", sep = "\t", row.names = F)
synStore(File("grouped_structures.txt", parentId = ""), ####"syn11678675"), change this to parentId for current build
         used = c("syn11672894","syn11673040", "syn11673095"), executed = this.file)


####assemble all the bits together
grouped_structures <- read.table(synGet("syn11678713")$path, 
                                 sep = "\t", 
                                 comment.char = "",
                                 header = T,
                                 colClasses = c("character", 
                                                "character",
                                                "character",
                                                "character")) %>% dplyr::select(2:4)

##Targets

##DrugBank 5.0.11
drugbank.targets <- read.table(synGet("syn11673549")$path) %>% 
  dplyr::select(drug, hugo_gene) %>% 
  filter(drug != "", hugo_gene != "") %>% 
  set_names(c("external_id", "hugo_gene")) %>% 
  group_by(external_id, hugo_gene) %>% 
  dplyr::add_tally() %>% 
  ungroup() %>% 
  distinct()

##DGIDB 3.0 ##this currently contains some likely redundant information from Chembl, remove
dgidb.targets <- read.table(synGet("syn11672978")$path, sep = "\t", quote = "",header = T) %>%
  filter(interaction_claim_source != "ChemblInteractions") %>% 
  dplyr::select(drug_claim_primary_name, gene_name) %>% 
  filter(drug_claim_primary_name != "", gene_name != "") %>% 
  set_names(c("external_id", "hugo_gene")) %>% 
  group_by(external_id, hugo_gene) %>% 
  dplyr::add_tally() %>% 
  ungroup() %>% 
  distinct()

######### ALL QUALITATIVE DATA - 24080 Qualitative Drug-Target Associations
qual.targets <- bind_rows(drugbank.targets, dgidb.targets) %>% 
  left_join(x=grouped_structures, y = .) %>% 
  dplyr::select(internal_id, hugo_gene, n) %>% 
  group_by(internal_id, hugo_gene) %>% 
  summarize("n_qualitative" = sum(n)) %>% 
  ungroup() %>% 
  filter(!is.na(hugo_gene))

######### Summarize Quantitative Data
##ChEMBL v23 - SQL query on synapse file page V
chembl.targets <- read.table(synGet("syn11672909")$path, sep = "\t",
                             header = T, comment.char = "", quote = "")

pchembl.summary <- chembl.targets %>% 
  mutate(external_id = as.character(molregno)) %>% 
  left_join(x=grouped_structures, y = .) %>% 
  dplyr::select(internal_id, component_synonym, pchembl_value) %>% 
  filter(pchembl_value != "NULL") %>% 
  mutate(pchembl_value = as.numeric(as.character(pchembl_value))) %>% 
  group_by(internal_id, component_synonym) %>%
  summarize("n_quantitative" = n(), 
            "mean_value" = mean(pchembl_value, na.rm = T)) %>% 
  set_names(c("internal_id", "hugo_gene", "n_quantitative", "mean_pchembl")) %>% 
  ungroup()

chembl.assaytype.summary <- chembl.targets %>% 
  mutate(external_id = as.character(molregno)) %>% 
  left_join(x=grouped_structures, y = .) %>% 
  dplyr::select(internal_id, component_synonym, standard_value, standard_type) %>% 
  filter(!is.na(standard_type)) %>% 
  group_by(internal_id, component_synonym, standard_type) %>% 
  summarize(mean_value = mean(standard_value)) %>% 
  spread(standard_type, mean_value) %>% 
  select(internal_id, component_synonym, IC50, AC50, EC50, C50, Potency, Ki, Kd, GI50) %>% 
  set_names(c("internal_id", "hugo_gene", "IC50_nM","AC50_nM",
              "EC50_nM", "C50_nM", "Potency_nM", "Ki_nM", "Kd_nM", 
              "GI50_nM")) %>% 
  ungroup()
  
all.chembl.summary <- full_join(pchembl.summary, chembl.assaytype.summary)

full.db <- full_join(all.chembl.summary, qual.targets)


##generate database name map, trying to retain useful IDs as well as human readable names for compounds that have them
grouped_structures <- read.table(synGet("syn11678713")$path, 
                                 sep = "\t", 
                                 comment.char = "",
                                 header = T,
                                 colClasses = c("character", 
                                                "character",
                                                "character",
                                                "character")) %>% dplyr::select(2:4)

dgidb.names <- read.table(synGet("syn11672978")$path, sep = "\t", quote = "", header = T) %>% 
  select(drug_claim_primary_name, drug_claim_name, drug_name, drug_chembl_id) %>% 
  distinct() %>% 
  mutate(drug_claim_primary_name2 = drug_claim_primary_name) %>% 
  gather("namesource", "name", -drug_claim_primary_name) %>% 
  select(-namesource) %>% 
  filter(!grepl("^\\d+$", name) & name != "") %>% 
  set_names(c("external_id", "common_name")) %>% 
  distinct()

db.names <- read.csv(synGet("syn11673040")$path, header = T) %>% 
  mutate(CAS.Number = paste0("CAS_",CAS.Number)) %>% 
  mutate(PubChem.Compound.ID = paste0("PubChemCID_",PubChem.Compound.ID)) %>% 
  mutate(ChEBI.ID = paste0("ChEBI_",ChEBI.ID)) %>% 
  mutate(ChemSpider.ID = paste0("ChemSpider_",ChemSpider.ID)) %>% 
  select(DrugBank.ID,Name,CAS.Number,PubChem.Compound.ID,ChEBI.ID,ChemSpider.ID) %>% 
  mutate(DrugBank.ID2 = DrugBank.ID) %>% 
  gather("namesource", "name", -DrugBank.ID) %>% 
  select(-namesource) %>% 
  set_names(c("external_id", "common_name"))

db.names$common_name[db.names$common_name=="CAS_NA"] <- NA
db.names$common_name[db.names$common_name=="PubChemCID_NA"] <- NA
db.names$common_name[db.names$common_name=="ChEBI_NA"] <- NA
db.names$common_name[db.names$common_name=="ChemSpider_NA"] <- NA

db.names <- filter(db.names, !is.na(common_name)) %>% distinct()

chembl.names <- read.table(synGet("syn11681825")$path, 
                           header = T, 
                           sep = "\t",
                           comment.char = "", 
                           quote = "", 
                           colClasses = c("character","character","character","character")) %>% 
  gather("namesource", "name", -molregno) %>% 
  select(-namesource) %>% 
  set_names(c("external_id", "common_name")) %>% 
  distinct() %>%
  filter(common_name != "NULL") %>% 
  mutate(external_id = as.character(external_id))

all.names <- bind_rows(dgidb.names, db.names, chembl.names) %>% 
  inner_join(grouped_structures)


#join single synonym to db for use in app, save full table of synonyms
syns <- all.names %>% 
  group_by(internal_id) %>% 
  top_n(1) %>% 
  mutate(count = n()) %>% slice(1) %>% 
  select(internal_id, common_name)

View(table(syns$common_name))

full.db <- left_join(full.db, syns)

write.table(full.db, "drug_target_associations_v1.txt", row.names = F)
synStore(File("drug_target_associations_v1.txt", parentId = "syn11678675"), executed = this.file, 
         used = c("syn11672909", "syn11672978", "syn11673549", "syn11678713"))

write.table(all.names, "NoGit/compound_names.txt", sep = '\t', row.names = F)
synStore(File("NoGit/compound_names.txt", parentId = "syn11678675"), executed = this.file, 
         used = c("syn11673040", "syn11681825", "syn11672978"))




####Generate fingerprints for database.

structures <- read.table(synGet("syn11678713")$path, header = T) 

structures.distinct <- structures %>% 
  group_by(internal_id) %>% 
  top_n(1) %>% 
  mutate(count = n()) %>% slice(1)

valid <- as.character(structures.distinct$smiles)

parseInputFingerprint <- function(input) {
  print("parsing smiles")
  input.mol <- parse.smiles(as.character(input))
  print("typing smiles")
  pblapply(input.mol, do.typing)
  print("imprinting aromaticity")
  pblapply(input.mol, do.aromaticity)
  print("identifying isotopes")
  pblapply(input.mol, do.isotopes)
  print("generating fingerprints")
  pblapply(input.mol, get.fingerprint, type = "circular")
}

foo <- list()
ct <- 1

for(i in 1:ceiling(length(valid)/5000)){
  if((length(valid)-(i*5000))>=0){
    print(ct)
    print(i*5000)
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, parseInputFingerprint(valid[ct:(i*5000)]))
    ct<-ct+5000
  }else{
    print(ct)
    print(length(valid))
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, parseInputFingerprint(valid[ct:length(valid)]))
  }
}

saveRDS(foo, "db_fingerprints.rds")
synStore(File("db_fingerprints.rds", parentId = "syn11678675"), used = c("syn11678713"), executed = this.file)
