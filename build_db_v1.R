options(java.parameters = "-Xmx12g" ) 
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

db_struct$database <- "drugbank"

##import CHeMBL v23 structure data downloaded from ftp site
chembl_struct <- read.csv(synGet("syn11672894")$path, header = T, comment.char = "", allowEscapes = FALSE) %>% distinct()
colnames(chembl_struct) <- c("external_id", "smiles")
chembl_struct$external_id <- as.factor(chembl_struct$external_id)

chembl_struct$database <- "chembl"

##import DGIDB 3.0
dgidb_struct <- read.table(synGet("syn11832827")$path, 
                           sep = "\t",
                           header = T, 
                           comment.char = "",
                           na.strings = NA,
                           allowEscapes = FALSE) %>% 
  set_names(c("external_id", "smiles"))

dgidb_struct$database <- "dgidb"

##import ChemicalProbes.org structures
cp_struct <- read.table(synGet("syn11685572")$path, 
                        sep = "\t",
                        quote = "",
                        header = T, 
                        comment.char = "",
                        na.strings = NA,
                        allowEscapes = FALSE) %>% 
  group_by(external.id) %>% 
  mutate(count = n()) %>% slice(1) %>% 
  ungroup() %>% 
  select(external.id, smiles, database) %>% 
  set_names(c("external_id", "smiles", "database")) %>% 
  mutate(source = "cp") %>% 
  unite(external_id, c("external_id", "source"), sep = "_")

##import 
klaeger_struct <- read.table(synGet("syn11685586")$path, 
                             sep = "\t",
                             quote = "",
                             header = T, 
                             comment.char = "",
                             na.strings = NA,
                             allowEscapes = FALSE) %>% 
  select(Drug, SMILES.code) %>% 
  set_names(c("external_id", "smiles")) %>% 
  mutate(database = "klaeger", source = "klaeger") %>% 
  unite(external_id, c("external_id", "source"), sep = "_")


structures <- bind_rows(chembl_struct, db_struct, dgidb_struct, cp_struct, klaeger_struct)

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

write.table(mat.df2, "NoGit/grouped_structures.txt", sep = "\t", row.names = F)
synStore(File("NoGit/grouped_structures.txt", parentId = "syn11678675"), ####"syn11678675"), change this to parentId for current build
         used = c("syn11672894","syn11673040","syn11673095","syn11685572","syn11685586"), executed = this.file)


####assemble all the bits together
grouped_structures <- read.table(synGet("syn11678713")$path, 
                                 sep = "\t", 
                                 comment.char = "",
                                 header = T,
                                 colClasses = c("character", 
                                                "character",
                                                "character",
                                                "character")) %>% 
  select(2:5) %>% 
  filter(smiles != "")

##Targets

##DrugBank 5.0.11
drugbank.targets <- read.table(synGet("syn11673549")$path, strip.white = TRUE) %>% 
  dplyr::select(drug, hugo_gene) %>% 
  filter(drug != "", hugo_gene != "") %>% 
  set_names(c("external_id", "hugo_gene")) %>% 
  group_by(external_id, hugo_gene) %>% 
  dplyr::add_tally() %>% 
  ungroup() %>% 
  distinct()

##DGIDB 3.0 ##this currently contains some likely redundant information from Chembl, remove
types <- c("inhibitor",  "antagonist", "agonist", "binder",  "modulator", "blocker", "channel blocker", 
           "positive allosteric modulator", "allosteric modulator", "activator", "inverse agonist", "partial agonist",
           "activator,channel blocker", "gating inhibitor", "agonist,antagonist", "agonist,allosteric modulator", "activator,antagonist",
           "stimulator", "negative modulator", "allosteric modulator,antagonist",  "channel blocker,gating inhibitor antagonist,inhibitor",            
           "inhibitory allosteric modulator")  

dgidb.targets <- read.table(synGet("syn11672978")$path, sep = "\t", quote = "",header = T, strip.white = TRUE) %>%
  filter(interaction_claim_source != "ChemblInteractions") %>% 
  filter(interaction_types %in% types) %>% 
  dplyr::select(drug_claim_primary_name, gene_name) %>% 
  filter(drug_claim_primary_name != "", gene_name != "") %>% 
  set_names(c("external_id", "hugo_gene")) %>% 
  group_by(external_id, hugo_gene) %>% 
  dplyr::add_tally() %>% 
  ungroup() %>% 
  distinct()


##CP Jan 2018 Targets

cp.targets <- read.table(synGet("syn11685574")$path, sep = "\t", quote = "", header = T, strip.white = TRUE) %>% 
  filter(Number.of.SAB.Reviews > 0) %>% 
  mutate(Protein.target = gsub("\"","",Protein.target), external_id = Probe.Name, source = c("cp"))%>% 
  unite(external_id, c("external_id", "source"), sep = "_") %>% 
  select(external_id, Protein.target) %>% 
  separate(Protein.target, into = c("1","2","3","4","5"), sep = ",") %>% 
  gather(key = "key", value = "hugo_gene", -external_id) %>%
  select(-key) %>% 
  filter(!is.na(hugo_gene)) %>% 
  mutate(hugo_gene = trimws(hugo_gene)) %>% 
  group_by(external_id, hugo_gene) %>% 
  add_tally() %>% 
  ungroup() %>% 
  distinct()

######### ALL QUALITATIVE DATA - 24080 Qualitative Drug-Target Associations
qual.targets <- bind_rows(drugbank.targets, dgidb.targets, cp.targets) %>% 
  left_join(x=grouped_structures, y = .) %>% 
  dplyr::select(internal_id, hugo_gene, n) %>% 
  group_by(internal_id, hugo_gene) %>% 
  summarize("n_qualitative" = sum(n)) %>% 
  ungroup() %>% 
  filter(!is.na(hugo_gene))

######### Summarize Quantitative Data
##Klaeger 2017
klaeger.drugnames <- t(read.table(synGet("syn11685587")$path, sep = "\t",
                                  header = F, comment.char = "", quote = "")[1,4:246]) %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  set_names(c("external_id")) %>% 
  mutate(Drug = make.names(external_id))

klaeger.targets <- read.table(synGet("syn11685587")$path, sep = "\t",
                              header =T, comment.char = "", quote = "", strip.white = TRUE) %>% 
  dplyr::select(-Kinase, -Direct.binder) %>% 
  gather(key = "Drug", value = "Kd", -Gene.name) %>% 
  filter(!is.na(Kd)) %>% 
  mutate(pchembl_value = -log10(Kd/1000000000)) %>%  #-log10(Kd molar)
  separate(Gene.name, into = c("1","2"), sep = ";") %>% ##some genes have 2 isoforms listed, separate and gather data
  gather(key = "key", value = "hugo_gene", -Drug, -Kd, -pchembl_value) %>% 
  select(-key) %>% 
  filter(!is.na(hugo_gene)) %>% 
  left_join(klaeger.drugnames) %>% 
  mutate(source = "klaeger") %>% 
  unite("external_id", external_id , source, sep = "_") %>% 
  left_join(grouped_structures) 

klaeger_pchembl <- klaeger.targets %>% 
  select(internal_id, hugo_gene, pchembl_value)

klaeger_assaytype_summary <- klaeger.targets %>% 
  select(internal_id, hugo_gene, Kd) %>% 
  set_names("internal_id", "hugo_gene", "standard_value") %>% 
  mutate(standard_type = "Kd")

##ChEMBL v23 - SQL query on synapse file page V
chembl.targets <- read.table(synGet("syn11672909")$path, sep = "\t",
                             header = T, comment.char = "", quote = "", strip.white = TRUE)

##filter out non-hugo genes introduced by chembl synonyms
hugo <- read.table(synGet("syn11695142")$path, header = T, sep = "\t", quote = "", comment.char = "", colClasses = "character")

chembl.targets <- dplyr::filter(chembl.targets, component_synonym %in% as.character(hugo$symbol))

pchembl.summary <- chembl.targets %>% 
  mutate(external_id = as.character(molregno)) %>% 
  left_join(x=grouped_structures, y = .) %>% 
  dplyr::select(internal_id, component_synonym, pchembl_value) %>% 
  set_names(c("internal_id", "hugo_gene", "pchembl_value")) %>% 
  filter(!is.na(pchembl_value)) %>% 
  filter(pchembl_value != "NULL") %>% 
  mutate(pchembl_value = as.numeric(as.character(pchembl_value))) %>% 
  bind_rows(klaeger_pchembl) %>% 
  group_by(internal_id, hugo_gene) %>% 
  summarize("n_quantitative" = n(), 
            "mean_pchembl" = mean(pchembl_value, na.rm = T)) %>% 
  ungroup()

chembl.assaytype.summary <- chembl.targets %>% ##add in klaeger_quant data too
  mutate(external_id = as.character(molregno)) %>% 
  left_join(x=grouped_structures, y = .) %>% 
  dplyr::select(internal_id, component_synonym, standard_value, standard_type) %>% 
  set_names("internal_id", "hugo_gene", "standard_value", "standard_type") %>% 
  filter(standard_type %in% c("IC50", "AC50", "EC50", "C50", "Potency", "Ki", "Kd", "GI50")) %>% 
  filter(!is.na(standard_type)) %>% 
  rbind(klaeger_assaytype_summary) %>% 
  group_by(internal_id, hugo_gene, standard_type) %>% 
  summarize(mean_value = mean(standard_value)) %>% 
  ungroup() %>% 
  spread(standard_type, mean_value) %>% 
  select(internal_id, hugo_gene, IC50, AC50, EC50, C50, Potency, Ki, Kd, GI50) %>% 
  set_names(c("internal_id", "hugo_gene", "IC50_nM","AC50_nM",
              "EC50_nM", "C50_nM", "Potency_nM", "Ki_nM", "Kd_nM", 
              "GI50_nM"))

quant.targets <- full_join(pchembl.summary, chembl.assaytype.summary) 

full.db <- full_join(quant.targets, qual.targets, by = c("hugo_gene", "internal_id"))



##generate database name map, trying to retain useful IDs as well as human readable names for compounds that have them
grouped_structures <- read.table(synGet("syn11678713")$path, 
                                 sep = "\t", 
                                 comment.char = "",
                                 header = T,
                                 colClasses = c("character", 
                                                "character",
                                                "character",
                                                "character")) %>% 
  select(2:5) %>% 
  filter(smiles != "")

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

cp.names <- cp_struct %>% select(external_id) %>% 
  mutate(common_name = gsub("_.+", "", external_id))

klaeger.names <- klaeger_struct %>% select(external_id) %>% 
  mutate(common_name = gsub("_.+", "", external_id))

all.names <- bind_rows(dgidb.names, db.names, chembl.names, cp.names, klaeger.names) %>% 
  inner_join(grouped_structures)


#join single synonym to db for use in app, save full table of synonyms
syns <- all.names %>% 
  group_by(internal_id) %>% 
  top_n(1) %>% 
  mutate(count = n()) %>% slice(1) %>% 
  select(internal_id, common_name)

full.db <- left_join(full.db, syns) %>% filter(!is.na(internal_id), !is.na(hugo_gene))

##add confidence metrics
full.db2 <- full.db %>% 
  group_by(internal_id) %>% 
  add_tally() %>% 
  ungroup() %>% 
  group_by(internal_id, hugo_gene) %>% 
  mutate(known_selectivity_index = 1/n) %>% 
  mutate(confidence = (sum(prod(n_qualitative,mean_pchembl,na.rm = T), n_quantitative, na.rm = T))) %>% 
  ungroup()

write.table(full.db2, "NoGit/drug_target_associations_v1.txt", row.names = F)
synStore(File("NoGit/drug_target_associations_v1.txt", parentId = "syn11678675"), executed = this.file, 
         used = c("syn11672909", "syn11672978", "syn11673549", "syn11678713"))

saveRDS(full.db2, "NoGit/drug_target_associations_v1.rds")
synStore(File("NoGit/drug_target_associations_v1.rds", parentId = "syn11678675"), executed = this.file, 
         used = c("syn11672909", "syn11672978", "syn11673549", "syn11678713"))

# library(fst)
# write_fst(full.db2,"NoGit/drug_target_associations_v1.fst")
# synStore(File("NoGit/drug_target_associations_v1.fst", parentId = "syn11678675"), executed = this.file, 
#          used = c("syn11672909", "syn11672978", "syn11673549", "syn11678713"))

write.table(all.names, "NoGit/compound_names.txt", sep = '\t', row.names = F)
synStore(File("NoGit/compound_names.txt", parentId = "syn11678675"), executed = this.file, 
         used = c("syn11673040", "syn11681825", "syn11672978"))

saveRDS(all.names, "NoGit/compound_names.rds")
synStore(File("NoGit/compound_names.rds", parentId = "syn11678675"), executed = this.file, 
         used = c("syn11673040", "syn11681825", "syn11672978"))

saveRDS(all.names, "NoGit/compound_names.fst")
synStore(File("NoGit/compound_names.fst", parentId = "syn11678675"), executed = this.file, 
         used = c("syn11673040", "syn11681825", "syn11672978"))


####Generate fingerprints for database.

structures <- read.table(synGet("syn11678713")$path, header = T) 

structures.distinct <- structures %>% 
  group_by(internal_id) %>% 
  top_n(1) %>% 
  mutate(count = n()) %>% 
  slice(1) %>% 
  select(-fp, -external_id, -database)

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
  pblapply(input.mol, get.fingerprint, type = "extended")
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

names(foo) <- structures.distinct$internal_id

saveRDS(foo, "NoGit/db_fingerprints.rds")
synStore(File("NoGit/db_fingerprints.rds", parentId = "syn11678675"), used = c("syn11678713"), executed = this.file)

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

names(foo) <- structures.distinct$internal_id

saveRDS(foo, "NoGit/db_fingerprints_circular.rds")
synStore(File("NoGit/db_fingerprints_circular.rds", parentId = "syn11678675"), used = c("syn11678713"), executed = this.file)

#### create igraph object from db
library(igraph)
db <- readRDS(synGet("syn11712148")$path) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(internal_id, hugo_gene, mean_pchembl, n_quantitative, n_qualitative) %>% 
  group_by(internal_id, hugo_gene) %>% 
  mutate(total_n = sum(n_quantitative, n_qualitative, na.rm = T)) %>% 
  ungroup() 

db.names <- readRDS(synGet("syn11712148")$path) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(internal_id, common_name) %>% 
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

names <- readRDS(synapser::synGet("syn11802195")$path)
write.table(names, "NoGit/drug-target_explorer_igraph_name_map.txt", sep = "\t", row.names = F)

results <- synapser::synTableQuery(sprintf("select * from %s", "syn11831632"))
x <- nrow(results$asDataFrame())/10000
for(i in 1:ceiling(x)){
  print(i)
  results <- synapser::synTableQuery(sprintf("select * from %s limit 10000", "syn11831632"))
  deleted <- synDelete(results$asRowSet())
}

results <- synapser::synGet("syn11831632")
tableToAppend <- Table(results, names)
table <- synStore(tableToAppend)
