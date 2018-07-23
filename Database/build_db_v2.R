options(java.parameters = "-Xmx8g" ) 
##set above parameter before loading rJavaif generating fingerprints
##should be on c4.8xlarge or equiv
library(plyr)
library(tidyverse)
library(rJava)
library(rcdk)
library(fingerprint)
library(parallel)
library(rlist)
library(pbapply)
library(reticulate)
use_python("/usr/local/bin/python3")
molvs <- import("molvs")
synapse <- import("synapseclient")
syn <- synapse$Synapse()
syn$login()
library(pbmcapply)


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

###CODE for DB build using ChEMBLv24.1, Drugbank 5.1.0, DGIDB 3.0.2

##import drugbank 5.1.0 structure data 
x <- read.csv(syn$get("syn12973257")$path, header = T, comment.char = "", sep = ",")
db_struct <- x %>% dplyr::select(DrugBank.ID, SMILES) %>% 
  set_names(c("external_id", "original_smiles")) %>% 
  filter(original_smiles != "")

db_struct$database <- "drugbank"

##import CHeMBL v24.1 structure data and chembl targets to limit structure space to the relevant mols
chembl.targets <- read.table(syn$get("syn12973247")$path, sep = "\t",
                             header = T, comment.char = "", quote = "\"", strip.white = TRUE)

chembl_struct <- read.table(syn$get("syn12973248")$path, sep = "\t", header = T, comment.char = "") %>% 
  select(1,4) %>% 
  distinct() %>% 
  filter(molregno %in% chembl.targets$molregno)

colnames(chembl_struct) <- c("external_id", "original_smiles")
chembl_struct$external_id <- as.factor(chembl_struct$external_id)

chembl_struct$database <- "chembl"

##import DGIDB 3.0.2
dgidb_struct <- read.table(syn$get("syn14721555")$path, 
                           sep = "\t",
                           header = T, 
                           comment.char = "",
                           na.strings = NA,
                           allowEscapes = FALSE) %>% 
  set_names(c("external_id", "original_smiles"))

dgidb_struct$database <- "dgidb"

##import ChemicalProbes.org structures
cp_struct <- read.csv(syn$get("syn12910000")$path, 
                      quote = "\"",
                      header = T, 
                      comment.char = "",
                      na.strings = NA,
                      allowEscapes = FALSE) %>% 
  group_by(external.id) %>% 
  mutate(count = n()) %>% slice(1) %>% 
  ungroup() %>% 
  select(external.id, smiles, database) %>% 
  set_names(c("external_id", "original_smiles", "database")) %>% 
  mutate(source = "cp") %>% 
  unite(external_id, c("external_id", "source"), sep = "_")

##import 
klaeger_struct <- read.table(syn$get("syn11685586")$path, 
                             sep = "\t",
                             quote = "",
                             header = T, 
                             comment.char = "",
                             na.strings = NA,
                             allowEscapes = FALSE) %>% 
  select(Drug, SMILES.code) %>% 
  set_names(c("external_id", "original_smiles")) %>% 
  mutate(database = "klaeger", source = "klaeger") %>% 
  unite(external_id, c("external_id", "source"), sep = "_")

structures <- bind_rows(chembl_struct, db_struct, dgidb_struct, cp_struct, klaeger_struct)

parser <- get.smiles.parser()
valid.smiles<-pbsapply(structures$original_smiles, is.smiles)

smiles <- names(valid.smiles)
valid.smiles <- as.data.frame(valid.smiles)
valid.smiles$original_smiles <- smiles
valid.smiles <- distinct(valid.smiles)

##remove structures that molvs cannot handle in standardization
##these were simply captured by running standardized_smiles iteratively with a print() inserted...
##better solution for future. not sure what the common theme between these mols is? they are all parseable

# write.molecules(list, "Data/db_prelim.sdf")
# list2 <- load.molecules("Data/db_prelim.sdf")
#' #'in order to standardize molecules, use library(camb)
#' #'However, due to recent OS updates, this can be challenging to install on macOS
#' #'and for us, required installation on an ubuntu EC2 instance
#' #'we used ami-fd2ffe87, built by Louis Aslett
#' #'In addition, the "StandardiseMolecules" function causes RStudio Server to crash,
#' #'so we recommend writing out molecule sdf file (as above) and then 
#' #'running the following two lines using command-line R, rather than in RStudio Server.
#' 
#' Note Jul 6 2018 - 
#' 
#' This approach does not appear to work on large SDFs (such one that includes all of the molecules in this database)
#' The dev for camb is looking into this, but for now I am using reticulate + python package molvs instead
#'
#' 
#' library(camb)
#' StandardiseMolecules("db_prelim.sdf", "db_standardized.sdf", removed.file = "db_removed.sdf")
#' 
#' 


#'after running this function above with base R, we can read the molecules using RStudio 
#'and continue building the db
#'
#'mols <- load.molecules(molfiles="Data/db_standardized.sdf", aromaticity = TRUE, typing = TRUE, isotopes = TRUE,
#'verbose=FALSE)
#
valid <- valid.smiles$original_smiles[valid.smiles$valid.smiles==TRUE]

standardized_smiles <- pbmclapply(unique(valid), function(x){
  tryCatch({
    molvs$standardize_smiles(r_to_py(x))
  }, warning = function(w) {
    print(paste("structure",x,"has an issue"))
  }, error = function(e) {
    print(paste("structure",x,"is invalid"))
  })
})

names(standardized_smiles) <- unique(valid)

valid.df <- standardized_smiles %>% 
  ldply() %>%
  set_names(c("original_smiles", "smiles"))

valid.smiles.2<-pbsapply(valid.df$smiles, is.smiles)

smiles <- names(valid.smiles.2)
valid.smiles.2 <- as.data.frame(valid.smiles.2)
valid.smiles.2$smiles <- smiles
valid.smiles.2 <- distinct(valid.smiles.2)

valid.df <- valid.df %>% 
  filter(!grepl("structure .+ is invalid", smiles)) %>% 
  filter(smiles %in% valid.smiles.2$smiles[valid.smiles.2$valid.smiles.2==TRUE]) %>% 
  distinct()

valid.df$internal_id <- group_indices(valid.df, smiles) 

# valid.std <- unique(valid.df$smiles)
# 
 parseInputFingerprint <- function(input, type) {
   print("parsing smiles")
   input.mol <- parse.smiles(as.character(input))
   print("doing typing")
   pblapply(input.mol, do.typing)
   print("doing aromaticity")
   pblapply(input.mol, do.aromaticity)
   print("doing isotopes")
   pblapply(input.mol, do.isotopes)
   print("generating fingerprints")
   pblapply(input.mol, get.fingerprint, type = type)
 }
# 
# ##takes a while
# mat <- matrix(ncol = 2)
# ct <- 1
# parser <- get.smiles.parser()
# 
# ##this loop parses fingerprints and converts them to strings to facilitate sorting later
# for(i in 1:ceiling(length(valid.std)/5000)){ 
#   if((length(valid.std)-(i*5000))>=0){
#     print(ct)
#     print(i*5000)
#     print(paste0("batch ", i," of ", ceiling(length(valid.std)/5000)))
#     foo <- parseInputFingerprint(valid.std[ct:(i*5000)], type = "circular")
#     mat <- rbind(mat, fp.to.simple.matrix(foo))
#     ct<-ct+5000
#   }else{
#     print(ct)
#     print(length(valid.std))
#     print(paste0("batch ", i," of ", ceiling(length(valid.std)/5000)))
#     foo <- parseInputFingerprint(valid.std[ct:length(valid.std)], type = "circular")
#     mat <- rbind(mat, fp.to.simple.matrix(foo))
#   }
# } 

##the above fingerprint generation is now done later in this script - will make this file MUCH smaller! 

# mat.df <- as.data.frame(mat)

# mat.df2 <- mat.df %>% 
#   set_names(c("fp", "smiles")) %>% 
#   distinct() %>% ##critical for mapping because there are repeated smiles
#   inner_join(valid.df) %>% 
#   inner_join(structures)

mat.df2 <- inner_join(valid.df, structures)

write.table(mat.df2, "Data/grouped_structures.txt", sep = "\t", row.names = F)
syn$store(synapse$File("Data/grouped_structures.txt", parentId = "syn12978846"), ####"syn11678675"), change this to parentId for current build
          used = c("syn12973257","syn12973247","syn12973248","syn11832827","syn12910000","syn11685586"), executed = this.file)


####assemble all the bits together
grouped_structures <- read.table(syn$get("syn12978848")$path, 
                                 sep = "\t", 
                                 comment.char = "",
                                 header = T,
                                 colClasses = c("character", 
                                                "character",
                                                "character",
                                                 "character")) %>% 
  # select(2:6) %>%
  filter(smiles != "")

##Targets

##DrugBank 5.0.11
drugbank.targets <- read.table(syn$get("syn12913659")$path, strip.white = TRUE) %>% 
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

dgidb.targets <- read.table(syn$get("syn12684108")$path, sep = "\t", quote = "",header = T, strip.white = TRUE) %>%
  filter(interaction_claim_source != "ChemblInteractions") %>% 
  filter(interaction_types %in% types) %>% 
  dplyr::select(drug_name, gene_name) %>% 
  filter(drug_name != "", gene_name != "") %>% 
  set_names(c("external_id", "hugo_gene")) %>% 
  group_by(external_id, hugo_gene) %>% 
  dplyr::add_tally() %>% 
  ungroup() %>% 
  distinct()


##CP Jan 2018 Targets

cp.targets <- read.csv(syn$get("syn12683504")$path, header = T, strip.white = TRUE) %>% 
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
klaeger.drugnames <- t(read.table(syn$get("syn11685587")$path, sep = "\t",
                                  header = F, comment.char = "", quote = "")[1,4:246]) %>% 
  as.data.frame() %>% 
  remove_rownames() %>% 
  set_names(c("external_id")) %>% 
  mutate(Drug = make.names(external_id))

klaeger.targets <- read.table(syn$get("syn11685587")$path, sep = "\t",
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
chembl.targets <- read.table(syn$get("syn12973247")$path, sep = "\t",
                             header = T, comment.char = "", quote = "\"", strip.white = TRUE)

##filter out non-hugo genes introduced by chembl synonyms
hugo <- read.table(syn$get("syn11695142")$path, header = T, sep = "\t", quote = "", comment.char = "", colClasses = "character") 

chembl.targets <- dplyr::filter(chembl.targets, component_synonym %in% as.character(hugo$symbol))%>% 
  filter(pchembl_value != "NULL") #remove null pchembl values (small percentage)

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
            "mean_pchembl" = mean(pchembl_value, na.rm = T),
            "cv" = raster::cv(pchembl_value),
            "sd" = sd(pchembl_value)) %>% 
  ungroup()

chembl.assaytype.summary <- chembl.targets %>% ##add in klaeger_quant data too
  mutate(external_id = as.character(molregno)) %>% 
  left_join(x=grouped_structures, y = .) %>% 
  dplyr::select(internal_id, component_synonym, standard_value, standard_type) %>% 
  set_names("internal_id", "hugo_gene", "standard_value", "standard_type") %>% 
  filter(standard_type %in% c("IC50", "AC50", "EC50", "C50", "Potency", "Ki", "Kd")) %>% 
  filter(!is.na(standard_type)) %>% 
  rbind(klaeger_assaytype_summary) %>% 
  group_by(internal_id, hugo_gene, standard_type) %>% 
  summarize(mean_value = mean(standard_value)) %>%
  ungroup() %>% 
  spread(standard_type, mean_value) %>% 
  select(internal_id, hugo_gene, IC50, AC50, EC50, Potency, Ki, Kd) %>% 
  set_names(c("internal_id", "hugo_gene", "IC50_nM","AC50_nM",
              "EC50_nM", "Potency_nM", "Ki_nM", "Kd_nM"))

quant.targets <- full_join(pchembl.summary, chembl.assaytype.summary)


full.db <- full_join(quant.targets, qual.targets, by = c("hugo_gene", "internal_id"))



##generate database name map, trying to retain useful IDs as well as human readable names for compounds that have them
grouped_structures <- read.table(syn$get("syn12978848")$path, 
                                 sep = "\t", 
                                 comment.char = "",
                                 allowEscapes = F,
                                 header = T,
                                 colClasses = c("character", 
                                                "character",
                                                "character",
                                                "character")) %>% 
  # dplyr::select(2:6) %>% 
  filter(smiles != "")

dgidb.names <- read.table(syn$get("syn12684108")$path, sep = "\t", quote = "", header = T) %>% 
  select(drug_claim_primary_name, drug_claim_name, drug_name, drug_chembl_id) %>% 
  distinct() %>% 
  mutate(drug_claim_primary_name2 = drug_claim_primary_name) %>% 
  gather("namesource", "name", -drug_claim_primary_name) %>% 
  select(-namesource) %>% 
  filter(!grepl("^\\d+$", name) & name != "") %>% 
  set_names(c("external_id", "common_name")) %>% 
  distinct()

db.names <- read.csv(syn$get("syn12973257")$path, header = T) %>% 
  mutate(CAS.Number = paste0("CAS_",CAS.Number)) %>% 
  mutate(PubChem.Compound.ID = paste0("PubChemCID_",PubChem.Compound.ID)) %>% 
  mutate(ChEBI.ID = paste0("ChEBI_",ChEBI.ID)) %>% 
  mutate(ChemSpider.ID = paste0("ChemSpider_",ChemSpider.ID)) %>% 
  select(DrugBank.ID,Name,CAS.Number,PubChem.Compound.ID,ChEBI.ID,ChemSpider.ID) %>% 
  mutate(DrugBank.ID2 = DrugBank.ID) %>% 
  gather("namesource", "name", -DrugBank.ID) %>% 
  select(-namesource) %>% 
  set_names(c("external_id", "common_name"))

db.names <- db.names %>%
  filter(!grepl("_$",db.names$common_name) | !grepl("_NA$",db.names$common_name)) %>% 
  distinct()

chembl.names <- read.table(syn$get("syn12972665")$path, 
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

grouped_structures <- mutate(grouped_structures, external_id = as.character(external_id))

all.names <- bind_rows(dgidb.names, db.names, chembl.names, cp.names, klaeger.names) %>% 
  full_join(grouped_structures) %>% 
  filter(!is.na(internal_id)) %>% 
  filter(!common_name %in% c("ChEBI_NA", "CAS_", "PubChemCID_NA","ChemSpider_NA"))

##drug_claim_primary_names are error prone in dgidb, use external ID instead
all.names$common_name[all.names$database=="dgidb"] <- all.names$external_id[all.names$database=="dgidb"]

syns <- all.names %>% 
  distinct() %>% 
  select(common_name, smiles, internal_id) %>% 
  mutate(tolowername = tolower(common_name)) %>% 
  group_by(tolowername, internal_id) %>% 
  add_tally() %>% 
  ungroup() %>% 
  group_by(internal_id) %>% 
  top_n(1, n) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(internal_id, common_name)

single.names <- all.names %>% 
  distinct() %>% 
  select(common_name, smiles, internal_id) %>% 
  mutate(tolowername = tolower(common_name)) %>% 
  group_by(tolowername, internal_id) %>% 
  add_tally() %>% 
  ungroup() %>% 
  group_by(tolowername) %>% 
  slice(1) %>%
  ungroup() %>% 
  select(-tolowername, -n)

#join single synonym to db for use in app, save full table of synonyms


full.db <- left_join(full.db, syns) %>% filter(!is.na(internal_id), !is.na(hugo_gene))

##add confidence metrics, round off data to make table much smaller

total_qualitative <- sum(full.db$n_qualitative, na.rm = T)
total_quantitative <- sum(full.db$n_quantitative, na.rm = T)

full.db2 <- full.db %>% 
  group_by(internal_id, hugo_gene) %>% 
  mutate(total_n = sum(n_quantitative,n_qualitative, na.rm=T)) %>% 
  ungroup()

sd.n <- sd(full.db2$total_n, na.rm = T)
mean.n <- mean(full.db2$total_n, na.rm = T)

full.db2 <- full.db2 %>% 
  mutate(confidence = (total_n-mean.n)/sd.n) %>% 
  group_by(internal_id) %>% 
  mutate(pchembl_d = sum(mean_pchembl, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(hugo_gene) %>% 
  mutate(pchembl_t = sum(mean_pchembl, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(internal_id, hugo_gene) %>% 
  # mutate(ksi_dt = mean_pchembl/((pchembl_t+pchembl_d)-mean_pchembl)) %>%  ##normalizes for target frequency, perhaps not useful for "drug" selectivity
  mutate(known_selectivity_index = mean_pchembl/pchembl_d) %>% 
  mutate(confidence = signif(confidence,3)) %>% 
  mutate(mean_pchembl = signif(mean_pchembl,3)) %>% 
  mutate(known_selectivity_index = signif(known_selectivity_index,3)) %>% 
  ungroup()

write.table(full.db2, "Data/drug_target_associations_v2.txt", row.names = F)
syn$store(synapse$File("Data/drug_target_associations_v2.txt", parentId = "syn12978846"), executed = this.file, 
          used = c("syn12978848","syn12913659","syn12684108","syn12683504","syn11685587","syn11695142","syn12973247",
                   "syn12684108","syn12973257","syn12972665"))

saveRDS(full.db2, "Data/drug_target_associations_v2.rds")
syn$store(synapse$File("Data/drug_target_associations_v2.rds", parentId = "syn12978846"), executed = this.file, 
          used = c("syn12978848","syn12913659","syn12684108","syn12683504","syn11685587","syn11695142","syn12973247",
                   "syn12684108","syn12973257","syn12972665"))

# library(fst)
# write_fst(full.db2,"Data/drug_target_associations_v1.fst")
# syn$store(synapse$File("Data/drug_target_associations_v1.fst", parentId = "syn11678675"), executed = this.file, 
#          used = c("syn11672909", "syn11672978", "syn11673549", "syn11678713"))


write.table(all.names, "Data/compound_names.txt", sep = '\t', row.names = F)
syn$store(synapse$File("Data/compound_names.txt", parentId = "syn12978846"), executed = this.file, 
          used = c("syn12978848","syn12684108","syn12973257","syn12972665","syn11685586","syn12910000"))

saveRDS(all.names, "Data/compound_names.rds")
syn$store(synapse$File("Data/compound_names.rds", parentId = "syn12978846"), executed = this.file, 
          used = c("syn12978848","syn12684108","syn12973257","syn12972665","syn11685586","syn12910000"))

write.table(single.names, "Data/distinct_compound_names.txt", sep = '\t', row.names = F)
syn$store(synapse$File("Data/distinct_compound_names.txt", parentId = "syn12978846"), executed = this.file, 
          used = c("syn12978848","syn12684108","syn12973257","syn12972665","syn11685586","syn12910000"))

saveRDS(single.names, "Data/distinct_compound_names.rds")
syn$store(synapse$File("Data/distinct_compound_names.rds", parentId = "syn12978846"), executed = this.file, 
          used = c("syn12978848","syn12684108","syn12973257","syn12972665","syn11685586","syn12910000"))

# write_fst(all.names, "Data/compound_names.fst")
# syn$store(synapse$File("Data/compound_names.fst", parentId = "syn11678675"), executed = this.file, 
#          used = c("syn11673040", "syn11681825", "syn11672978"))

####Generate fingerprints for database.

full.db2 <- readRDS(syn$get("syn12978910")$path)

structures <- read.table(syn$get("syn12978848")$path, 
                                 sep = "\t", 
                                 comment.char = "",
                                 allowEscapes = F,
                                 header = T,
                                 colClasses = c("character", 
                                                "character",
                                                "character",
                                                "character")) %>% 
  # dplyr::select(2:6) %>% 
  filter(smiles != "")

structures.distinct <- structures %>% 
  filter(internal_id %in% full.db2$internal_id) %>% 
  group_by(internal_id) %>% 
  select(-external_id, -database, -original_smiles) %>% 
  distinct()

valid <- as.character(structures.distinct$smiles)

parseInputFingerprint <- function(input, type) {
  print("parsing smiles")
  input.mol <- parse.smiles(as.character(input))
  print("doing typing")
  pblapply(input.mol, do.typing)
  print("doing aromaticity")
  pblapply(input.mol, do.aromaticity)
  print("doing isotopes")
  pblapply(input.mol, do.isotopes)
  print("generating fingerprints")
  pblapply(input.mol, get.fingerprint, type = type)
}

parser <- get.smiles.parser()
foo <- list()
ct <- 1

for(i in 1:ceiling(length(valid)/5000)){
  if((length(valid)-(i*5000))>=0){
    print(ct)
    print(i*5000)
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, parseInputFingerprint(valid[ct:(i*5000)], type = "extended"))
    ct<-ct+5000
  }else{
    print(ct)
    print(length(valid))
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, parseInputFingerprint(valid[ct:length(valid)], type = "extended"))
  }
}

names(foo) <- structures.distinct$internal_id

saveRDS(foo, "Data/db_fingerprints_extended.rds")
syn$store(synapse$File("Data/db_fingerprints_extended.rds", parentId = "syn12978846"), used = c("syn12978910","syn12978848"), executed = this.file)

foo <- list()
ct <- 1

for(i in 1:ceiling(length(valid)/5000)){
  if((length(valid)-(i*5000))>=0){
    print(ct)
    print(i*5000)
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, parseInputFingerprint(valid[ct:(i*5000)], type = "circular"))
    ct<-ct+5000
  }else{
    print(ct)
    print(length(valid))
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, parseInputFingerprint(valid[ct:length(valid)], type = "circular"))
  }
}

names(foo) <- structures.distinct$internal_id

saveRDS(foo, "Data/db_fingerprints_circular.rds")
syn$store(synapse$File("Data/db_fingerprints_circular.rds", parentId = "syn12978846"), used = c("syn12978910","syn12978848"), executed = this.file)

# foo <- list()
# ct <- 1

# for(i in 1:ceiling(length(valid)/5000)){ ##note - kr fingerprints take a while to generate (nearly 24hr on laptop for 282k mol)
#   if((length(valid)-(i*5000))>=0){
#     print(ct)
#     print(i*5000)
#     print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
#     foo <- append(foo, parseInputFingerprint(valid[ct:(i*5000)], type = "kr"))
#     ct<-ct+5000
#   }else{
#     print(ct)
#     print(length(valid))
#     print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
#     foo <- append(foo, parseInputFingerprint(valid[ct:length(valid)], type = "kr"))
#   }
# }
# 
# names(foo) <- structures.distinct$internal_id
# 
# saveRDS(foo, "Data/db_fingerprints_kr.rds")
# syn$store(synapse$File("Data/db_fingerprints_kr.rds", parentId = "syn11678675"), used = c("syn11678713"), executed = this.file)
# 
foo <- list()
ct <- 1

for(i in 1:ceiling(length(valid)/5000)){
  if((length(valid)-(i*5000))>=0){
    print(ct)
    print(i*5000)
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, parseInputFingerprint(valid[ct:(i*5000)], type = "maccs"))
    ct<-ct+5000
  }else{
    print(ct)
    print(length(valid))
    print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
    foo <- append(foo, parseInputFingerprint(valid[ct:length(valid)], type = "maccs"))
  }
}

names(foo) <- structures.distinct$internal_id

saveRDS(foo, "Data/db_fingerprints_maccs.rds")
syn$store(synapse$File("Data/db_fingerprints_maccs.rds", parentId = "syn12978846"), used = c("syn12978910","syn12978848"), executed = this.file)

# foo <- list()
# ct <- 1
# 
# for(i in 1:ceiling(length(valid)/5000)){
#   if((length(valid)-(i*5000))>=0){
#     print(ct)
#     print(i*5000)
#     print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
#     foo <- append(foo, parseInputFingerprint(valid[ct:(i*5000)], type = "pubchem"))
#     ct<-ct+5000
#   }else{
#     print(ct)
#     print(length(valid))
#     print(paste0("batch ", i," of ", ceiling(length(valid)/5000)))
#     foo <- append(foo, parseInputFingerprint(valid[ct:length(valid)], type = "pubchem"))
#   }
# }
# 
# names(foo) <- structures.distinct$internal_id
# 
# saveRDS(foo, "Data/db_fingerprints_pubchem.rds")
# syn$store(synapse$File("Data/db_fingerprints_pubchem.rds", parentId = "syn11678675"), used = c("syn11678713"), executed = this.file)


#### create igraph object from db 
###for fendr 
# library(igraph)
# db <- readRDS(syn$get("syn11712148")$path) %>% 
#   filter(!is.na(hugo_gene)) %>% 
#   select(internal_id, hugo_gene, mean_pchembl, n_quantitative, n_qualitative) %>% 
#   group_by(internal_id, hugo_gene) %>% 
#   mutate(total_n = sum(n_quantitative, n_qualitative, na.rm = T)) %>% 
#   ungroup() 
# 
# db.names <- readRDS(syn$get("syn11712148")$path) %>% 
#   filter(!is.na(hugo_gene)) %>% 
#   select(internal_id, common_name) %>% 
#   distinct()
# 
# db.igraph <- graph.data.frame(db)
# 
# saveRDS(db.igraph, "drug-target_explorer_igraph.rds")
# saveRDS(db.names, "drug-target_explorer_igraph_name_map.rds")
# 
# syn$store(synapse$File("drug-target_explorer_igraph.rds", parentId = "syn11802193"),
#          used = "syn11712148",
#          executed = this.file)
# 
# syn$store(synapse$File("drug-target_explorer_igraph_name_map.rds", parentId = "syn11802193"),
#          used = "syn11712148",
#          executed = this.file)
# 
# names <- readRDS(synapser::syn$get("syn11802195")$path)
# write.table(names, "Data/drug-target_explorer_igraph_name_map.txt", sep = "\t", row.names = F)
# 
# results <- synapser::synTableQuery(sprintf("select * from %s", "syn11831632"))
# x <- nrow(results$asDataFrame())/10000
# for(i in 1:ceiling(x)){
#   print(i)
#   results <- synapser::synTableQuery(sprintf("select * from %s limit 10000", "syn11831632"))
#   deleted <- synDelete(results$asRowSet())
# }
# 
# results <- synapser::syn$get("syn11831632")
# tableToAppend <- Table(results, names)
# table <- syn$store(tableToAppend)


###external links file
db.names <- readRDS(syn$get("syn12978912")$path)

chembl.internal.ids <- db.names %>% 
  filter(database == "chembl") %>%
  select(external_id, internal_id, database) %>% 
  distinct()

chembl.links <- read.table(syn$get("syn12972665")$path, sep = "\t", quote = "", comment.char = "", header = T) %>% 
  select(molregno, chembl_id) %>% 
  distinct() %>% 
  set_names(c("external_id", "link")) %>% 
  mutate(external_id = as.character(external_id)) %>% 
  right_join(chembl.internal.ids) %>% 
  mutate(link = paste0("<a href='https://www.ebi.ac.uk/chembl/compound/inspect/",link,"', target = '_blank'>",database,"</a>")) %>% 
  select(internal_id, link, external_id, database) %>% 
  distinct()

dgidb.links <- read.table(syn$get("syn12684108")$path, sep = "\t", quote = "", header = T) %>% 
  select(drug_claim_primary_name, drug_claim_name, drug_name, drug_chembl_id) %>% 
  distinct() %>% 
  mutate(drug_claim_primary_name2 = drug_claim_primary_name) %>% 
  gather("namesource", "name", -drug_claim_primary_name) %>% 
  select(-namesource) %>% 
  filter(!grepl("^\\d+$", name) & name != "") %>% 
  set_names(c("external_id", "common_name")) %>% 
  distinct() %>% 
  left_join(db.names) %>% 
  mutate(link = sapply(external_id, function(x) URLencode(x))) %>%
  mutate(link = paste0("<a href='http://www.dgidb.org/interaction_search_results/",link,"', target = '_blank'>",database,"</a>")) %>% 
  filter(!is.na(internal_id)) %>% 
  select(internal_id, link, external_id, database) %>% 
  distinct()

drugbank.links <- db.names %>% 
  filter(database == "drugbank") %>% 
  mutate(link = sapply(external_id, function(x) URLencode(x))) %>%
  mutate(link = paste0("<a href='https://www.drugbank.ca/drugs/",link,"', target = '_blank'>",database,"</a>")) %>% 
  select(internal_id, link, external_id, database) %>% 
  distinct()

db.links <- bind_rows(chembl.links, drugbank.links, dgidb.links)
db.links$link <- gsub(">chembl<", ">ChEMBL<", db.links$link)
db.links$link <- gsub(">dgidb<", ">DGIdb<", db.links$link)
db.links$link <- gsub(">drugbank<", ">DrugBank<", db.links$link)

write.table(db.links, "Data/db_external_links.txt", sep = "\t", row.names = F)
syn$store(synapse$File("Data/db_external_links.txt", parentId = "syn12978846"), used = c("syn12978912","syn12972665","syn12684108"), executed = this.file)

### gene links
db <- readRDS(syn$get("syn12978910")$path)
genes <- unique(db$hugo_gene)
link <- sapply(genes, function(x){
  paste0("<a href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=",x,"', target = '_blank'>GeneCards</a>")
})

link <- as.data.frame(link) %>% rownames_to_column("hugo_gene")

write.table(link, "Data/gene_external_links.txt", sep = "\t", row.names = F)
syn$store(synapse$File("Data/gene_external_links.txt", parentId = "syn12978846"), used = c("syn12978910"), executed = this.file)



##drug response dataset preparation
ctrp.structures <- read.table(syn$get("syn5632193")$path, header = T, sep = "\t", quote = "", comment.char = "") 
ctrp.structures$makenames <- make.names(ctrp.structures$cpd_name)

drug.resp <- read.table(syn$get("syn7466611")$path, sep = "\t", header = TRUE) #%>% 
#rownames_to_column("cellLine") %>% 
#gather(makenames, auc, -cellLine)

ctrp.structures$cpd_smiles <- pbmclapply(unique(ctrp.structures$cpd_smiles), function(x){
  tryCatch({
    molvs$standardize_smiles(r_to_py(x))
  }, warning = function(w) {
    print(paste("structure",x,"has an issue"))
  }, error = function(e) {
    print(paste("structure",x,"is invalid"))
  })
}) %>% unlist()


fp.ctrp <- parseInputFingerprint(as.character(unique(ctrp.structures$cpd_smiles)), type = "circular")

saveRDS(drug.resp, "Data/drugresp.rds")
saveRDS(ctrp.structures, "Data/ctrpstructures.rds")

fp.ctrp.extended <- parseInputFingerprint(as.character(unique(ctrp.structures$cpd_smiles)), type = "extended")
fp.ctrp.circular <- parseInputFingerprint(as.character(unique(ctrp.structures$cpd_smiles)), type = "circular")
fp.ctrp.maccs <- parseInputFingerprint(as.character(unique(ctrp.structures$cpd_smiles)), type = "maccs")

saveRDS(fp.ctrp.extended, "Data/fpctrp_extended.rds")
saveRDS(fp.ctrp.circular, "Data/fpctrp_circular.rds")
saveRDS(fp.ctrp.maccs, "Data/fpctrp_maccs.rds")


####Sanger Data

sang.structures <- read.table("Data/sanger_structures.txt", header = T, sep = "\t", quote = "", comment.char = "") %>% filter(smiles != "none found")
sang.structures$makenames <- make.names(sang.structures$sanger_names)

drug.resp.sang <- read.table(syn$get("syn9987866")$path, sep = "\t", header = TRUE) %>% rownames_to_column("cellLine")
drug.resp.sang <- drug.resp.sang %>% gather(drug, auc, -cellLine)
drug.resp.sang$drug <- gsub("_.+$", "", drug.resp.sang$drug)
drug.resp.sang <- drug.resp.sang %>% 
  group_by(cellLine, drug) %>% 
  summarize(medianAUC = median(auc)) %>% 
  ungroup() %>% 
  spread(drug, medianAUC) %>% 
  remove_rownames() %>% 
  column_to_rownames("cellLine")

sang.structures$smiles <- pbmclapply(sang.structures$smiles, function(x){
  tryCatch({
    molvs$standardize_smiles(r_to_py(x))
  }, warning = function(w) {
    print(paste("structure",x,"has an issue"))
  }, error = function(e) {
    print(paste("structure",x,"is invalid"))
  })
}) %>% unlist()


fp.sang.extended <- parseInputFingerprint(as.character(unique(sang.structures$smiles)), type = "extended")
fp.sang.circular <- parseInputFingerprint(as.character(unique(sang.structures$smiles)), type = "circular")
fp.sang.maccs <- parseInputFingerprint(as.character(unique(sang.structures$smiles)), type = "maccs")

saveRDS(drug.resp.sang, "Data/drugresp_sang.rds")
saveRDS(sang.structures, "Data/sangstructures.rds")

saveRDS(fp.sang.extended, "Data/fpsang_extended.rds")
saveRDS(fp.sang.circular, "Data/fpsang_circular.rds")
saveRDS(fp.sang.maccs, "Data/fpsang_maccs.rds")

