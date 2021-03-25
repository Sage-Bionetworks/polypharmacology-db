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
library(synapser)
library(pbmcapply)
synLogin()

this.file <- "https://raw.githubusercontent.com/Sage-Bionetworks/polypharmacology-db/master/build_db_v3_0_1_map_v2_targets.R"


used.these <- c("syn17091006","syn12913659","syn12684108",
                "syn12683504","syn11685587","syn12973247",
                "syn11695142","syn17090819")


v2_to_v3 <- synTableQuery('select * from syn17091006')$filepath %>% 
  read.csv() 

##Targets

##DrugBank 5.0.11
drugbank.targets <- read.table(synGet("syn12913659")$path, strip.white = TRUE) %>% 
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

dgidb.targets <- read.table(synGet("syn12684108")$path, sep = "\t", quote = "",header = T, strip.white = TRUE) %>%
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

cp.targets <- read.csv(synGet("syn12683504")$path, header = T, strip.white = TRUE) %>% 
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
  left_join(x=v2_to_v3, y = .) %>% 
  dplyr::select(internal_id, hugo_gene, n) %>% 
  group_by(internal_id, hugo_gene) %>% 
  summarize("n_qualitative" = sum(n)) %>% 
  ungroup() %>% 
  filter(!is.na(hugo_gene))

######### Summarize Quantitative Data
##Klaeger 2017
klaeger.drugnames <- t(read.table(synGet("syn11685587")$path, sep = "\t",
                                  header = F, comment.char = "", quote = "")[1,4:246]) %>% ##removing metadata names
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
  left_join(v2_to_v3) 

klaeger_pchembl <- klaeger.targets %>% 
  select(internal_id, hugo_gene, pchembl_value)

klaeger_assaytype_summary <- klaeger.targets %>% 
  select(internal_id, hugo_gene, Kd) %>% 
  set_names("internal_id", "hugo_gene", "standard_value") %>% 
  mutate(standard_type = "Kd")

##ChEMBL v23 - SQL query on synapse file page V
chembl.targets <- read.table(synGet("syn12973247")$path, sep = "\t",
                             header = T, comment.char = "", quote = "\"", strip.white = TRUE)

##filter out non-hugo genes introduced by chembl synonyms
hugo <- read.table(synGet("syn11695142")$path, header = T, sep = "\t", quote = "", comment.char = "", colClasses = "character") 

chembl.targets <- dplyr::filter(chembl.targets, component_synonym %in% as.character(hugo$symbol))%>% 
  filter(pchembl_value != "NULL") #remove null pchembl values (small percentage)

pchembl.summary <- chembl.targets %>% 
  mutate(external_id = as.character(molregno)) %>% 
  left_join(x=v2_to_v3, y = .) %>% 
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
  left_join(x=v2_to_v3, y = .) %>% 
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

##add names back

v3_names <- synTableQuery("select internal_id, std_name from syn17090819")$filepath %>% 
  read.csv() %>% 
  select(-ROW_ID, -ROW_VERSION)

full.db <- left_join(full.db, v3_names)
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

write.table(full.db2, "Data/drug_target_associations_v3_0_1.txt", row.names = F)
synStore(File("Data/drug_target_associations_v3_0_1.txt", parentId = "syn17091504"), executed = this.file, 
          used = used.these)

saveRDS(full.db2, "Data/drug_target_associations_v3_0_1.rds")
synStore(File("Data/drug_target_associations_v3_0_1.rds", parentId = "syn17091504"), executed = this.file, 
          used = used.these)


####Generate fingerprints for database.

full.db2 <- readRDS(syn$get("syn17091507")$path)

structures <- syn$tableQuery("SELECT * FROM syn17090819")$filepath %>% 
  read.csv() %>% 
  select(internal_id, std_smiles) %>% 
  distinct()

valid <- as.character(structures$std_smiles)

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

names(foo) <- structures$internal_id

saveRDS(foo, "Data/db_fingerprints_extended_v3_0_1.rds")
syn$store(synapse$File("Data/db_fingerprints_extended_v3_0_1.rds", parentId = "syn17091504"), used = c("syn17091507","syn17090819"), executed = this.file)

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

names(foo) <- structures$internal_id

saveRDS(foo, "Data/db_fingerprints_circular_v3_0_1.rds")
syn$store(synapse$File("Data/db_fingerprints_circular_v3_0_1.rds", parentId = "syn17091504"), used = c("syn17091507","syn17090819"), executed = this.file)

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

names(foo) <- structures$internal_id

saveRDS(foo, "Data/db_fingerprints_maccs_v3_0_1.rds")
syn$store(synapse$File("Data/db_fingerprints_maccs_v3_0_1.rds", parentId = "syn17091504"), used = c("syn17091507","syn17090819"), executed = this.file)

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
db.names <- syn$tableQuery("SELECT * FROM syn17091006")$filepath %>% 
  read.csv()

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
  mutate(link = sapply(as.character(external_id), function(x) URLencode(x))) %>%
  mutate(link = paste0("<a href='https://www.drugbank.ca/drugs/",link,"', target = '_blank'>",database,"</a>")) %>% 
  select(internal_id, link, external_id, database) %>% 
  distinct()

db.links <- bind_rows(chembl.links, drugbank.links, dgidb.links)
db.links$link <- gsub(">chembl<", ">ChEMBL<", db.links$link)
db.links$link <- gsub(">dgidb<", ">DGIdb<", db.links$link)
db.links$link <- gsub(">drugbank<", ">DrugBank<", db.links$link)

write.table(db.links, "Data/db_external_links.txt", sep = "\t", row.names = F)
syn$store(synapse$File("Data/db_external_links.txt", parentId = "syn17091504"), used = c("syn12978912","syn12972665","syn12684108"), executed = this.file)

### gene links
db <- readRDS(syn$get("syn17091507")$path)
genes <- unique(db$hugo_gene)
link <- sapply(genes, function(x){
  paste0("<a href = 'http://www.genecards.org/cgi-bin/carddisp.pl?gene=",x,"', target = '_blank'>GeneCards</a>")
})

link <- as.data.frame(link) %>% rownames_to_column("hugo_gene")

write.table(link, "Data/gene_external_links.txt", sep = "\t", row.names = F)
syn$store(synapse$File("Data/gene_external_links.txt", parentId = "syn17091504"), used = c("syn12978910"), executed = this.file)



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

