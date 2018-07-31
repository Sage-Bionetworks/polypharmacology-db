library(tidyverse)
library(synapser)
synLogin()

db <- readRDS(synGet("syn12978910", version = 8)$path) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(internal_id, common_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative)

sum(db$n_quantitative, na.rm = T)+sum(db$n_qualitative, na.rm = T)
length(unique(db$internal_id))
length(unique(db$hugo_gene))

CHEMBL3356433 <- fp.db[[130105]]
compound30 <- fp.db[[59664]]
limk.sim <- distance(CHEMBL3356433, compound30, method = "tanimoto")

###
#'Response to reviewers
#'Summary of DTE

dte.metrics <-
  tibble(
    "metric" = c("molecules",
                 "targets",
                 "associations",
                 "quantitative interactions",
                 "qualitative interactions"),
    "number" = c(length(unique(db$internal_id)),
                 length(unique(db$hugo_gene)),
                 nrow(db),
                 sum(db$n_quantitative, na.rm = T),
                 sum(db$n_qualitative,na.rm=T))
  )

###
#'Response to reviewers
#'Summary of source databases

grouped_structures <- read.table(synGet("syn12978848")$path, 
                                 sep = "\t", 
                                 comment.char = "",
                                 header = T,
                                 colClasses = c("character", 
                                                "character",
                                                "character",
                                                "character")) %>% 
  # select(2:6) %>%
  filter(smiles != "") %>% 
  filter(internal_id %in% db$internal_id)

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
chembl.targets <- read.table(synGet("syn12973247")$path, sep = "\t",
                             header = T, comment.char = "", quote = "\"", strip.white = TRUE)

##filter out non-hugo genes introduced by chembl synonyms
hugo <- read.table(synGet("syn11695142")$path, header = T, sep = "\t", quote = "", comment.char = "", colClasses = "character") 

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


klaeger.filt <- filter(klaeger.targets, external_id %in% grouped_structures$external_id)
chembl.filt <- filter(chembl.targets, molregno %in% grouped_structures$external_id)
dgidb.filt <- filter(dgidb.targets, external_id %in% grouped_structures$external_id)
drugbank.filt <- filter(drugbank.targets, external_id %in% grouped_structures$external_id)
cp.filt <- filter(cp.targets, external_id %in% grouped_structures$external_id)

dte.source.metrics <-
  tibble(
    "metric" = c("molecules",
                 "targets",
                 "associations"),
    "ChEMBL v 24.1" = c(length(unique(chembl.filt$molregno)),
                  length(unique(chembl.filt$component_synonym)),
                  nrow(chembl.filt)), 
    "ChemicalProbes" = c(length(unique(cp.filt$external_id)),
                         length(unique(cp.filt$hugo_gene)),
                         nrow(cp.filt)), 
    "DGIdb" = c(length(unique(dgidb.filt$external_id)),
                length(unique(dgidb.filt$hugo_gene)),
                nrow(dgidb.filt)),
    "DrugBank" = c(length(unique(drugbank.filt$external_id)),
                   length(unique(drugbank.filt$hugo_gene)),
                   nrow(drugbank.filt)),
    "Klaeger et al" = c(length(unique(klaeger.filt$external_id)),
                        length(unique(klaeger.filt$hugo_gene)),
                        nrow(klaeger.filt))
  )

write.csv(dte.source.metrics, "DTE_source_metrics.csv", row.names = F)


