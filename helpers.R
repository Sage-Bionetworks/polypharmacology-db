
db <- readRDS("Data/compound_target_associations_v4.rds")
db.genes <- db$hugo_gene %>% unique()

db_structures <- readRDS("Data/compound_structures_v4.rds")

db.names <- readRDS("Data/distinct_compound_names_v4.rds")

db$inchikey <- as.character(db$inchikey)

db <- select(db.names, inchikey, pref_name) %>% 
  distinct() %>% 
  full_join(db)
  
fp.extended <- readRDS("Data/db_fingerprints_extended_v4.rds")

fp.circular <- readRDS("Data/db_fingerprints_circular_v4.rds")

fp.maccs <- readRDS("Data/db_fingerprints_maccs_v4.rds")

db.links <- readRDS("Data/cmpd_links_v4.rds")
db.gene.links <- readRDS("Data/gene_links_v4.rds")

##get CTRP data for heatmap

drug.resp <- readRDS("Data/drugresp.rds")
ctrp.structures <- readRDS("Data/ctrpstructures.rds")
fp.ctrp.extended <- readRDS("Data/fpctrp_extended.rds")
fp.ctrp.circular <- readRDS("Data/fpctrp_circular.rds")
fp.ctrp.maccs <- readRDS("Data/fpctrp_maccs.rds")

##same as above but for prepared sanger data

drug.resp.sang<-readRDS("Data/drugresp_sang.rds")
sang.structures<-readRDS("Data/sangstructures.rds")
fp.sang.extended<-readRDS("Data/fpsang_extended.rds")
fp.sang.circular<-readRDS("Data/fpsang_circular.rds")
fp.sang.maccs<-readRDS("Data/fpsang_maccs.rds")

