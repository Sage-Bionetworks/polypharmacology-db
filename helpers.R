

db <- readRDS("Data/drug_target_associations_v2.rds") %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(-n)


db.names <- readRDS("Data/compound_names.rds")
db$internal_id <- as.character(db$internal_id)

fp.extended <- readRDS("Data/db_fingerprints_extended.rds")
fp.extended <- fp.extended[names(fp.extended) %in% unique(db$internal_id)]

fp.circular <- readRDS("Data/db_fingerprints_circular.rds")
fp.circular <- fp.circular[names(fp.circular) %in% unique(db$internal_id)]

# fp.kr <- readRDS(synGet("syn11898431")$path)
# fp.kr <- fp.kr[names(fp.kr) %in% unique(db$internal_id)]

fp.maccs <- readRDS("Data/db_fingerprints_maccs.rds")
fp.maccs <- fp.maccs[names(fp.maccs) %in% unique(db$internal_id)]

# fp.pubchem <- readRDS(synGet("syn12031332")$path)
# fp.pubchem <- fp.pubchem[names(fp.pubchem) %in% unique(db$internal_id)]


db.links <- read.table("Data/db_external_links.txt", sep = "\t", header = T)
db.gene.links <- read.table("Data/gene_external_links.txt", sep = "\t", header = T)
db.genes <- db.gene.links$hugo_gene

##get CTRP data for heatmap

drug.resp <- readRDS("Data/drugresp.rds")
ctrp.structures <- readRDS("Data/ctrpstructures.rds")
fp.ctrp <- readRDS("Data/fpctrp.rds")

##same as above but for prepared sanger data

drug.resp.sang<-readRDS("Data/drugresp_sang.rds")
sang.structures<-readRDS("Data/sangstructures.rds")
fp.sang<-readRDS("Data/fpsang.rds")

