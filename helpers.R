
db <- readRDS("Data/drug_target_associations_v2.rds")


db.names <- readRDS("Data/distinct_compound_names.rds")
db$internal_id <- as.character(db$internal_id)

fp.extended <- readRDS("Data/db_fingerprints_extended.rds")

fp.circular <- readRDS("Data/db_fingerprints_circular.rds")

# fp.kr <- readRDS(synGet("syn11898431")$path)
# fp.kr <- fp.kr[names(fp.kr) %in% unique(db$internal_id)]

fp.maccs <- readRDS("Data/db_fingerprints_maccs.rds")

# fp.pubchem <- readRDS(synGet("syn12031332")$path)
# fp.pubchem <- fp.pubchem[names(fp.pubchem) %in% unique(db$internal_id)]


db.links <- read.table("Data/db_external_links.txt", sep = "\t", header = T)
db.gene.links <- read.table("Data/gene_external_links.txt", sep = "\t", header = T)
db.genes <- db.gene.links$hugo_gene

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

