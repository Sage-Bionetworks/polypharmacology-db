
db <- readRDS(synGet("syn11712148")$path) %>% 
  filter(!is.na(hugo_gene)) 

db.names <- readRDS(synGet("syn11712154")$path)
db$internal_id <- as.character(db$internal_id)

fp.db <- readRDS(synGet("syn11683261")$path)
fp.db <- fp.db[names(fp.db) %in% unique(db$internal_id)]

db.genes <- unique(db$hugo_gene)

##get CTRP data for heatmap

drug.resp <- readRDS("Data/drugresp.rds")
ctrp.structures <- readRDS("Data/ctrpstructures.rds")
fp.ctrp <- readRDS("Data/fpctrp.rds")

##same as above but for prepared sanger data

drug.resp.sang<-readRDS("Data/drugresp_sang.rds")
sang.structures<-readRDS("Data/sangstructures.rds")
fp.sang<-readRDS("Data/fpsang.rds")

