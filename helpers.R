
db <- readRDS(synGet("syn11712148")$path) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(-n)

db.names <- readRDS(synGet("syn11712154")$path)
db$internal_id <- as.character(db$internal_id)

fp.extended <- readRDS(synGet("syn11889426")$path)
fp.extended <- fp.extended[names(fp.extended) %in% unique(db$internal_id)]

fp.circular <- readRDS(synGet("syn11808628")$path)
fp.circular <- fp.circular[names(fp.circular) %in% unique(db$internal_id)]

# fp.kr <- readRDS(synGet("syn11898431")$path)
# fp.kr <- fp.kr[names(fp.kr) %in% unique(db$internal_id)]

fp.maccs <- readRDS(synGet("syn11899073")$path)
fp.maccs <- fp.maccs[names(fp.maccs) %in% unique(db$internal_id)]

# fp.pubchem <- readRDS(synGet("syn12031332")$path)
# fp.pubchem <- fp.pubchem[names(fp.pubchem) %in% unique(db$internal_id)]


db.links <- read.table(synGet("syn11932224")$path, sep = "\t", header = T)
db.gene.links <- read.table(synGet("syn11941368")$path, sep = "\t", header = T)
db.genes <- db.gene.links$hugo_gene

##get CTRP data for heatmap

drug.resp <- readRDS("Data/drugresp.rds")
ctrp.structures <- readRDS("Data/ctrpstructures.rds")
fp.ctrp <- readRDS("Data/fpctrp.rds")

##same as above but for prepared sanger data

drug.resp.sang<-readRDS("Data/drugresp_sang.rds")
sang.structures<-readRDS("Data/sangstructures.rds")
fp.sang<-readRDS("Data/fpsang.rds")

