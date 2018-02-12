
db <- readRDS(synGet("syn11712148")$path) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(internal_id, common_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative, known_selectivity_index, confidence)

mini.db <- db %>% group_by(internal_id) %>% 
  mutate(n = sum(n_qualitative, n_quantitative, na.rm = T)) %>% 
  select(internal_id, common_name, n) %>% 
  distinct() %>% 
  filter(n>=10)

db.names <- readRDS(synGet("syn11712154")$path)

db$internal_id <- as.character(db$internal_id)

fp.db <- readRDS(synGet("syn11693143")$path)
fp.db <- fp.db[names(fp.db) %in% unique(db$internal_id)]

fp.snappy <- fp.db[names(fp.db) %in% unique(mini.db$internal_id)]

db.genes <- unique(db$hugo_gene)

##get CTRP data for heatmap

drug.resp <- readRDS("Data/drugresp.rds")
ctrp.structures <- readRDS("Data/ctrpstructures.rds")
fp.ctrp <- readRDS("Data/fpctrp.rds")


##same as above but for prepared sanger data

drug.resp.sang<-readRDS("Data/drugresp_sang.rds")
sang.structures<-readRDS("Data/sangstructures.rds")
fp.sang<-readRDS("Data/fpsang.rds")

