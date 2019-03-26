
### create igraph object from db
##for fendr

this.file <- 
library(synapser)
library(igraph)
db <- readRDS(synGet("syn17091507")$path) %>%
  filter(!is.na(hugo_gene)) %>%
  select(internal_id, hugo_gene, mean_pchembl, n_quantitative, n_qualitative) %>%
  group_by(internal_id, hugo_gene) %>%
  mutate(total_n = sum(n_quantitative, n_qualitative, na.rm = T)) %>%
  ungroup()

db.names <- readRDS(synGet("syn17091507")$path) %>%
  filter(!is.na(hugo_gene)) %>%
  select(internal_id, std_name) %>%
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

names <- readRDS(synTableQuery("SELECT internal_id, std_name FROM syn17090819")$asDataFrame())
write.table(names, "Data/drug-target_explorer_igraph_name_map.txt", sep = "\t", row.names = F)

results <- synTableQuery(sprintf("select * from %s", "syn11831632"))
x <- nrow(results$asDataFrame())/10000
for(i in 1:ceiling(x)){
  print(i)
  results <- synTableQuery(sprintf("select * from %s limit 10000", "syn11831632"))
  deleted <- synDelete(results$asRowSet())
}

results <- synGet("syn11831632")
tableToAppend <- Table(results, names)
table <- synStore(tableToAppend)