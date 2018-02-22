library()

db <- readRDS(synGet("syn11712148")$path) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(internal_id, common_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative)


CHEMBL3356433 <- fp.db[[130105]]
compound30 <- fp.db[[59664]]
limk.sim <- distance(CHEMBL3356433, compound30, method = "tanimoto")
