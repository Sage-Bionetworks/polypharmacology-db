library(tidyverse)
library(synapser)
synLogin()

db <- readRDS(synGet("syn11712148")$path) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(internal_id, common_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative)

hist(db$mean_pchembl)
hist(log(db$n_quantitative))
hist(log(db$n_qualitative))

sum(db$n_quantitative, na.rm = T)+sum(db$n_qualitative, na.rm = T)
length(unique(db$internal_id))
length(unique(db$hugo_gene))

CHEMBL3356433 <- fp.db[[130105]]
compound30 <- fp.db[[59664]]
limk.sim <- distance(CHEMBL3356433, compound30, method = "tanimoto")



