library(tidyverse)
library(synapser)
synLogin()

v1.structures <- read.table(synGet("syn11685572")$path, sep ="\t", header = T, comment.char = "") %>% 
  distinct() %>% 
  group_by(external.id) %>% 
  top_n(1) %>% 
  ungroup

v2.names <- read.table(synGet("syn12683504")$path, sep = ",", header = T) %>% 
  select(Probe.Name) %>% 
  set_names("external.id") %>% 
  full_join(v1.structures)

write.csv(v2.names, "Data/chemicalprobes_structures_07_02_2018.csv")
synStore(File("Data/chemicalprobes_structures_07_02_2018.csv", parentId = "syn12683460"))
