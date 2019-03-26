library(dplyr)
library(tidyr)
library(tibble)
library(reticulate)
library(pbapply)
use_python("/usr/local/bin/python2")
molvs <- import("molvs")
synapse <- import("synapseclient")
syn <- synapse$Synapse()
syn$login()

##NOTE - smiles column in this script is the std_smiles from V2
##NOTE - std_smiles column in this script is the standardized output smiles from the V3 python script,
##but the standardization approach should be identical to V2, so these columns are expected to be identical!!!!

##syn.tableQuery is giving me an issue with reticulate, so just manually download and read in syn17008823
db <-read.csv("Data/v2_table.csv", comment.char = "", quote = "\"", allowEscapes = T)

##I'll be ignoring the manually withdrawn structure from V2 (i.e. not withdrawing them), because many of these should be taken care of 
##by switch to inchi. I also don't want to reassign original_smiles values for backwards compatibility

##import python-standardized smiles 
# query <- syn$get("syn17011195")$path

std_smiles<-read.csv("Data/standardized_structures.csv", comment.char = "", quote = "\"", allowEscapes = F, header = T)

database<- db %>% ##update database with newly processed structures
  select(internal_id, std_name, smiles, common_name) %>% 
  distinct() %>% 
  left_join(std_smiles) %>%  ##add standardized smiles to database
  group_by(inchi) %>% ##group by these standardized smiles
  mutate(min_id = min(internal_id)) %>% ##set aside the minimum id for each group
  mutate(internal_id = min_id) %>% ##assign all ids to min id within group
  ungroup() %>% 
  select(-min_id) %>% 
  distinct() ##eliminate identical rows caused by combining

names <- database %>% 
  group_by(internal_id) %>% ##group by these standardized smiles
  slice(1) %>% 
  ungroup() %>% 
  select(internal_id, std_name)

database <- database %>%
  select(-std_name) %>% 
  left_join(names) %>% 
  distinct()

withdrawn_1<-db %>% ##build table of withdrawn structures to separate out into own 
  select(internal_id, std_name, smiles) %>% 
  distinct() %>% 
  left_join(std_smiles) %>% 
  group_by(inchi) %>%
  mutate(min_id = min(internal_id)) %>% ##all of these steps are same as previous piped command
  mutate(withdrawn = case_when(internal_id==min_id ~ FALSE, ## in cases where internal_id is the min_id, withdrawn == FALSE
                               internal_id!=min_id ~ TRUE)) %>%  ## in cases where internal_id is NOT the min_id, withdrawn == TRUE
  mutate(withdrawn_to = case_when(internal_id!=min_id ~ min_id)) %>% ##in cases where internal_id is not the min_id, assign withdrawn_to to min_id to track
  ungroup() %>% 
  filter(withdrawn == TRUE) %>% ## get rid of non withdrawn
  select(-min_id) %>% 
  distinct()

###there are some (102) std_names that are duplicated or more, so presumably there are duplicate compounds

database_count <- database %>%  ##this will quickly evaluate replicated std_names
  select(internal_id, std_name) %>% 
  distinct() %>% 
  add_count(std_name) %>% 
  filter(n >1)

##so there are 102 names that are associated with more than one internal_id. Let's group and 
##withdraw the duplicates here, too. 
database_dupes <- database %>% filter(internal_id %in% database_count$internal_id) %>% 
  group_by(std_name) %>% ##group by standardized names
  mutate(internal_id = min(internal_id)) %>%
  mutate(inchi = nth(inchi, 1)) %>%
  mutate(std_smiles = nth(std_smiles, 1)) %>%
  mutate(inchikey = nth(inchikey, 1)) %>%
  ungroup() %>% 
  distinct() ##eliminate identical rows caused by combining
  
withdrawn_2 <- database %>% filter(internal_id %in% database_count$internal_id) %>% 
  group_by(std_name) %>% ##group by these standardized smiles
  mutate(min_id = min(internal_id)) %>%
  mutate(withdrawn = case_when(internal_id==min_id ~ FALSE, ## in cases where internal_id is the min_id, withdrawn == FALSE
                               internal_id!=min_id ~ TRUE)) %>%  ## in cases where internal_id is NOT the min_id, withdrawn == TRUE
  mutate(withdrawn_to = case_when(internal_id!=min_id ~ min_id)) %>% 
  ungroup() %>%
  filter(withdrawn == TRUE) %>% 
  select(-min_id, -common_name) %>% 
  distinct() #eliminate identical rows caused by combining
    
  
all_withdrawn <- bind_rows(withdrawn_1, withdrawn_2) %>% 
  mutate_all(as.character)

final_db <- database %>% 
  filter(!internal_id %in% database_count$internal_id) %>% ##remove all dupes from db 
  bind_rows(database_dupes)  ##add back df with combined duplicates

final_db_fix_smiles <- final_db %>%  ##need to extract std_smiles and harmonize to inchi groups
  ##because there are still in some cases more than one smiles per inchikey
  select(std_smiles, inchi) %>% 
  group_by(inchi) %>% 
  slice(1) %>% ##just pick first smiles per group to eliminate variant std_smiles 
  distinct() %>% 
  ungroup()

final_db <- final_db %>% ##remove smiles column and replace with uniformized std_smiles :)
  select(-std_smiles) %>% 
  left_join(final_db_fix_smiles) 

final_db_structures <- final_db %>% 
  select(internal_id, std_name, inchikey, inchi, std_smiles) %>% 
  distinct() %>% 
  mutate_all(as.character) %>% 
  mutate(internal_id = as.numeric(internal_id)) %>% 
  arrange(internal_id)

#######Structures Table Unit Tests
if(nrow(table(final_db_structures$internal_id) %>% as.data.frame() %>% filter(Freq>1)) > 0 |
   nrow(table(final_db_structures$std_name) %>% as.data.frame() %>% filter(Freq>1)) > 0 |
   nrow(table(final_db_structures$inchikey) %>% as.data.frame() %>% filter(Freq>1)) > 0 |
   nrow(table(final_db_structures$inchi) %>% as.data.frame() %>% filter(Freq>1)) > 0 |
   nrow(table(final_db_structures$std_smiles) %>% as.data.frame() %>% filter(Freq>1)) > 0){
  print("structures table invalid")
}else{
  print("structures table valid")
}

final_db_synonyms <- final_db %>% 
  select(internal_id, common_name) %>% 
  distinct() %>% 
  mutate_all(as.character) %>% 
  mutate(internal_id = as.numeric(internal_id)) %>% 
  arrange(internal_id) %>% 
  distinct()

if(nrow(table(final_db_synonyms$internal_id) %>% as.data.frame() %>% filter(is.na(Var1)))|
   nrow(table(final_db_synonyms$common_name) %>% as.data.frame() %>% filter(is.na(Var1)))){
  print("synonyms table invalid")
}else{
  print("synonyms table valid")
}

grouped_structures <- syn$get("syn12978848")$path %>% 
  read.table(., header = T, sep  = "\t", quote = "\"", comment.char = "") %>% 
  select(internal_id, original_smiles, smiles, external_id, database) %>% 
  purrr::set_names(c("internal_id_v2", "original_smiles", "v2_smiles", "external_id", "database"))

final_db_previous_structures <- final_db %>% 
  select(internal_id, smiles) %>% 
  distinct() %>% 
  purrr::set_names(c("internal_id", "v2_smiles")) %>% 
  full_join(grouped_structures) %>% 
  filter(!is.na(internal_id)) %>% 
  distinct() %>% 
  mutate_all(as.character) %>% 
  mutate(internal_id = as.numeric(internal_id))

project = syn$get('syn11672851')
    
cols = list(
  synapse$Column(name='internal_id', columnType='INTEGER'),
  synapse$Column(name='std_name', columnType='LARGETEXT'),
  synapse$Column(name='inchikey', columnType='STRING'),
  synapse$Column(name='inchi', columnType='LARGETEXT'),
  synapse$Column(name='std_smiles', columnType='LARGETEXT'))

schema = synapse$Schema(name='v3 IDs, Structures, Standard Name', columns=cols, parent=project)
table = synapse$Table(schema, final_db_structures)
syn$store(table)

cols = list(
  synapse$Column(name='internal_id', columnType='INTEGER'),
  synapse$Column(name='common_name', columnType='LARGETEXT')
)
schema = synapse$Schema(name='v3 IDs, Synonyms', columns=cols, parent=project)
table = synapse$Table(schema, final_db_synonyms)
syn$store(table)

cols = list(
  synapse$Column(name='internal_id', columnType='INTEGER'),
  synapse$Column(name='v2_smiles', columnType='LARGETEXT'),
  synapse$Column(name='internal_id_v2', columnType='STRING'),
  synapse$Column(name='original_smiles', columnType='LARGETEXT'),
  synapse$Column(name='external_id', columnType='LARGETEXT'),
  synapse$Column(name='database', columnType='STRING'))

schema = synapse$Schema(name='v3 IDs, v2 IDs, and Source Structures', columns=cols, parent=project)
table = synapse$Table(schema, final_db_previous_structures)
syn$store(table)

cols = list(
  synapse$Column(name='internal_id', columnType='INTEGER'),
  synapse$Column(name='std_name', columnType='LARGETEXT'),
  synapse$Column(name='smiles', columnType='LARGETEXT'),
  synapse$Column(name='std_smiles', columnType='LARGETEXT'),
  synapse$Column(name='inchi', columnType='LARGETEXT'),
  synapse$Column(name='inchikey', columnType='STRING'),
  synapse$Column(name='common_name', columnType='LARGETEXT'),
  synapse$Column(name='withdrawn', columnType='STRING'),
  synapse$Column(name='withdrawn_to', columnType='STRING'))

schema = synapse$Schema(name='v3 Withdrawn and Deprecated Rows', columns=cols, parent=project)
table = synapse$Table(schema, all_withdrawn)
syn$store(table)

