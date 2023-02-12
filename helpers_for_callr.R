db_structures <- readRDS("Data/compound_structures_v4.rds")

db.names <- readRDS("Data/distinct_compound_synonyms_v4.rds")

fp.extended <- readRDS("Data/db_fingerprints_extended_v4.rds")

fp.circular <- readRDS("Data/db_fingerprints_circular_v4.rds")

fp.maccs <- readRDS("Data/db_fingerprints_maccs_v4.rds")


#----- repeated from helpers.R ---------
db <- readRDS("Data/compound_target_associations_v4.rds")

db.names <- readRDS("Data/distinct_compound_synonyms_v4.rds")

db$inchikey <- as.character(db$inchikey)

db <- select(db.names, inchikey, pref_name) %>% 
  distinct() %>% 
  full_join(db)

#---------------------------------------------------
