library(ChemmineR)
library(dplyr)
library(rJava)
library(rcdk)
library(fingerprint)
library(parallel)

##import drugbank 5.0 structure data from Drugbank 5.0
##to read in sdf file in sensible way
sdfset<-read.SDFset("NoGit/structures 2.sdf")
valid <- validSDF(sdfset)
sdfset<-sdfset[valid]

struct <- datablock(sdfset) 

x<-as.data.frame(t(sapply(names(struct), function(i){
  p<-unlist(struct[i], use.names = F)
  p[1:3]
})))

colnames(x) <- c("external_id","x", "smiles")

db_struct <- x %>% dplyr::select(external_id, smiles)

##import CHeMBL v23 structure data downloaded from ftp site
chembl_struct <- read.csv("NoGit/chembl_structures.csv", comment.char = "") %>% distinct()
colnames(chembl_struct) <- c("external_id", "smiles")
chembl_struct$external_id <- as.factor(chembl_struct$external_id)


##import DGIDB
dgidb_struct <- read.table("Data/dgidb_smiles_from_pubchem.txt", 
                           sep = "\t",
                           quote = "",
                           header = F, 
                           comment.char = "",
                           na.strings = NA)
dgidb_struct$count <- c(1:nrow(dgidb_struct))
dgidb_struct <- dgidb_struct %>% 
  group_by(V1) %>% 
  top_n(1, -count) %>% 
  select(V1,V2) %>% 
  filter(V2 != "") %>% 
  ungroup()

colnames(dgidb_struct) <- c("external_id", "smiles")

structures <- bind_rows(chembl_struct, db_struct, dgidb_struct)

##save fingerprints as RData, calculate fp.sim.matrix on beefy ec2 instance with many threads

fps <- parseInputFingerprint(structures$smiles[1:100])
fpnames <- as.data.frame(names(fps))
fpnames$flag <- rep(T, nrow(fpnames))
fpnames$same <- rep(NA, nrow(fpnames))

list <- list()

for(i in fpnames$`names(fps)`[fpnames$flag==T]) {
  sim <- mclapply(fpnames$`names(fps)`[fpnames$flag==T & fpnames$`names(fps)`!= i], function(j) {
    distance(fps[[i]], fps[[j]])
  }, mc.cores = detectCores())
  names(sim) <- fpnames$`names(fps)`[fpnames$flag==T]
  sim <- sim[sim==1]
  list[[i]] <- sim
  fpnames$flag[fpnames$`names(fps)` %in% names(sim)] <- F
}
