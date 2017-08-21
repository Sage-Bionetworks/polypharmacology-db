library(shiny)
library(DT)
library(plyr)
library(dplyr)
library(rJava)
library(rcdk)
library(fingerprint)
library(pheatmap)
library(enrichR)
library(webchem)

evo <- readRDS("Data/evotec_dgidb.RDS")
evo$Structure_ID <- as.character(evo$Structure_ID)
evo$Common_Name <- as.character(evo$Common_Name)
evo <- evo %>% filter(N_quantitative >= N_inactive | N_qualitative >= N_inactive | N_DGIDB > 0)

fp.evo <- readRDS("Data/fpevo.rds")[unique(evo$Original_molecule_SMILES)]
syns <- readRDS("Data/commname.rds") %>% filter(Original_molecule_SMILES %in% unique(evo$Original_molecule_SMILES)) %>% 
  filter(!is.na(Common_Name))


parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
}

convertDrugToSmiles <- function(input) {
  filt <- filter(syns, Common_Name == input & Original_molecule_SMILES != "") %>% dplyr::select(Original_molecule_SMILES)
  filt
}

getTargetList <- function(selectdrugs) {
  targets <- filter(evo, Common_Name %in% selectdrugs) %>% dplyr::select(Common_Name, Hugo_Gene, MedianActivity_nM, N_quantitative, N_qualitative, 
                                                                          N_inactive, N_DGIDB) %>% arrange(-N_quantitative)
  if (nrow(targets) > 1) {
    targets
  } else {
    print("none found")
  }
}

getSimMols <- function(input, sim.thres) {
  input <- input
  fp.inp <- parseInputFingerprint(input)
  
  sims <- lapply(fp.inp, function(i) {
    sim <- sapply(fp.evo, function(j) {
      distance(i, j)
    })
    bar <- as.data.frame(sim)
    bar$match <- rownames(bar)
    bar
  })
  sims <- ldply(sims)
  sims2 <- sims %>% filter(sim >= sim.thres) %>% arrange(-sim)
  sims2$Original_molecule_SMILES <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$sim, 3)
  targets <- left_join(sims2, evo) %>% dplyr::select(Common_Name, `Tanimoto Similarity`) %>% distinct()
}

getNetwork <- function(input, sim.thres, selectdrugs) {
  input <- input
  
  fp.inp <- parseInputFingerprint(input)
  
  sims <- lapply(fp.inp, function(i) {
    sim <- sapply(fp.evo, function(j) {
      distance(i, j)
    })
    bar <- as.data.frame(sim)
    bar$match <- rownames(bar)
    bar
  })
  sims <- ldply(sims)
  sims2 <- sims %>% filter(sim >= sim.thres) %>% arrange(-sim)
  sims2$Original_molecule_SMILES <- as.character(sims2$match)
  targets <- left_join(sims2, evo) %>% dplyr::select(Common_Name, sim) %>% 
    distinct() %>% filter(Common_Name %in% selectdrugs)
  targets$from <- "input"
  targets$to <- as.character(targets$Common_Name)
  targets$width <- ((targets$sim)^2) * 10
  targets$color <- "red"
  targets <- select(targets, from, to, width, color)
}

getTargetNetwork <- function(input, sim.thres, selectdrugs) {
  input <- input
  fp.inp <- parseInputFingerprint(input)
  
  sims <- lapply(fp.inp, function(i) {
    sim <- sapply(fp.evo, function(j) {
      distance(i, j)
    })
    bar <- as.data.frame(sim)
    bar$match <- rownames(bar)
    bar
  })
  sims <- ldply(sims)
  sims2 <- sims %>% filter(sim >= sim.thres) %>% arrange(-sim)
  sims2$Original_molecule_SMILES <- as.character(sims2$match)
  targets <- left_join(sims2, evo) %>% distinct() %>% filter(Common_Name %in% selectdrugs)
  targets$from <- targets$Common_Name
  targets$to <- as.character(targets$Hugo_Gene)
  targets$width <- 5
  targets$color <- "green"
  targets <- select(targets, from, to, width, color) %>% 
    filter(from !="NA" & to != "NA")
}

dbs <- c("GO_Molecular_Function_2017", "GO_Cellular_Component_2017", "GO_Biological_Process_2017", 
         "KEGG_2016")

getGeneOntologyfromTargets <- function(selectdrugs) {
  selectdrugs <- selectdrugs
  targets <- getTargetList(selectdrugs)
  
  if (nrow(targets) > 1) {
    enriched <- enrichr(as.vector(targets$Hugo_Gene), dbs)
  } else {
    print("no targets")
  }
  
}


getMolsFromGenes <- function(inp.gene) {
  mols <- filter(evo, Hugo_Gene == inp.gene) %>% select(Structure_ID, 
                                                        Supplier_Molname, Common_Name, MedianActivity_nM, N_quantitative, N_qualitative, N_inactive, Original_molecule_SMILES) %>% 
    distinct()
  mols
}

getMolsFromGeneNetworks.edges <- function(inp.gene) {
  mols <- filter(evo, Hugo_Gene == inp.gene & N_quantitative > 1) %>% 
    select(Structure_ID, N_quantitative) %>% distinct()
  
  net <- filter(evo, Structure_ID %in% mols$Structure_ID & N_quantitative > 
                  1) %>% group_by(Structure_ID) %>% top_n(10, N_quantitative) %>% 
    ungroup()
  
  net$from <- as.character(net$Structure_ID)
  net$to <- as.character(net$Hugo_Gene)
  net$width <- net$N_quantitative/10
  net$color <- "red"
  net <- net %>% select(from, to, width, color)
  as.data.frame(net)
}

getMolsFromGeneNetworks.nodes <- function(inp.gene) {
  mols <- filter(evo, Hugo_Gene == inp.gene & N_quantitative > 1) %>% 
    select(Structure_ID, N_quantitative) %>% distinct()
  
  net <- filter(evo, Structure_ID %in% mols$Structure_ID & N_quantitative > 
                  1) %>% group_by(Structure_ID) %>% top_n(10, N_quantitative) %>% 
    ungroup()
  
  id <- c(unique(as.character(net$Structure_ID)), unique(as.character(net$Hugo_Gene)))
  label <- c(unique(as.character(net$Structure_ID)), unique(as.character(net$Hugo_Gene)))
  color <- c(rep("blue", length(unique(as.character(net$Structure_ID)))), 
             rep("green", length(unique(as.character(net$Hugo_Gene)))))
  
  net <- as.data.frame(cbind(id, label, color))
  
}

getSmiles <- function(input.name) {
  input.name <- input.name
  input.name <- URLencode(input.name)
  query <- as.vector(cir_query(input.name, representation = "smiles", first = TRUE))
  query
}