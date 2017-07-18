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

evo <- readRDS("Data/evotec.rds")
evo$Structure_ID <- as.character(evo$Structure_ID)
evo <- evo %>% filter(N_quantitative >= N_inactive | N_qualitative >= N_inactive)
fp.evo <- readRDS("Data/fpevo.rds")

parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
}

getTargetList <- function(input, sim.thres) {
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
  sims2 <- sims %>% filter(sim >= sim.thres)
  targets <- filter(evo, Structure_ID %in% sims2$match) %>% dplyr::select(Structure_ID, 
                                                                          Hugo_Gene, Protein_names, MedianActivity_nM, N_quantitative, N_qualitative, 
                                                                          N_inactive) %>% arrange(-N_quantitative)
  
  if (nrow(targets) > 1) {
    print(targets)
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
  sims2$Structure_ID <- as.character(sims2$match)
  targets <- left_join(sims2, evo) %>% dplyr::select(Structure_ID, Supplier_Molname, 
                                                     Common_Name, sim, Original_molecule_SMILES) %>% distinct()
}

getNetwork <- function(input, sim.thres) {
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
  sims2$Structure_ID <- as.character(sims2$match)
  targets <- left_join(sims2, evo) %>% dplyr::select(Structure_ID, sim) %>% 
    distinct()
  targets$from <- "input"
  targets$to <- as.character(targets$Structure_ID)
  targets$width <- ((targets$sim)^2) * 10
  targets$color <- "red"
  targets <- select(targets, from, to, width, color)
  print(targets)
}

getSimMolMatrix <- function(input, sim.thres) {
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
  sims2$Structure_ID <- as.character(sims2$match)
  
  sim.mols <- filter(evo, Structure_ID %in% sims2$Structure_ID)
  sim.mols <- unique(sim.mols$Original_molecule_SMILES)
  
  fp.sim <- parse.smiles(c(sim.mols, input))
  fp.sim <- lapply(fp.sim, get.fingerprint, type = "extended")
  
  mat <- fp.sim.matrix(fp.sim, method = "tanimoto")
  colnames(mat) <- names(fp.sim)
  rownames(mat) <- names(fp.sim)
  print(mat)
  plot <- pheatmap(mat)
  plot
}

getTargetNetwork <- function(input, sim.thres) {
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
  sims2$Structure_ID <- as.character(sims2$match)
  targets <- left_join(sims2, evo) %>% distinct()
  targets$from <- targets$Structure_ID
  targets$to <- as.character(targets$Uniprot_accession_numbers)
  targets$width <- 5
  targets$color <- "green"
  targets <- select(targets, from, to, width, color) %>% filter(from != 
                                                                  "NA" & to != "NA")
  # targets <- sample_n(targets, 50) ## this is temporary limit on number
  # of targets
  print(targets)
}

dbs <- c("GO_Molecular_Function_2017", "GO_Cellular_Component_2017", "GO_Biological_Process_2017", 
         "PPI_Hub_Proteins")

getGeneOntologyfromTargets <- function(input, sim.thres) {
  input <- input
  sim.thres <- sim.thres
  targets <- getTargetList(input, sim.thres)
  
  if (nrow(targets) > 1) {
    enriched <- enrichr(as.vector(targets$Hugo_Gene), dbs)
  } else {
    print("no targets")
  }
  
}


getMolsFromGenes <- function(inp.gene) {
  mols <- filter(evo, Hugo_Gene == inp.gene) %>% select(Structure_ID, 
                                                        Supplier_Molname, Common_Name, N_quantitative, N_qualitative, N_inactive) %>% 
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
  print(as.data.frame(net))
}

getMolsFromGeneNetworks.nodes <- function(inp.gene) {
  mols <- filter(evo, Hugo_Gene == inp.gene & N_quantitative > 1) %>% 
    select(Structure_ID, N_quantitative) %>% distinct()
  
  net <- filter(evo, Structure_ID %in% mols$Structure_ID & N_quantitative > 
                  1) %>% group_by(Structure_ID) %>% top_n(10, N_quantitative) %>% 
    ungroup()
  
  print(net)
  id <- c(unique(as.character(net$Structure_ID)), unique(as.character(net$Hugo_Gene)))
  label <- c(unique(as.character(net$Structure_ID)), unique(as.character(net$Hugo_Gene)))
  color <- c(rep("blue", length(unique(as.character(net$Structure_ID)))), 
             rep("green", length(unique(as.character(net$Hugo_Gene)))))
  
  net <- as.data.frame(cbind(id, label, color))
  print(net)
  
}

getSmiles <- function(input.name) {
  input.name <- input.name
  input.name <- URLencode(input.name)
  query <- cir_query(input.name, representation = "smiles", first = TRUE)
  query
}