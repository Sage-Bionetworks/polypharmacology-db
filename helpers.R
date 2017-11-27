library(shiny)
library(rJava)
library(rcdk)
library(fingerprint)
library(enrichR)
library(webchem)
library(plyr)
library(dplyr)
library(tidyr)
library(grid)
library(plotly)
library(tibble)
library(heatmaply)

evo <- readRDS("Data/evotec_dgidb.RDS")
evo$Structure_ID <- as.character(evo$Structure_ID)
evo$Common_Name <- as.character(evo$Common_Name)
evo <- evo %>% filter(N_quantitative >= N_inactive | N_qualitative >= N_inactive | N_DGIDB > 0)

db.genes <- unique(evo$Hugo_Gene)

fp.evo <- readRDS("Data/fpevo.rds")[unique(evo$Original_molecule_SMILES)]
syns <- readRDS("Data/commname.rds")
syns$Original_molecule_SMILES <- as.character(syns$Original_molecule_SMILES)


##converts SMILES string to fingerprint
parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")

}

convertDrugToSmiles <- function(input) {
  filt <- filter(syns, Common_Name == input) %>% dplyr::select(Original_molecule_SMILES)
  filt
}

getTargetList <- function(selectdrugs) {
  targets <- filter(evo, Common_Name %in% selectdrugs) %>% dplyr::select(Common_Name, Hugo_Gene, MedianActivity_nM, N_quantitative, N_qualitative, 
                                                                          N_inactive, N_DGIDB, Confidence_Score) %>% arrange(-N_quantitative)
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

getMolImage <- function(input) {
  smi <- parse.smiles(input)
  view.image.2d(smi[[1]])
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
  mols <- filter(evo, Hugo_Gene == inp.gene) %>% 
    select(Structure_ID,Supplier_Molname, Common_Name, MedianActivity_nM, N_quantitative, N_qualitative, N_inactive, Original_molecule_SMILES) %>% 
    distinct()

}

getMolsFromGeneNetworks.edges <- function(inp.gene) {
  mols <- filter(evo, Hugo_Gene == inp.gene & N_quantitative > 1) %>% 
    select(Common_Name, N_quantitative) %>% distinct() %>% group_by(Common_Name) %>% top_n(5, N_quantitative) %>% ungroup()
  print(mols)
  
  net <- filter(evo, Common_Name %in% mols$Common_Name & N_quantitative >1)
  
  net$from <- as.character(net$Common_Name)
  net$to <- as.character(net$Hugo_Gene)
  net$width <- net$N_quantitative/10
  net$color <- "orange"
  net <- net %>% select(from, to, width, color)
  as.data.frame(net)
  print(net)
}

getMolsFromGeneNetworks.nodes <- function(inp.gene) {
  mols <- filter(evo, Hugo_Gene == inp.gene & N_quantitative > 1) %>% 
    select(Common_Name, N_quantitative) %>% 
    distinct() %>% group_by(Common_Name) %>% top_n(5, N_quantitative) %>% ungroup()
  
  net <- filter(evo, Common_Name %in% mols$Common_Name & N_quantitative >1) 
  
  id <- c(unique(as.character(net$Common_Name)), unique(as.character(net$Hugo_Gene)))
  label <- c(unique(as.character(net$Common_Name)), unique(as.character(net$Hugo_Gene)))
  color <- c(rep("blue", length(unique(as.character(net$Common_Name)))), 
             rep("green", length(unique(as.character(net$Hugo_Gene)))))
  
  net <- as.data.frame(cbind(id, label, color))
}

getSmiles <- function(input.name) {
  input.name <- input.name
  input.name <- URLencode(input.name)
  query <- as.vector(cir_query(input.name, representation = "smiles", first = TRUE))
  query
}

##get CTRP data for heatmap

drug.resp <- readRDS("Data/drugresp.rds")
ctrp.structures <- readRDS("Data/ctrpstructures.rds")
fp.ctrp <- readRDS("Data/fpctrp.rds")


plotSimCTRPDrugs <- function(input, sim.thres) {
  input <- input
  fp.inp <- parseInputFingerprint(input)
  
  sims <- lapply(fp.inp, function(i) {
    sim <- sapply(fp.ctrp, function(j) {
      distance(i, j)
    })
    bar <- as.data.frame(sim)
    bar$match <- rownames(bar)
    bar
  })
  sims <- ldply(sims)
  sims2 <- sims %>% filter(sim >= sim.thres) %>% arrange(-sim)
  sims2$cpd_smiles <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$sim, 3)
  drugs <- left_join(sims2, ctrp.structures) %>% dplyr::select(makenames, cpd_name, `Tanimoto Similarity`) %>% distinct()

  foo <- select(drug.resp, one_of(drugs$makenames)) %>% na.omit() %>% data.matrix(., rownames.force = T)
  #foo <- foo[complete.cases(foo),]
  print(foo)
  list(foo, drugs)

}

##same as above but for prepared sanger data

drug.resp.sang<-readRDS("Data/drugresp_sang.rds")
sang.structures<-readRDS("Data/sangstructures.rds")
fp.sang<-readRDS("Data/fpsang.rds")

plotSimSangDrugs <- function(input, sim.thres) {
  input <- input
  fp.inp <- parseInputFingerprint(input)
  
  sims <- lapply(fp.inp, function(i) {
    sim <- sapply(fp.sang, function(j) {
      distance(i, j)
    })
    bar <- as.data.frame(sim)
    bar$match <- rownames(bar)
    bar
  })
  sims <- ldply(sims)
  sims2 <- sims %>% filter(sim >= sim.thres) %>% arrange(-sim)
  sims2$smiles <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$sim, 3)
  drugs <- left_join(sims2, sang.structures) %>% dplyr::select(makenames, sanger_names, `Tanimoto Similarity`) %>% distinct()
  
  foo <- select(drug.resp.sang, one_of(drugs$makenames)) %>% na.omit() %>% data.matrix(., rownames.force = T)
  #foo <- foo[complete.cases(foo),]
  print(foo)
  list(foo, drugs)
}
