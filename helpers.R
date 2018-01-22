library(shiny)
library(rJava)
library(rcdk)
library(fingerprint)
library(enrichR)
library(webchem)
library(plyr)
library(tidyverse)
library(plotly)
library(synapser)

is.smiles <- function(x, verbose = TRUE) { ##corrected version from webchem
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("rcdk needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # x <- 'Clc(c(Cl)c(Cl)c1C(=O)O)c(Cl)c1Cl'
  if (length(x) > 1) {
    stop('Cannot handle multiple input strings.')
  }
  out <- try(rcdk::parse.smiles(x), silent = TRUE)
  if (inherits(out[[1]], "try-error") | is.null(out[[1]])) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

synLogin()

db <- read.table(synGet("syn11681928")$path, header = T) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(internal_id, common_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative)

mini.db <- db %>% group_by(internal_id) %>% 
  mutate(n = sum(n_qualitative, n_quantitative, na.rm = T)) %>% 
  select(internal_id, common_name, n) %>% 
  distinct() %>% 
  filter(n>=10)

db.names <- read.table(synGet("syn11681849")$path, header = T)

db$internal_id <- as.character(db$internal_id)

fp.db <- readRDS(synGet("syn11693143")$path)
fp.db <- fp.db[names(fp.db) %in% unique(db$internal_id)]

fp.snappy <- fp.db[names(fp.db) %in% unique(mini.db$internal_id)]

db.genes <- unique(db$hugo_gene)

parseInputFingerprint <- function(input) {
  test_smiles <- is.smiles(input)
  if(is.smiles(input==TRUE)){
    input.mol <- parse.smiles(as.character(input))
    lapply(input.mol, do.typing)
    lapply(input.mol, do.aromaticity)
    lapply(input.mol, do.isotopes)
    fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
  }else{
    print('Please input a valid SMILES string.')
  }
}


convertDrugToSmiles <- function(input) {
  filt <- filter(db.names, common_name == input) %>% dplyr::select(smiles)
  filt
}

getTargetList <- function(selectdrugs) {
  targets <- filter(db, common_name %in% selectdrugs) %>% 
    arrange(-n_quantitative) %>%
    select(common_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative)
  
  if (nrow(targets) > 1) {
    targets
  } else {
    print("none found")
  }
}

getSimMols <- function(input, sim.thres, snappy) {
  input <- input
  fp.inp <- parseInputFingerprint(input)

  if(snappy == FALSE){
  sims <- lapply(fp.inp, function(i) {
    sim <- lapply(fp.db, function(j) {
      distance(i, j)
    })
    bar <- ldply(sim)
    colnames(bar) <- c("match", "similarity")
    bar
  })
  }

  if(snappy == TRUE){
    sims <- lapply(fp.inp, function(i) {
      sim <- lapply(fp.snappy, function(j) {
        distance(i, j)
      })
      bar <- ldply(sim)
      colnames(bar) <- c("match", "similarity")
      bar
    })
  }
  sims <- ldply(sims)
  sims2 <- sims %>% dplyr::filter(similarity >= sim.thres) %>% arrange(-similarity)
  sims2$internal_id <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$similarity, 3)
  targets <- left_join(sims2, db) %>% 
    dplyr::select(common_name, `Tanimoto Similarity`) %>% 
    distinct()
}

getMolImage <- function(input) {
  smi <- parse.smiles(input)
  view.image.2d(smi[[1]])
}

getNetwork <- function(drugsfound, selectdrugs) {
  targets <- drugsfound %>% 
    distinct() %>% filter(common_name %in% selectdrugs)
  targets$from <- "input"
  targets$to <- as.character(targets$common_name)
  targets$width <- ((targets$`Tanimoto Similarity`)^2) * 10
  targets$color <- "red"
  targets <- select(targets, from, to, width, color)
}

getTargetNetwork <- function(selectdrugs) {
  selectdrugs <- selectdrugs
  targets <- getTargetList(selectdrugs) %>% distinct() %>% filter(common_name %in% selectdrugs)
  targets$from <- targets$common_name
  targets$to <- as.character(targets$hugo_gene)
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
    enriched <- enrichr(as.vector(targets$hugo_gene), dbs)
  } else {
    print("no targets")
  }
  
}


getMolsFromGenes <- function(inp.gene) {
  mols <- filter(db, hugo_gene == inp.gene) %>% 
    distinct()

}

getMolsFromGeneNetworks.edges <- function(inp.gene) {
  mols <- filter(db, hugo_gene == inp.gene & n_quantitative > 1) %>% 
    select(common_name, n_quantitative) %>% distinct() %>% group_by(common_name) %>% top_n(5, n_quantitative) %>% ungroup()
  
  net <- filter(db, common_name %in% mols$common_name & n_quantitative >1)
  
  net$from <- as.character(net$common_name)
  net$to <- as.character(net$hugo_gene)
  net$width <- net$n_quantitative/10
  net$color <- "orange"
  net <- net %>% select(from, to, width, color)
  as.data.frame(net)
}

getMolsFromGeneNetworks.nodes <- function(inp.gene) {
  mols <- filter(db, hugo_gene == inp.gene & n_quantitative > 1) %>% 
    select(common_name, n_quantitative) %>% 
    distinct() %>% group_by(common_name) %>% top_n(5, n_quantitative) %>% ungroup()
  
  net <- filter(db, common_name %in% mols$common_name & n_quantitative >1) 
  
  id <- c(unique(as.character(net$common_name)), unique(as.character(net$hugo_gene)))
  label <- c(unique(as.character(net$common_name)), unique(as.character(net$hugo_gene)))
  color <- c(rep("blue", length(unique(as.character(net$common_name)))), 
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

plotSimCTRPDrugs <- function(input) {
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
  sims2 <- sims %>% arrange(-sim)
  sims2$cpd_smiles <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$sim, 3)
  drugs <- left_join(sims2, ctrp.structures) %>% dplyr::select(makenames, cpd_name, `Tanimoto Similarity`) %>% distinct()
  
  top_drug <- top_n(drugs, 1, `Tanimoto Similarity`)
  
  drug.resp.single <- drug.resp[[top_drug$makenames]]

    cors<-sapply(colnames(drug.resp), function(x){
    test <- data.frame(drug.resp.single, drug.resp[[x]])
    if(nrow(test[complete.cases(test),])>1){
      cor<-cor.test(drug.resp.single, drug.resp[[x]], method = "spearman", use = "complete.obs")
      res <- c("p.val" = cor$p.value, cor$estimate)
    }else{
      res <- c("p.val" = -1, "rho" = 0)
    }
  })

  cors <- cors %>% 
    t() %>%
    as.data.frame() %>% 
    rownames_to_column("makenames") %>%
    inner_join(drugs) %>% 
    filter(p.val != -1)
  
  cors$Correlation <- cors$rho

  cors$`BH adj p.val` <- p.adjust(cors$p.val, method = "BH")
  cors
}

##same as above but for prepared sanger data

drug.resp.sang<-readRDS("Data/drugresp_sang.rds")
sang.structures<-readRDS("Data/sangstructures.rds")
fp.sang<-readRDS("Data/fpsang.rds")


plotSimSangDrugs <- function(input) {
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
  sims2 <- sims %>% arrange(-sim)
  sims2$smiles <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$sim, 3)
  drugs <- left_join(sims2, sang.structures) %>% dplyr::select(makenames, sanger_names, `Tanimoto Similarity`) %>% distinct()
  
  top_drug <- top_n(drugs, 1, `Tanimoto Similarity`)
  
  drug.resp.single <- drug.resp.sang[[top_drug$makenames]]
    
  cors<-sapply(colnames(drug.resp.sang), function(x){
    test <- data.frame(drug.resp.single, drug.resp.sang[[x]])
    if(nrow(test[complete.cases(test),])>1){
      cor<-cor.test(drug.resp.single, drug.resp.sang[[x]], method = "spearman", use = "complete.obs")
      res <- c("p.val" = cor$p.value, cor$estimate)
    }else{
      res <- c("p.val" = -1, "rho" = 0)
    }
  })
  
  cors <- cors %>% 
    t() %>%
    as.data.frame() %>% 
    rownames_to_column("makenames") %>%
    inner_join(drugs) %>% 
    filter(p.val != -1)
  
  cors$Correlation <- cors$rho
  
  cors$`BH adj p.val` <- p.adjust(cors$p.val, method = "BH")
  cors
}
