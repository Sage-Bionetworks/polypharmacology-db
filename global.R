options(java.parameters = "-Xmx8g" ) 
library(shiny)
library(DT)
library(png)
library(rJava)
library(rcdk)
library(fingerprint)
library(enrichR)
library(webchem)
library(plyr)
library(tidyverse)
library(plotly)
library(synapser)
library(shinyBS)
library(shinythemes)
library(visNetwork)
library(igraph)
library(shinyjs)
library(shinycssloaders)

loading <- function() {
  shinyjs::hide("loading_page")
  shinyjs::show("main_content")
}


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

similarityFunction <- function(input, snappy) {
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
}

getSimMols <- function(sims, sim.thres) {
  sims2 <- sims %>% dplyr::filter(similarity >= sim.thres) %>% arrange(-similarity)
  sims2$internal_id <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$similarity, 3)
  targets <- left_join(sims2, db) %>% 
    dplyr::select(common_name, `Tanimoto Similarity`) %>% 
    distinct() %>% 
    as.data.frame()
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
  genes <- trimws(unlist(strsplit(inp.gene,",")))
  if(length(genes)>1){
  mols <- db %>% 
    mutate(hugo_gene = as.character(hugo_gene)) %>% 
    mutate(keep = hugo_gene %in% genes) %>% 
    group_by(internal_id, keep) %>% 
    mutate(count = n()) %>% 
    filter(keep == TRUE, count >= length(genes)) %>% 
    ungroup() %>% 
    distinct() %>% 
    select(-keep, -count)
  }
  if(length(genes)==1){
    mols <- filter(db, hugo_gene == inp.gene)
  }
}

getMolsFromGeneNetworks.edges <- function(inp.gene, genenetmols) {
  mols <- genenetmols %>% top_n(10, confidence)
  
  net <- filter(db, common_name %in% mols$common_name)
  
  net$from <- as.character(net$common_name)
  net$to <- as.character(net$hugo_gene)
  net$width <- net$n_quantitative/10
  net$color <- "black"
  net <- net %>% select(from, to, width, color)
  as.data.frame(net)
}

getMolsFromGeneNetworks.nodes <- function(inp.gene, genenetmols) {
  mols <- genenetmols %>% top_n(10, confidence)
  
  net <- filter(db, common_name %in% mols$common_name) 
  
  id <- c(unique(as.character(net$common_name)), unique(as.character(net$hugo_gene)))
  label <- c(unique(as.character(net$common_name)), unique(as.character(net$hugo_gene)))
  color <- c(rep("blue", length(unique(as.character(net$common_name)))), 
             rep("green", length(unique(as.character(net$hugo_gene)))))
  
  net <- as.data.frame(cbind(id, label, color))
}

getSmiles <- function(input.name) {
  input.name <- input.name
  input.name <- URLencode(input.name)
  query <- as.vector(cir_query(input.name, representation = "smiles", first = TRUE))
  query
}

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
