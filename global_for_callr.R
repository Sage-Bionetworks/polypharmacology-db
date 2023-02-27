#
# This file contains the parts of `global.R` which use `rcdk` and which will be
# invoked from a separate process using `callr`.
#

options(java.parameters = "-Xmx8g" ) 
options(shiny.trace = TRUE)
library(shiny)
library(DT)
library(png)
library(rJava)
library(rcdk)
library(fingerprint)
library(enrichR)
library(webchem)
library(plyr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(plotly)
library(shinyBS)
library(shinythemes)
library(visNetwork)
library(igraph)
library(shinyjs)
library(rjson)
library(shinycssloaders)
library(conflicted)
library(callr)
conflict_prefer("filter", "dplyr")
conflict_prefer("is.connected", "rcdk")
conflict_prefer("count", "fingerprint")
conflict_prefer("renderDataTable", "DT")
conflict_prefer("arrange", "dplyr")
conflict_prefer("mutate", "dplyr")

source("helpers_for_callr.R")

smiles_parser <- get.smiles.parser()

parse_smiles<-function(input) {
    tryCatch(rcdk::parse.smiles(input), error = function(e) {
		stop(sprintf("In 'parse_smiles' rcdk::parse.smiles raised an exception: %s", e))
  	})
}

is.smiles <- function(x, verbose = TRUE) { ##corrected version from webchem
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("rcdk needed for this function to work. Please install it.",
         call. = FALSE)
  }
  if (!is.character(x)) {
  	return(FALSE)
  }
  if (length(x) > 1) {
    stop('Cannot handle multiple input strings.')
  }
  out <- tryCatch(parse_smiles(x), error = function(e) {
	message(e)
	return(FALSE)
  })
  if (is.null(out[[1]])) {
    return(FALSE)
  } else {
  	return(TRUE)
  }
}

parseInputFingerprint <- function(input, fp.type) {
  if(is.smiles(input)==TRUE){
    input.mol <- parse_smiles(input)
    if (is.null(input.mol[[1]])) {
    	stop("rcdk::parse.smiles failed to return a result.")
    }
    lapply(input.mol, set.atom.types)
    lapply(input.mol, do.aromaticity)
    lapply(input.mol, do.isotopes)
    fp.inp <- lapply(input.mol, get.fingerprint, type = fp.type)
  }else{
    print('Please input a valid SMILES string.')
  }
}

distance.minified <- function(fp1,fp.list){ #this function is a stripped down fingerprint::distance that runs about 2-3x faster; big gains for the app, but not as feature rich
  n <- length(fp1)
  f1 <- numeric(n)
  f2 <- numeric(n)
  f1[fp1@bits] <- 1
  
  sapply(fp.list, function(x){ 
    f2[x@bits] <- 1
    sim <- 0.0
    ret <- .C("fpdistance", as.double(f1), as.double(f2),
              as.integer(n), as.integer(1),
              as.double(sim),
              PACKAGE="fingerprint")
    return (ret[[5]])
  })
}

convertDrugToSmiles <- function(input) {
  filt <- filter(db.names, synonym == input) %>% 
    dplyr::select(inchikey) %>% 
    dplyr::inner_join(db_structures) %>% 
    dplyr::select(std_smiles)
  as.character(filt[1,1])
}

similarityFunction <- function(input, fp.type) {
  input <- input
  fp.type <- fp.type
  fp.inp <- parseInputFingerprint(input, fp.type)
  
  if(fp.type=="extended"){ sim <- distance.minified(fp.inp[[1]], fp.extended) }
  if(fp.type=="circular"){ sim <- distance.minified(fp.inp[[1]], fp.circular) }
  if(fp.type=="maccs"){ sim <- distance.minified(fp.inp[[1]], fp.maccs) }
  # if(fp.type=="kr"){ sim <- distance.minified(fp.inp[[1]], fp.kr) }
  # if(fp.type=="pubchem"){ sim <- distance.minified(fp.inp[[1]], fp.pubchem) 
  
  bar <- enframe(sim) %>% 
    set_names(c("match", "similarity")) %>% 
    top_n(50, similarity) ##hard cutoff to avoid overloading the app - large n of compounds can cause sluggish response wrt visualizations
  
}

getSimMols <- function(sims, sim.thres) {
  sims2 <- sims %>% dplyr::filter(similarity >= sim.thres) %>% arrange(-similarity)
  sims2$inchikey <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$similarity, 3)
  targets <- left_join(sims2, db) %>% 
    dplyr::select(inchikey, pref_name, `Tanimoto Similarity`) %>% 
    distinct() %>% 
    as.data.frame()
  
}

getMolImage <- function(input, outfile) {
  smi <- parse_smiles(input)
  img <- view.image.2d(smi[[1]])
  writePNG(img, target = outfile, dpi = 600)
}



plotSimSangDrugs <- function(input, fp.type) {
  
  fp.inp <- parseInputFingerprint(input, fp.type = fp.type)
  
  if(fp.type == "circular"){fp.sang <- fp.sang.circular}
  if(fp.type == "extended"){fp.sang <- fp.sang.extended}
  if(fp.type == "maccs"){fp.sang <- fp.sang.maccs}
  
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
