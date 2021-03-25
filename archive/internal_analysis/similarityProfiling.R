library(synapser)
library(plyr)
library(tidyverse)
library(rcdk)
library(fingerprint)
library(pbapply)
library(pheatmap)
library(tibble)
library(ggplot2)
library(ggrepel)
library(plotly)
library(diffusr)
library(annotables)
synLogin()

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
    fp.inp <- lapply(input.mol, get.fingerprint, type = "circular")
  }else{
    print('Please input a valid SMILES string.')
  }
}

fp.char <- sapply(fp.ncats, as.character)
mat <- 1-stringdistmatrix("110","1111", method = "jaccard", nthread=8, q = 1, useBytes = T) %>% t()


ncats.structures <- read.table(synGet("syn11559906")$path, sep = "\t", comment.char = "", header = T, quote = "") %>% 
  filter(smiles2 != "")

fp.ncats <- sapply(ncats.structures$smiles2, parseInputFingerprint)


distance.minified <- function(fp1,fp.list, method=c('tanimoto')){
  method <- match.arg(method)
  n <- length(fp1)
  f1 <- numeric(n)
  f2 <- numeric(n)
  f1[fp1@bits] <- 1
  
  sapply(fp.list, function(x){
    f2[x@bits] <- 1
    sim <- 0.0
    ret <-  .C("fpdistance", as.double(f1), as.double(f2),
             as.integer(n), as.integer(1),
             as.double(sim),
             PACKAGE="fingerprint")
    return (ret[[5]])
    })
}

fp <- c(rep(fp.ncats,200))

library(profvis)

profvis({

distances <- distance.minified(fp.ncats[[1]],fp)

sim <- sapply(fp, function(j) {
  distance(fp[[1]],j)
})
})

