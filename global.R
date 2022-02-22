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
library(shinyBS)
library(shinythemes)
library(visNetwork)
library(igraph)
library(shinyjs)
library(shinycssloaders)
library(flexdashboard)

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

parseInputFingerprint <- function(input, fp.type) {
  if(is.smiles(input)==TRUE){
    input.mol <- parse.smiles(as.character(input))
    lapply(input.mol, do.typing)
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
  filt
}

getTargetList <- function(selectdrugs) {
  
  targets <- db %>% 
    filter(inchikey %in% selectdrugs) %>% 
    as.data.frame() %>% 
    dplyr::select(inchikey, pref_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative, known_selectivity_index, confidence) %>% 
    arrange(-n_quantitative) 
  
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

getMolImage <- function(input) {
  smi <- parse.smiles(input)
  view.image.2d(smi[[1]])
}


getExternalDrugLinks <- function(inchikey_input) {
  links <- filter(db.links, inchikey %in% inchikey_input)
  links <- as.character(links$link)
  links <- paste(links, collapse = ", ")
}

getExternalGeneLinks <- function(gene) {
  links <- filter(db.gene.links, hugo_gene %in% gene)
  links <- as.character(links$link)
}

getNetwork <- function(drugsfound, selectdrugs) {
  targets <- drugsfound %>% 
    distinct() %>% filter(inchikey %in% selectdrugs)
  
  targets$from <- "input"
  targets$to <- as.character(targets$pref_name)
  targets$width <- ((targets$`Tanimoto Similarity`)^2) * 10
  targets$color <- "tomato"
  
  links <- sapply(selectdrugs, function(x){
    links <- getExternalDrugLinks(x)
  })
  
  targets$title <- links
  targets <- dplyr::select(targets, from, to, width, color, title)
}

getTargetNetwork <- function(selectdrugs, edge.size) {
  targets <- getTargetList(selectdrugs)
  targets$from <- targets$pref_name
  targets$to <- as.character(targets$hugo_gene)
  if(edge.size==TRUE){
    targets$width <- scales::rescale(targets$confidence, to = c(1,10))
  }
  if(edge.size==FALSE){
    targets$width <- 5
  }
  targets$color <- "tomato"

  targets <- dplyr::select(targets, from, to, width, color, inchikey) %>% 
    filter(from !="NA" & to != "NA")
}

dbs <- c("GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", 
         "KEGG_2019_Human")

getGeneOntologyfromTargets <- function(selectdrugs) {
  selectdrugs <- selectdrugs
  targets <- getTargetList(selectdrugs) %>% as.data.frame()
  target.list <- targets$hugo_gene
  
  if (length(target.list) > 0) {
    enriched <- enrichr(target.list, dbs)
  } else {
    print("no targets")
  }
  
}

getMolsFromGenes <- function(genes) {

  if(length(genes)>1){
  mols <- db %>% 
    filter(hugo_gene %in% genes) %>% 
    group_by(inchikey) %>% 
    mutate(count = n()) %>% 
    filter(count == length(genes)) %>% 
    ungroup() %>% 
    distinct() %>% 
    dplyr::select(-count) 
  }else{
    mols <- filter(db, hugo_gene == genes)
  }
  
  mols %>% 
    select(inchikey, pref_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative, known_selectivity_index, confidence) 
}


getMolsFromGeneNetworks.edges <- function(inp.gene, genenetmols, edge.size, gene.filter.metric) {
  mols <- genenetmols %>% top_n(15, !!sym(gene.filter.metric))
  
  net <- filter(db, inchikey %in% mols$inchikey) %>% distinct()
  
  net$from <- as.character(net$inchikey)
  net$to <- as.character(net$hugo_gene)
  if(edge.size==TRUE){
    net$width <- (net$confidence)/10
  }
  if(edge.size==FALSE){
    net$width <- 5
  }
  net$color <- "tomato"
  net <- net %>% dplyr::select(from, to, width, color)
  as.data.frame(net)
}


getMolsFromGeneNetworks.nodes <- function(inp.gene, genenetmols, gene.filter.metric) {
  mols <- genenetmols %>% top_n(15, !!sym(gene.filter.metric))

  net <- filter(db, inchikey %in% mols$inchikey) %>% 
      distinct() # %>% 
    # group_by(common_name) %>% 
    # top_n(20, confidence) %>% 
    # ungroup()
   
  id <- c(unique(as.character(net$inchikey)), 
          unique(as.character(net$hugo_gene)))
  label <- c(unique(as.character(net$pref_name)), 
             unique(as.character(net$hugo_gene)))
  color <- c(rep("blue", length(unique(as.character(net$pref_name)))), 
             rep("green", length(unique(as.character(net$hugo_gene)))))
  
  druglinks <- sapply(unique(as.character(net$inchikey)), function(x){
    druglinks <- getExternalDrugLinks(x)
  })

  genelinks <- sapply(unique(as.character(net$hugo_gene)), function(x){
    getExternalGeneLinks(x)
  })

  title <- c(druglinks, genelinks)

  net <- as.data.frame(cbind(id, label, color, title))
  
}

convert_id_to_structure_pubchem <- function(input_id, id_type = c("name", "inchikey"), output_type = c("InChI", "InChIKey", "CanonicalSMILES", "IsomericSMILES")){
  Sys.sleep(0.25) ##to prevent requests from happening too fast
  
  input <- URLencode(input_id)
  
  statement <- glue::glue('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{id_type}/{input}/property/{output_type}/xml')
  
  res <- httr::GET(statement)
  if(res$status_code==200){
    res_2 <- XML::xmlToList(rawToChar(res$content))
    struct <- purrr::pluck(res_2, "Properties", 2)
    if(is.null(struct)){struct <- NA}
  }else{
    message(glue::glue('input "{input_id}" appears to be invalid'))
    struct <- NA
  }
  
  struct
  
}

plotSimCTRPDrugs <- function(input, fp.type) {

  fp.inp <- parseInputFingerprint(input, fp.type = fp.type)
  
  if(fp.type == "circular"){fp.ctrp <- fp.ctrp.circular}
  if(fp.type == "extended"){fp.ctrp <- fp.ctrp.extended}
  if(fp.type == "maccs"){fp.ctrp <- fp.ctrp.maccs}
  
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
