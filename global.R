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

sessionInfo()

loading <- function() {
  shinyjs::hide("loading_page")
  shinyjs::show("main_content")
}

#-----------------------------------------------------------------------------------------

withRetries<-function(fcn, context) {
	cntr<-0
	backoff<-as.difftime("00:00:01") # one second
	startTime<-Sys.time()
	
	repeat {
		cntr<-cntr+1
		fcnResult<-try(fcn())
		
		# retry
		if (is(fcnResult, "try-error")) {
			reportableResult<-fcnResult[[1]]
		} else {
			return(fcnResult)
		}
		sleepTime<-sleepTime(startTime, Sys.time(), backoff)
		if (sleepTime>0) {
			message(sprintf("%s Error encountered: %s. Will wait for %.1f seconds then retry. Press CTRL+C to quit.", 
							context, reportableResult, sleepTime))
			Sys.sleep(sleepTime)
		}
		backoff <- increaseBackoff(backoff)
		if (maxWaitTimeExceeded(startTime, Sys.time())) break
	} # end 'repeat' loop
	stop("Max backoff exceeded.")
}

MAX_WAIT_TIME<-as.difftime("00:05:00")
BACKOFF_MULTIPLIER <- 2
MAX_BACKOFF<-as.difftime("00:01:00")

maxWaitTimeExceeded<-function(startTime, currentTime) {
	currentTime-startTime>=MAX_WAIT_TIME
}

sleepTime<-function(startTime, currentTime, backoff) {
	maxWaitTimeRemaining<-max(as.difftime("00:00:00"),MAX_WAIT_TIME-(currentTime-startTime))
	min(backoff, maxWaitTimeRemaining)
}

increaseBackoff<-function(currentBackOffSeconds) {
	min(MAX_BACKOFF, BACKOFF_MULTIPLIER*currentBackOffSeconds)
}

#-----------------------------------------------------------------------------------------

remote_session<-callr::r_session$new(wait=TRUE)

call_remote_function<-function(function_name, arguments=list()) {
	if (FALSE) { # retries and debug statements
		fcn<-function() {
			remote_session$run(do.call, list(what=function_name, args=arguments))
		}
		message(sprintf("About to call %s.  Remote session has state %s.", function_name, remote_session$get_state()))
		result<-withRetries(fcn, function_name)
		message(sprintf("Done calling %s.  Remote session has state %s.", function_name, remote_session$get_state()))
		result
	} else {
		remote_session$run(do.call, list(what=function_name, args=arguments))
	}
}

message("Loading global_for_callr.R into remote session...")
remote_session$run(function(){source("global_for_callr.R")}, list())
message("...done.")


is.smiles <- function(x, verbose = TRUE) {
	call_remote_function("is.smiles", list(x, verbose))
}

convertDrugToSmiles <- function(input) {
  call_remote_function("convertDrugToSmiles", list(input))
}

getTargetList <- function(selectdrugs) {
  targets <- db %>% 
    filter(inchikey %in% selectdrugs) %>% 
    as.data.frame() %>% 
    dplyr::select(inchikey, pref_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative, known_selectivity_index, confidence) %>% 
    arrange(-n_quantitative) 
  
}

similarityFunction <- function(input, fp.type) {
  call_remote_function("similarityFunction", list(input, fp.type))
}

getSimMols <- function(sims, sim.thres) {
  call_remote_function("getSimMols", list(sims, sim.thres))
}

getMolImage <- function(input, outfile) {
  call_remote_function("getMolImage", list(input, outfile))
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


plotSimSangDrugs <- function(input, fp.type) {
  call_remote_function("plotSimSangDrugs", list(input, fp.type))
}
