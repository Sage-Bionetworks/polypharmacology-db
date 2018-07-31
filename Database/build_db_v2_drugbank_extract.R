library(XML)
library(biomaRt)
library(synapser)
synLogin()

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/polypharmacology-db/master/build_db_v2_drugbank_extract.R"
db <- synGet("syn12685088")

drugbank<-xmlToList(db$path) ##takes a bit of time to read in (20-30min?)

drugbank.summary <- lapply(1:(length(drugbank)-1), function(x){
  print(x)
  targets <- c()
  for(i in 1:length(drugbank[[x]]$targets)){
    targets<-c(targets, drugbank[[x]]$targets[[i]]$polypeptide$`external-identifiers`$`external-identifier`$identifier)
  }
  # 
  # hgnc <- c()
  # for(i in 1:length(drugbank[[x]]$targets)){
  #   if(!is.null(drugbank[[x]]$targets[[i]]$polypeptide$`external-identifiers`$`external-identifier`$identifier)){
  #     hgnc<-c(hgnc, drugbank[[x]]$targets[[i]]$polypeptide$`external-identifiers`$`external-identifier`$identifier)
  #   }else{
  #     hgnc<-c(hgnc,"null")
  #   }
  # }

  drug <- c(rep(drugbank[[x]]$`drugbank-id`$text, length(targets)))
  
  data.frame(drug, targets)
})

drugbank.df <- ldply(drugbank.summary) 

ensembl = useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")

map <- getBM(attributes = c("hgnc_id", "hgnc_symbol"), filters = "hgnc_id", values = unique(drugbank.df$targets), mart = ensembl) %>% 
  set_names(c("targets", "hugo_gene"))

drugbank.df <- inner_join(drugbank.df, map)

write.table(drugbank.df, "NoGit/drugbank_5-1-0_targets.txt", sep = "\t")

synStore(File("NoGit/drugbank_5-1-0_targets.txt", parentId = "syn12683828"),
         executed = this.file, used = "syn12685088")
