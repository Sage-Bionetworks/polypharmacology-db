library(visNetwork)
library(shinythemes)
library(igraph)
library(shinyBS)

#syns <- readRDS("Data/commname.rds")

shinyUI(navbarPage("Polypharmacology DB", 
                   theme = shinytheme("flatly"),
  tabPanel("About",
           p("Welcome to the PPDB app, powered by Sage Bionetworks.
             The purpose of this app is to facilitate exploration of drug-target interaction databases.
             This app currently contains the Children's Tumor Foundation Drug-Target Database, licensed from Evotec, which
             summarizes activity data deposited in ChEMBL and inactivity data deposited in Pubchem."),
           br(),
           p(strong("How does PPDB work?")),
           p("PPDB leverages structural information of molecules and the associated target annotations to build a drug-target map 
             based on chemical similarity between molecules. Examples of use-cases for this include:"),
           p(" - prediction of molecular targets for novel molecules based on structural similarity"),
           p(" - identification of off targets for molecules of interest"), 
           p(" - facilitating polypharmacologic drug discovery"),
           br(),
           p("To begin, please click on 'SMILES lookup' to find the SMILES string for your chemical of interest,
             or look it up in an external database."),
           p("Then enter that SMILES string on the 'Molecules' tab to find related molecules in the database."),
           p("Alternatively, if you have a target in mind, please enter the HUGO Gene Symbol on the 'Genes' tab.")
           ),
  tabPanel("Molecules",
  sidebarLayout(
    sidebarPanel(
      bsCollapse(
        bsCollapsePanel("Pubchem Search",
                 fluidRow(textInput("input.name",
                                    "input text", 
                                    label = NULL, 
                                    value = "lenalidomide",
                                    width = "90%"),
                                 bsTooltip("inpu", "Type molecule name here to search this database for SMILES strings.",
                                           "right", options = list(container = "body")), align = "center"),
                 fluidRow(actionButton("pubchembutton", "Find Pubchem Mols"), align = "center"), style = "danger"),
        bsCollapsePanel("PPDB Lookup",
                 fluidRow(selectizeInput("drugnames",
                                         choices = NULL,
                                         label = "",
                                         selected = "",
                                         width = "90%"),
                          bsTooltip("drugnames", "Type molecule name here to search this database for SMILES strings.",
                                    "right", options = list(container = "body")), align = "center"),
                 fluidRow(actionButton("ppdbsearchbutton", "Find PPDB Mols", align = "center"), align = "center"), style = "warning"),
        bsCollapsePanel("Direct Structure Input",fluidRow(textInput("smiles",
                                     "SMILES string", 
                                     label = "", 
                                     value = "Nc1cccc2C(=O)N(Cc12)C3CCC(=O)NC3=O",
                                     width = "90%"
                                     ),
                 bsTooltip("smiles", "Input the structural string (SMILES) here.",
                           "right", options = list(container = "body")), align = "center"), style = "success"), open = "Pubchem Search"),
                 sliderInput("sim.thres", "Similarity Threshold", 
                             min=0.3, 
                             max=1,
                             value=0.90,
                             step = 0.01,
                             ticks = FALSE),
                 bsTooltip("sim.thres", "Set the Tanimoto similarity (1 being identical) here.",
                           "right", options = list(container = "body")),
                 uiOutput("sims"),
                 textOutput("pubchemsearchNA")), 
    mainPanel(
      tabsetPanel(
        tabPanel(
          strong("Similar Molecules"),
          DT::dataTableOutput("simmoltab")),
        tabPanel(
          strong("Similarity Net"),
          visNetworkOutput("net")),
        tabPanel(
              strong("Targets"),
              DT::dataTableOutput("value")),
        tabPanel(
          strong("Target Net"),
          visNetworkOutput("targetnet")),
        tabPanel(strong("Enrichr"),
          tabsetPanel(
            tabPanel(
              strong("GO Molecular Function"),
              DT::dataTableOutput("GOMF.mol")),
            tabPanel(
              strong("GO Cellular Component"),
              DT::dataTableOutput("GOCC.mol")),
            tabPanel(
              strong("GO Biological Process"),
              DT::dataTableOutput("GOBP.mol")),
            tabPanel(
              strong("KEGG Pathways"),
              DT::dataTableOutput("kegg"))))
        )
    )
  )
  ),
  tabPanel("Genes",
           sidebarLayout(
             sidebarPanel("Genes",
                          fluidRow(textInput("inp.gene",
                                             "input text", 
                                             label = NULL, 
                                             value = "HDAC6",
                                             width = "90%"
                          ))
                          ),
             mainPanel(
               tabsetPanel(
                 tabPanel(
                   strong("Matched Mols"),
                   dataTableOutput("genetargets")),
                 tabPanel(
                   strong("Target Network"),
                   visNetworkOutput("genetargetnet"))
               )
             )
           )
         )
))
     