library(visNetwork)
library(shinythemes)
library(igraph)

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
  tabPanel("SMILES Lookup",
           
           sidebarLayout(
             sidebarPanel("Pubchem lookup",
                          fluidRow(textInput("input.name",
                                             "input text", 
                                             label = NULL, 
                                             value = "lenalidomide",
                                             width = "100%"
                          )),
                          div(),
                          p("Database lookup"),
                          fluidRow(selectizeInput("drugnames",
                                                  choices = NULL,
                                                  label = "",
                                                  selected = ""))),
             mainPanel(p("Pubchem results:"),
                       fluidRow(textOutput("smileslookup")),
                       p("Database results:"),
                       fluidRow(textOutput("smileslookup2"))))),
  tabPanel("Molecules",
  sidebarLayout(
    sidebarPanel("Smiles Input",
        fluidRow(textInput("smiles",
                                     "input text", 
                                     label = "", 
                                     value = "CC(C)(C)C1=CN=C(CSC2=CN=C(NC(C3CCNCC3)=O)S2)O1",
                                     width = "100%"
                                     ),
                 sliderInput("sim.thres", "Similarity Threshold", 
                             min=0.3, 
                             max=1,
                             value=0.90,
                             step = 0.01,
                             ticks = FALSE),
                 uiOutput("sims"))), 
    mainPanel(
      tabsetPanel(
        tabPanel(
          strong("Similar Molecules"),
          DT::dataTableOutput("simmoltab")),
        tabPanel(
              strong("Targets"),
              DT::dataTableOutput("value")),
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
              DT::dataTableOutput("kegg")))),
        tabPanel(
              strong("Similarity Net"),
              visNetworkOutput("net")),
        tabPanel(
          strong("Target Net"),
          visNetworkOutput("targetnet"))
       # tabPanel(
      #        strong("Matrix"),
       #       visNetworkOutput("simmat"))
      #  )
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
                                             width = "100%"
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
     