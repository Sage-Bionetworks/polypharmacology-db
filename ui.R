library(shiny)
library(shinyBS)
library(shinythemes)
library(visNetwork)
library(igraph)
library(shinyjs)


shinyUI(
  fluidPage(
    useShinyjs(),
    div(id = "loading_page",
      img(src = "Ellipsis.svg"),
      h3("LOADING..."),
      align = "center"),
    tags$head(tags$style(type='text/css', "#loading_page { width:100%; margin-top: 20%;}")),
    hidden(
      div(id = "main_content",
  navbarPage("Polypharmacology DB", theme = shinytheme("flatly"),
  tabPanel("About",
           p("Welcome to the PPDB app, powered by Sage Bionetworks.
             The purpose of this app is to facilitate exploration of drug-target interaction databases.
             This app currently contains the Children's Tumor Foundation Drug-Target Database, licensed from Evotec, which
             summarizes activity data deposited in ChEMBL and inactivity data deposited in Pubchem."),
           br(),
           p(strong("How does PPDB work?")),
           p("PPDB leverages structural information of molecules and the associated target annotations to build a drug-target map 
             based on chemical similarity between molecules. PPDB includes drug-target interactions collated by Evotec, as well as a subset of those available in the DGIdb app.
             Examples of use-cases for this include:"),
           p(" - prediction of molecular targets for novel molecules based on structural similarity"),
           p(" - identification of off targets for molecules of interest"), 
           p(" - facilitating polypharmacologic drug discovery"),
           br(),
           p("Instructions:"),
           p("Click on the 'Molecules' tab to find targets associated with your molecule of interest."),
           p("Alternatively, if you have a target in mind, please enter the HUGO Gene Symbol on the 'Genes' tab."), 
           br(),
           p(strong("This app exists thanks to the following excellent R packages and organizations:")),
           br(),
           img(src='CTF_Logo.png'),
           br(),
           img(src= "sage_logo.png"),
           br(),
           tags$ul(
             tags$li("shiny"), 
             tags$li("shinyBS"),
             tags$li("shinythemes"),
             tags$li("rcdk"),
             tags$li("fingerprint"),
             tags$li("rJava"),
             tags$li("plyr"),
             tags$li("dplyr"),
             tags$li("DT"),
             tags$li("enrichR"),
             tags$li("webchem"),
             tags$li("visNetwork"),
             tags$li("igraph"))
           ), 
  
  tabPanel("Molecules",
  sidebarLayout(
    sidebarPanel(
      bsCollapse(
        bsCollapsePanel("Quick Start Guide",
                        p("You have three options for molecule lookup. 
                          You can 1) look up by molecule name in our database, 
                          2) search Pubchem for structures associated with your molecule name, or
                          3) directly enter the structure as represented by a SMILES string.
                          Then, pick a similarity using the slider, where 0 = highly dissimilar, and 1 = identical."), style = "primary"),
        bsCollapsePanel("Molecule Lookup",
                        fluidRow(selectizeInput("drugnames",
                                                choices = NULL,
                                                label = "",
                                                selected = "",
                                                width = "90%"),
                                 bsTooltip("drugnames", "Type molecule name here to search this database for SMILES strings.",
                                           "right", options = list(container = "body")), align = "center"),
                        fluidRow(actionButton("ppdbsearchbutton", "Find PPDB Mols", align = "center"), align = "center"),
                        div(),
                        p("Search this database for structures by compound name. Can't find what you're looking for? Move to the next panel to search Pubchem."), style = "warning"),
        bsCollapsePanel("Pubchem Search",
                 fluidRow(textInput("input.name",
                                    "input text", 
                                    label = NULL, 
                                    value = "lenalidomide",
                                    width = "90%"),
                                 bsTooltip("inpu", "Type molecule name here to search this database for SMILES strings.",
                                           "right", options = list(container = "body")), align = "center"),
                 fluidRow(actionButton("pubchembutton", "Find Pubchem Mols"), align = "center"), 
                 div(),
                 p("Input a compound name in this box to search PubChem's structure database."), style = "danger"),
        bsCollapsePanel("Direct Structure Input", fluidRow(textInput("smiles",
                                     "SMILES string", 
                                     label = "", 
                                     value = "Nc1cccc2C(=O)N(Cc12)C3CCC(=O)NC3=O",
                                     width = "90%"
                                     ),
                 bsTooltip("smiles", "Input the structural string (SMILES) here.",
                           "right", options = list(container = "body")), align = "center"),
                 div(),
                 p("Know your SMILES string already? Search by structure directly here."), style = "success"), open = "Quick Start Guide"),
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
              DT::dataTableOutput("kegg")))
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
  )
  )
  )
  )
  )
     