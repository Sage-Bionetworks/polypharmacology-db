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
           h4(strong("Welcome to the PPDB app, powered by Sage Bionetworks.")),
           h5("The purpose of this app is to facilitate exploration of drug-target interaction databases.
             This app currently contains the Children's Tumor Foundation Drug-Target Database, licensed from Evotec, which
             summarizes activity data deposited in ChEMBL and inactivity data deposited in Pubchem."),
           br(),
           h4(strong("How does PPDB work?")),
           p("PPDB leverages structural information of molecules and the associated target annotations to build a drug-target map 
             based on chemical similarity between molecules. PPDB includes drug-target interactions collated by Evotec, as well as a subset of those available in the DGIdb app.
             Examples of use-cases for this include:"),
           p(" - prediction of molecular targets for novel molecules based on structural similarity"),
           p(" - identification of off targets for molecules of interest"), 
           p(" - facilitating polypharmacologic drug discovery"),
           br(),
           h4(strong("Instructions:")),
           p("Click on the 'Molecules' tab to find targets associated with your molecule of interest."),
           p("Alternatively, if you have a target in mind, please enter the HUGO Gene Symbol on the 'Genes' tab."), 
           br(),
           h4(strong("This app exists thanks to the following excellent R packages and organizations:")),
           br(),
           fluidRow(
           img(src='CTF_Logo.png'),
           img(src= "sage_logo.png"), align = "center"),
           br(),
           fluidRow(
           h4("shiny"), 
             h4("shinyBS"),
             h4("shinythemes"),
             h4("rcdk"),
             h4("fingerprint"),
             h4("rJava"),
             h4("plyr"),
             h4("dplyr"),
             h4("DT"),
             h4("enrichR"),
             h4("webchem"),
             h4("visNetwork"),
             h4("igraph"), align = "center")
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
                        p("Search this database for structures by compound name. Can't find what you're looking for? Move to the next panel to search Pubchem."), style = "info"),
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
                 p("Input a compound name in this box to search PubChem's structure database."), style = "info"),
        bsCollapsePanel("Direct Structure Input", fluidRow(textInput("smiles",
                                     "SMILES string", 
                                     label = "", 
                                     value = "Nc1cccc2C(=O)N(Cc12)C3CCC(=O)NC3=O",
                                     width = "90%"
                                     ),
                 bsTooltip("smiles", "Input the structural string (SMILES) here.",
                           "right", options = list(container = "body")), align = "center"),
                 div(),
                 p("Know your SMILES string already? Search by structure directly here."), style = "info"), open = "Quick Start Guide"),
        bsCollapse(bsCollapsePanel("Similarity Threshold",  
      sliderInput("sim.thres", label = "", 
                             min=0.3, 
                             max=1,
                             value=0.90,
                             step = 0.01,
                             ticks = FALSE),
                 bsTooltip("sim.thres", "Set the Tanimoto similarity (1 being identical) here.",
                           "right", options = list(container = "body")), style = "warning"), open = "Similarity Threshold"),
                 uiOutput("sims"),
                 textOutput("pubchemsearchNA")), 
    mainPanel(
      tabsetPanel(
        tabPanel(title = img("Similar Molecules  ", id = "similarmolstab", src = "help.png", align = "right"),
          bsTooltip(id = "similarmolstab", title = "This table displays all similar molecules to the input molecule.", placement = "bottom", trigger = "hover"),
          DT::dataTableOutput("simmoltab")),
        tabPanel(title = img("Similarity Net  ", id = "similarmolnettab", src = "help.png", align = "right"),
          bsTooltip(id = "similarmolnettab", title = "This is a graphical representation of the previous table, where edge thickness is the similarity to the input.", placement = "bottom", trigger = "hover"),
          visNetworkOutput("net")),
        tabPanel(title = img("Targets  ", id = "targetstab", src = "help.png", align = "right"), 
          bsTooltip(id = "targetstab", title = "These are all of the targets associated with the input and similar molecules. Filter using the checkboxes in the sidebar.", placement = "bottom", trigger = "hover"),
          DT::dataTableOutput("value")),
        tabPanel(title = img("Target Net  ", id = "targetnettab", src = "help.png", align = "right"),
                 bsTooltip(id = "targetnettab", title = "This a network of all similar molecules (blue) and their targets (green).", placement = "bottom", trigger = "hover"),
          visNetworkOutput("targetnet")),
        tabPanel(title = img("Enrichr  ", id = "enrichrtab", src = "help.png", align = "right"),
                 bsTooltip(id = "enrichrtab", title = "This tab interactively queries Enrichr for enriched gene ontology and KEGG terms using your target list.", placement = "bottom", trigger = "hover"),
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
             sidebarPanel("Gene Input:",
                          fluidRow(textInput("inp.gene",
                                             "input text", 
                                             label = NULL, 
                                             value = "HDAC6",
                                             width = "90%"
                          )),
                          fluidRow(actionButton("genesearchbutton", "Search"), align = "center")
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
     