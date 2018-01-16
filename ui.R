library(shiny)
library(shinyBS)
library(shinythemes)
library(visNetwork)
library(igraph)
library(shinyjs)
library(plotly)

shinyUI(
  fluidPage( 
    tags$head(includeScript("https://www.googletagmanager.com/gtag/js?id=UA-109127366-1"),
              includeScript("www/google_analytics.js")),
    useShinyjs(),
    div(id = "loading_page",
      img(src = "Ellipsis.svg"),
      h3("LOADING..."),
      align = "center"),
    tags$head(tags$style(type='text/css', "#loading_page { width:100%; margin-top: 20%;}")),
    hidden(
      div(id = "main_content",
  navbarPage("Drug-Target Explorer", theme = shinytheme("flatly"),
  tabPanel("About",
           h4(strong("Welcome to the Drug-Target Explorer, powered by Sage Bionetworks.")),
           h5("The purpose of this app is to facilitate exploration of drug-target interaction databases.
             This app currently contains the Children's Tumor Foundation Drug-Target Database, licensed from Evotec, which
             summarizes activity data deposited in ChEMBL and inactivity data deposited in Pubchem."),
           br(),
           h4(strong("How does this app work?")),
           h5("PPDB leverages structural information of molecules and the associated target annotations to build a drug-target map 
             based on chemical similarity between molecules. PPDB includes drug-target interactions collated by Evotec, as well as a subset of those available in the DGIdb app.
             Examples of use-cases for this include:"),
           h5(" - prediction of molecular targets for novel molecules based on structural similarity"),
           h5(" - identification of off targets for molecules of interest"), 
           h5(" - facilitating polypharmacologic drug discovery"),
           br(),
           h4(strong("Instructions:")),
           p("Click on the 'Molecules' tab to find targets associated with your molecule of interest."),
           p("Alternatively, if you have a target in mind, please enter the HUGO Gene Symbol on the 'Genes' tab."), 
           br(),
           fluidRow(
           img(src='CTF_Logo.png'),
           img(src= "sage_logo.png"), align = "center"),
           br(),
           h4(strong("This app relies on the following excellent R packages:")),
           fluidRow(
           h5("shiny, shinyBS, shinythemes, rcdk, fingerprint, rJava"),
           h5("plyr, dplyr, DT, enrichR, webchem, visNetwork, igraph"), align = "center")), 
  
  tabPanel("Molecules",
  sidebarLayout(
    sidebarPanel(
      bsCollapse(
        bsCollapsePanel("Quick Start Guide",
                        p("You have three options for molecule lookup. 
                          You can:"),
                          div(),
                          p("1) look up by molecule name in our database,"),
                          div(),
                          p("2) search Pubchem for structures associated with your molecule name, or"),
                          div(),
                          p("3) directly enter the structure as represented by a SMILES string."),
                          div(),
                          p(strong("After entering a molecule, pick a similarity using the slider, where 0 = highly dissimilar, and 1 = identical.")),
                          div(),
                          checkboxInput("snappy", label = "Snappy mode? Uses smaller dataset.", value = TRUE), style = "primary"),
        bsCollapsePanel("Molecule Lookup",
                        fluidRow(selectizeInput("drugnames",
                                                choices = NULL,
                                                label = "",
                                                selected = "",
                                                multiple = FALSE,
                                                width = "90%"),
                                 bsTooltip("drugnames", "Type molecule name here to search this database for SMILES strings.",
                                           "right", options = list(container = "body")), align = "center"),
                        fluidRow(actionButton("ppdbsearchbutton", "Find PPDB Mols", align = "center"), align = "center"),
                        br(),
                        p("Search this database for structures by compound name. Can't find what you're looking for? Move to the next panel to search Pubchem."), style = "info"),
        bsCollapsePanel("Pubchem Search",
                 fluidRow(textInput("input.name",
                                    "input text", 
                                    label = NULL, 
                                    value = "lenalidomide",
                                    width = "90%"),
                                 bsTooltip("input.name", "Type molecule name here to search this database for SMILES strings.",
                                           "right", options = list(container = "body")), align = "center"),
                 fluidRow(actionButton("pubchembutton", "Find Pubchem Mols"), align = "center"), 
                 br(),
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
                 textOutput("pubchemsearchNA"),
      bsCollapse(
      bsCollapsePanel("2D Structure (input)", fluidRow(imageOutput("structureimage"), align = "center"), style = "danger"), open = "2D Structure (input)")
      ),
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
                 bsTooltip(id = "enrichrtab", title = "This tab interactively queries Enrichr for enriched gene ontology and KEGG terms using your target list. May take a few seconds, please be patient!", placement = "bottom", trigger = "hover"),
            tabsetPanel(
                 tabPanel("GO Molecular Function",
              DT::dataTableOutput("GOMF.mol")),
            tabPanel("GO Cellular Component",
              DT::dataTableOutput("GOCC.mol")),
            tabPanel("GO Biological Process",
              DT::dataTableOutput("GOBP.mol")),
            tabPanel("KEGG Pathways",
              DT::dataTableOutput("kegg")))),
        tabPanel(title = img("CCLE  ", id = "ccletab", src = "help.png", align = "right"),
                 bsTooltip(id = "ccletab", title = "This tab searches CCLE drug response data for drugs related to your query compound.", placement = "bottom", trigger = "hover"),
                 h4(textOutput("ccle_mol")),
                 div(),
                 plotlyOutput("ccle_2", height = "400px", width = "90%"),
                 plotlyOutput("ccle_1", height = "400px", width = "90%")),
        tabPanel(title = img("Sanger  ", id = "sangertab", src = "help.png", align = "right"),
                 bsTooltip(id = "sangertab", title = "This tab searches Sanger cell line response data for drugs related to your query compound.", placement = "bottom", trigger = "hover"),
                 h4(textOutput("sang_mol")),
                 div(),
                 plotlyOutput("sang_2", height = "400px", width = "90%"),
                 plotlyOutput("sang_1", height = "400px", width = "90%"))
                 
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
                          ), actionButton("genebutton", "Search", align = "center"), align = "center")
                          ),
             mainPanel(
               tabsetPanel(
                 tabPanel(title = img("Matched Molecules  ", id = "matchmolstab", src = "help.png", align = "right"),
                          bsTooltip(id = "matchmolstab", title = "This tab displays molecules with a known interaction with the input gene", placement = "right", trigger = "hover"),
                          DT::dataTableOutput("genetargets")),
                 tabPanel(title = img("Target Network  ", id = "targnetworktab", src = "help.png", align = "right"),
                            bsTooltip(id = "targnetworktab", title = "This tab graphically displays the input gene, associated drugs, and the other targets associated with those drugs.", placement = "right", trigger = "hover"),
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
     