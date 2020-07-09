
options(spinner.color="#0dc5c1", spinner.color.background="#FFFFFF", spinner.type = 2)

shinyUI(
  fluidPage( 
    tags$head(
      singleton(
        includeScript("www/readCookie.js")
      )
    ),
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
  navbarPage("Drug-Target EXplorer", theme = shinytheme("flatly"),
  tabPanel("About",
           h4(strong("Welcome to the Drug-Target Explorer.")),
           h5("The purpose of this app is to facilitate exploration of drug-target interaction databases.
             The underlying database for this app is a harmonized dataset which summarizes quantitative and qualitative 
             small molecule activity data for human targets from ChEMBL, the Drug-Gene Interaction Database, DrugBank, ChemicalProbes.org, and Klaeger et al, Science, 2017.
             The current database (v2) contains summarizes evidence for >507,000 small-molecule-target interactions (>304,000 chemical entities and 3650 human targets)."),
           br(),
           h4(strong("How does this tool work?")),
           h5("This tool leverages structural information of molecules and the associated target annotations to build a drug-target map 
             based on chemical similarity between molecules.
             Examples of use-cases for this include:"),
           h5(" - prediction of molecular targets for novel molecules based on structural similarity"),
           h5(" - identification of off targets for molecules of interest"), 
           h5(" - facilitating polypharmacologic drug discovery"),
           h5(" - identifying Gene Ontology terms and KEGG pathways associated with small molecule target lists"),
           br(),
           h4(strong("How do I use the app to search by chemical entity?")),
           p("Click on the 'Molecules' tab to find targets associated with your molecule of interest. 
             Search for the drug using the 'Molecule Lookup' panel on the left. 
             Alternatively, use the 'CIR Search' panel to use the NCI CACTUS Chemical Identifier Resolver, or directly enter the SMILES structural string in the 'Direct Structure Input' panel.
             Set the Tanimoto similarity threshold using the slider, where 1 represents an identical molecular fingerprint."),
           h4(strong("How do I use the app to search by target?")),
           p("If you have a target in mind, please enter the HUGO Gene Symbol on the 'Genes' tab."), 
           br(),
           fluidRow(a("Feedback? Click here.", href = "https://goo.gl/forms/EoyI3da7Y0X50jah2", target="_blank"), align = "center"),
           br(),
           fluidRow(
           img(src= "sage_logo.png"),
           img(src='CTF_Logo.png'), align = "center"),
           br(),
           h4(strong("This app relies on these excellent R packages:")),
           fluidRow(
           h5("shiny, shinyBS, shinythemes, rcdk, fingerprint, rJava"),
           h5("tidyverse, DT, enrichR, webchem, visNetwork, igraph"), align = "center")), 
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
                          p("2) search CIR for structures associated with your molecule name, or"),
                          div(),
                          p("3) directly enter the structure as represented by a SMILES string."),
                          div(),
                          p(strong("After entering a molecule, pick a similarity using the slider, where 0 = highly dissimilar, and 1 = identical.")),
                          div()),
        bsCollapsePanel("Molecule Lookup",
                        fluidRow(selectizeInput("drugnames",
                                                choices = NULL,
                                                label = "",
                                                selected = "",
                                                multiple = FALSE,
                                                width = "90%"),
                                 bsTooltip("drugnames", "Type molecule name here to search this database for SMILES strings.",
                                           "right"), align = "center"),
                        fluidRow(actionButton("ppdbsearchbutton", "Find Molecules", align = "center"), align = "center"),
                        br(),
                        p("Search this database for structures by compound name. Can't find what you're looking for? Move to the next panel to search CIR."), style = "info"),
        bsCollapsePanel("CIR Search",
                 fluidRow(textInput("input.name",
                                    "input text", 
                                    label = NULL, 
                                    value = "brigatinib",
                                    width = "90%"),
                                 bsTooltip("input.name", "Type molecule name here to search this database for SMILES strings.",
                                           "right"), align = "center"),
                 fluidRow(actionButton("cirbutton", "Find CIR Mols"), align = "center"), 
                 br(),
                 p("Input a compound name in this box to search CIR for the structure."),
                 textOutput("cirsearchNA"), style = "info"),
        bsCollapsePanel("Direct Structure Input", fluidRow(textInput("smiles",
                                     "SMILES string", 
                                     label = "", 
                                     value = "",
                                     width = "90%"
                                     ),
                 bsTooltip("smiles", "Input the structural string (SMILES) here.",
                           "right"), align = "center"),
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
                           "right"), style = "warning"), open = "Similarity Threshold"),
                 uiOutput("sims"),
      bsCollapse(
      bsCollapsePanel("2D Structure (input)", fluidRow(imageOutput("structureimage"), align = "center"), style = "danger"), open = "2D Structure (input)")
      ),
    mainPanel(
      tabsetPanel(
        tabPanel(title = img("Similar Molecules  ", id = "similarmolstab", src = "help.png", align = "right"),
          bsTooltip(id = "similarmolstab", title = "This table displays all similar molecules to the input molecule.", placement = "bottom", trigger = "hover"),
          DT::dataTableOutput("simmoltab") %>% withSpinner()),
        tabPanel(title = img("Similarity Net  ", id = "similarmolnettab", src = "help.png", align = "right"),
          bsTooltip(id = "similarmolnettab", title = "This is a graphical representation of the previous table, where edge thickness is the similarity to the input.", placement = "bottom", trigger = "hover"),
          visNetworkOutput("net", height = "600px", width = "100%") %>% withSpinner(),
         DT::dataTableOutput("drug_net_info")),
        tabPanel(title = img("Targets  ", id = "targetstab", src = "help.png", align = "right"), 
          bsTooltip(id = "targetstab", title = "These are all of the targets associated with the input and similar molecules. Filter using the checkboxes in the sidebar.", placement = "bottom", trigger = "hover"),
          DT::dataTableOutput("value") %>% withSpinner()),
        tabPanel(title = img("Target Net  ", id = "targetnettab", src = "help.png", align = "right"),
                 bsTooltip(id = "targetnettab", title = "This a network of all similar molecules (blue) and their targets (green).", placement = "bottom", trigger = "hover"),
          visNetworkOutput("targetnet", height = "800px", width = "100%") %>% withSpinner()),
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
        tabPanel(title = img("in vitro Comparison  ", id = "celltab", src = "help.png", align = "right"),
                 bsTooltip(id = "celltab", title = "This tab searches CCLE and Sanger drug response data for drugs related to your query compound.", placement = "bottom", trigger = "hover"),
                 h4(textOutput("ccle_mol")),
                 div(),
                 plotlyOutput("ccle_2", height = "400px", width = "90%") %>% withSpinner(),
                 # plotlyOutput("ccle_1", height = "400px", width = "90%") %>% withSpinner(),
                 div(),
                 h4(textOutput("sang_mol")),
                 plotlyOutput("sang_2", height = "400px", width = "90%") %>% withSpinner()
                 # plotlyOutput("sang_1", height = "400px", width = "90%") %>% withSpinner()
                 ))
        
    )
  )
  ),
  tabPanel("Genes",
           sidebarLayout(
             sidebarPanel("Gene Input:",
                          fluidRow(selectizeInput(
                            'inp.gene', 
                            label = NULL,
                            choices = NULL, 
                            multiple = TRUE
                          ),
                          bsTooltip(id = "inp.gene", title = "Input 1 or more HUGO gene symbols.", placement = "right", trigger = "hover"),
                          actionButton("genebutton", "Search", align = "center"), align = "center")
                          ),
             mainPanel(
               tabsetPanel(
                 tabPanel(title = img("Matched Molecules  ", id = "matchmolstab", src = "help.png", align = "right"),
                          bsTooltip(id = "matchmolstab", title = "This tab displays molecules with a known interaction with the input gene", placement = "right", trigger = "hover"),
                          DT::dataTableOutput("genetargets") %>% withSpinner()),
                 tabPanel(title = img("Target Network  ", id = "targnetworktab", src = "help.png", align = "right"),
                            bsTooltip(id = "targnetworktab", title = "This tab graphically displays the input gene, associated drugs, and the other targets associated with those drugs.", placement = "right", trigger = "hover"),
                            visNetworkOutput("genetargetnet",  height = "800px", width = "100%") %>% withSpinner())
               )
             )
           )
           ),
  tabPanel("Settings",
           sidebarLayout(
             sidebarPanel(),
             mainPanel(
               fluidRow(h1("Molecules Tab")),
               fluidRow(
                 radioButtons("fp.type", "Fingerprint type:",
                                     c("Extended" = "extended",
                                       "Circular (ECFP/FCFP-like, default)" = "circular",
                                       # "Pubchem" = "pubchem",
                                       "MACCS" = "maccs"
                                       # ,"Klekota and Roth" = "kr"
                                       ))),
               fluidRow(p("The molecules were grouped for this database using circular fingerprints. Therefore, any other selection may result in multiple compounds having a Tanimoto similarity of 1. 
                          Also, Tanimoto similarity is greatly impacted by fingerprint choice. For example, circular fingerprints will generally have a lower Tanimoto similarity than extended fingerprints for a given chemical pair.
                          We chose extended as the default for here as it seemed to work best for commonly-used Tanimoto thresholds for similarity (e.g. >0.85 indicating highly similar.")),
               fluidRow(
                 radioButtons("edge.size", "Render edges with confidence score? (values are scaled uniquely for each plot from 0 to 10)",
                              c("Yes (default)" = TRUE,
                                "No" = FALSE))),
               fluidRow(h1("Gene Tab")),
               fluidRow(radioButtons("gene.filter.metric", "Target networks on gene panel are restricted to top 15 molecules for performance. Select on:",
                                     c("largest pChEMBL" = "mean_pchembl",
                                       "largest confidence score" = "confidence",
                                       "largest known selectivity index" = "known_selectivity_index"))))
             )
           ))
  )
  )
  )  
)

