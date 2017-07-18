library(visNetwork)
library(shinythemes)
library(igraph)

shinyUI(navbarPage("Polypharmacology DB", theme = shinytheme("flatly"),
  tabPanel("Molecules",
  sidebarLayout(
    sidebarPanel("SMILES Input",
                  fluidRow(textInput("smiles",
                                     "input text", 
                                     label = NULL, 
                                     value = "CC(C)(C)C1=CN=C(CSC2=CN=C(NC(C3CCNCC3)=O)S2)O1",
                                     width = "100%"
                                     )), 
                 sliderInput("sim.thres", "Similarity Threshold", 
                             min=0.3, 
                             max=1,
                             value=0.90,
                             step = 0.01,
                             ticks = FALSE),
                 submitButton("Get Targets")),
    mainPanel(
      tabsetPanel(
        tabPanel(
              strong("Targets"),
              DT::dataTableOutput("value")),
        tabPanel(
              strong("Similar Molecules"),
              tableOutput("sims")),
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
              strong("PPI Hub Proteins"),
              DT::dataTableOutput("PPI.hub.genes")))),
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
                          )),
                          submitButton("Get Molecules")),
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
         ),
  tabPanel("SMILES Lookup",sidebarLayout(
    sidebarPanel("Molecule Names",
                 fluidRow(textInput("input.name",
                                    "input text", 
                                    label = NULL, 
                                    value = "brefeldin A",
                                    width = "100%"
                 )),
                 submitButton("Get Molecules")),
    mainPanel(fluidRow(textOutput("smileslookup")))))
  
))
     