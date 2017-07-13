library(visNetwork)

shinyUI(navbarPage(strong("Polypharmacology DB"),
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
                             min=0.3, max=1, value=0.90),
                 submitButton("Get Targets")),
    mainPanel(
      tabsetPanel(
        tabPanel(
              strong("Targets"),
              dataTableOutput("value")),
        tabPanel(
              strong("Similar Molecules"),
              tableOutput("sims")),
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
           
           )
))