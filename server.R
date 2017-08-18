library(DT)
source("helpers.R")

shinyServer(function(input, output, session) {
  
  simmols <- reactive({getSimMols(input$smiles, input$sim.thres)})
  
  output$sims <- renderUI({
    mols <- simmols()
    choices <- mols$Common_Name
    checkboxGroupInput(inputId = "selectdrugs", "Molecules (similarity)", choices = choices, selected = choices)
  })
  
  output$simmoltab <- renderDataTable({
    simmols()
  })
  
  ##molecule (SMILES) Lookput
  
  observeEvent(input$pubchembutton, {
    pc.smiles <- getSmiles(input$input.name)
    updateTextInput(session, "smiles", value = pc.smiles)
    if(is.na(pc.smiles)){ 
      output$pubchemsearchNA <- renderText({
      "No Pubchem Molecule Found"
      })
    }
  })
  
  
  updateSelectizeInput(session, "drugnames", choices = syns$Common_Name, server = TRUE) 

  observeEvent(input$ppdbsearchbutton, {
    pp.smiles<-as.character(convertDrugToSmiles(input$drugnames)[1,1])
    updateTextInput(session, "smiles", value = pp.smiles)
  })
  
  ## molecule tab
  output$value <- DT::renderDataTable({
    targ <- getTargetList(input$selectdrugs)
    DT::datatable(targ, options = list(dom = "Bfrtip", 
                                       buttons = c("copy", 
                                                   "excel", 
                                                   "pdf", 
                                                   "print", 
                                                   "colvis")), extensions = "Buttons")
  }, server = FALSE)
  

  output$net <- renderVisNetwork({
    edges <- getNetwork(input$smiles, input$sim.thres, input$selectdrugs)
    nodes <- distinct(data.frame(id = as.character(c("input", edges$to)), 
                                 label = c("INPUT", edges$to)))
    visNetwork(nodes = nodes, edges = edges) %>% visIgraphLayout()
  })
  
  output$targetnet <- renderVisNetwork({
    edges <- getTargetNetwork(input$smiles, input$sim.thres, input$selectdrugs)
    nodes <- distinct(data.frame(id = as.character(c(as.character(edges$from), edges$to)), 
                                 label = c(edges$from, edges$to), color = c(rep("blue", length(edges$from)), 
                                                                            rep("green", length(edges$to)))))
    visNetwork(nodes = nodes, edges = edges) %>% visEdges(smooth = FALSE) %>% 
      visPhysics(stabilization = FALSE) %>% visLayout(randomSeed = 123) %>% 
      visIgraphLayout()
  })
  
  gene.ont.mol <- reactive({
    getGeneOntologyfromTargets(input$selectdrugs)
  })
  
  output$GOMF.mol <- DT::renderDataTable({
    foo <- gene.ont.mol()[["GO_Molecular_Function_2017"]] %>% 
      select(Term, Overlap, Adjusted.P.value, Z.score) %>%
      filter(Adjusted.P.value < 0.05) %>%
      arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy", "excel", "pdf", "print", "colvis")), extensions = "Buttons")
  }, server = FALSE)
  
  output$GOCC.mol <- DT::renderDataTable({
    foo <- gene.ont.mol()[["GO_Cellular_Component_2017"]] %>% select(Term, 
                                                                     Overlap, Adjusted.P.value, Z.score) %>% filter(Adjusted.P.value < 
                                                                                                                      0.05) %>% arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy", 
                                                                  "excel", "pdf", "print", "colvis")), extensions = "Buttons")
  }, server = FALSE)
  
  output$GOBP.mol <- DT::renderDataTable({
    foo <- gene.ont.mol()[["GO_Biological_Process_2017"]] %>% select(Term, 
                                                                     Overlap, Adjusted.P.value, Z.score) %>% filter(Adjusted.P.value < 
                                                                                                                      0.05) %>% arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy","excel", "pdf", "print", "colvis")), extensions = "Buttons")
  }, server = FALSE)
  
  output$kegg <- DT::renderDataTable({
    foo <- gene.ont.mol()[["KEGG_2016"]] %>% select(Term, Overlap, Adjusted.P.value, Z.score) %>% 
      filter(Adjusted.P.value <0.05) %>% 
      arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy", 
                                                                  "excel", "pdf", "print", "colvis")), extensions = "Buttons")
  }, server = FALSE)
  
  ## gene tab
  output$genetargets <- renderDataTable({
    mol<-getMolsFromGenes(input$inp.gene)
    DT::datatable(mol, options = list(dom = "Bfrtip", buttons = c("copy", 
                                                                  "excel", "pdf", "print", "colvis")), extensions = "Buttons")
  }, server = FALSE)
  
  output$genetargetnet <- renderVisNetwork({
    edges <- getMolsFromGeneNetworks.edges(input$inp.gene)
    nodes <- getMolsFromGeneNetworks.nodes(input$inp.gene)
    visNetwork(nodes = nodes, edges = edges, height = "1000px") %>% 
      visEdges(smooth = FALSE) %>% visPhysics(stabilization = FALSE) %>% 
      visLayout(randomSeed = 123) %>% visIgraphLayout()
  })
})