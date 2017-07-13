source("helpers.R")

shinyServer(function(input, output){
  ##molecule tab
  output$value <- renderDataTable({ getTargetList(input$smiles, input$sim.thres) })
  
  output$sims <- renderTable({ getSimMols(input$smiles, input$sim.thres) })
  
  output$net <- renderVisNetwork({ 
    edges<-getNetwork(input$smiles, input$sim.thres)
    nodes<-data.frame(id=as.character(c("input",edges$to)), label = c("INPUT",edges$to))
    visNetwork(nodes=nodes,edges=edges)})
  
  output$targetnet <- renderVisNetwork({ 
    edges<-getTargetNetwork(input$smiles, input$sim.thres)
    nodes<-data.frame(id=as.character(c("input",edges$to)), label = c("INPUT",edges$to))
    visNetwork(nodes=nodes,edges=edges) %>%
      visEdges(smooth = FALSE) %>% 
      visPhysics(stabilization = FALSE) %>% 
      visLayout(randomSeed = 123)})
  
  gene.ont.mol <- reactive( {getGeneOntologyfromTargets(input$smiles, input$sim.thres)} )

  output$GOMF.mol <- renderDataTable({ gene.ont.mol()[['GO_Molecular_Function_2017']]  %>% 
      select(Term, Overlap, Adjusted.P.value, Z.score) %>% 
      filter(Adjusted.P.value < 0.05) %>% 
      arrange(Adjusted.P.value) })
  
  output$GOCC.mol <- renderDataTable({ gene.ont.mol()[['GO_Cellular_Component_2017']]  %>% 
      select(Term, Overlap, Adjusted.P.value, Z.score) %>% 
      filter(Adjusted.P.value < 0.05) %>% 
      arrange(Adjusted.P.value) })
 
  output$GOBP.mol <- renderDataTable({ gene.ont.mol()[['GO_Biological_Process_2017']]  %>% 
      select(Term, Overlap, Adjusted.P.value, Z.score) %>% 
      filter(Adjusted.P.value < 0.05) %>% 
      arrange(Adjusted.P.value) })
  
  output$MSigDB.onc.mol <- renderDataTable({ gene.ont.mol()[['MSigDB_Oncogenic_Signatures']]  %>% 
      select(Term, Overlap, Adjusted.P.value, Z.score) %>% 
      filter(Adjusted.P.value < 0.05) %>% 
      arrange(Adjusted.P.value) })
  
   #output$simmat <- renderImage({ getSimMolMatrix(input$smiles, input$sim.thres) })
  
  ##gene tab
  output$genetargets <- renderDataTable({ getMolsFromGenes(input$inp.gene) })
  
  output$genetargetnet <- renderVisNetwork({ 
    edges<-getMolsFromGeneNetworks.edges(input$inp.gene)
    nodes<-getMolsFromGeneNetworks.nodes(input$inp.gene)
    visNetwork(nodes=nodes,edges=edges) %>%
      visEdges(smooth = FALSE) %>% 
      visPhysics(stabilization = FALSE) %>% 
      visLayout(randomSeed = 123)})
  
})