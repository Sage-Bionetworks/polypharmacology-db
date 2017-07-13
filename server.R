source("helpers.R")

shinyServer(function(input, output){
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
  #output$simmat <- renderImage({ getSimMolMatrix(input$smiles, input$sim.thres) })
  output$genetargets <- renderDataTable({ getMolsFromGenes(input$inp.gene) })
  output$genetargetnet <- renderVisNetwork({ 
    edges<-getMolsFromGeneNetworks.edges(input$inp.gene)
    nodes<-getMolsFromGeneNetworks.nodes(input$inp.gene)
    visNetwork(nodes=nodes,edges=edges) %>%
      visEdges(smooth = FALSE) %>% 
      visPhysics(stabilization = FALSE) %>% 
      visLayout(randomSeed = 123)})
})