library(DT)
source("helpers.R")

shinyServer(function(input, output) {
  ## molecule tab
  output$value <- DT::renderDataTable({
    targ <- getTargetList(input$smiles, input$sim.thres)
    DT::datatable(targ, options = list(dom = "Bfrtip", buttons = c("copy", 
                                                                   "excel", "pdf", "print", "colvis")), extensions = "Buttons")
  }, server = FALSE)
  
  output$sims <- renderTable({
    getSimMols(input$smiles, input$sim.thres)
  })
  
  output$net <- renderVisNetwork({
    edges <- getNetwork(input$smiles, input$sim.thres)
    nodes <- distinct(data.frame(id = as.character(c("input", edges$to)), 
                                 label = c("INPUT", edges$to)))
    visNetwork(nodes = nodes, edges = edges) %>% visIgraphLayout()
  })
  
  output$targetnet <- renderVisNetwork({
    edges <- getTargetNetwork(input$smiles, input$sim.thres)
    nodes <- distinct(data.frame(id = as.character(c(edges$from, edges$to)), 
                                 label = c(edges$from, edges$to), color = c(rep("blue", length(edges$from)), 
                                                                            rep("green", length(edges$to)))))
    visNetwork(nodes = nodes, edges = edges) %>% visEdges(smooth = FALSE) %>% 
      visPhysics(stabilization = FALSE) %>% visLayout(randomSeed = 123) %>% 
      visIgraphLayout()
  })
  
  gene.ont.mol <- reactive({
    getGeneOntologyfromTargets(input$smiles, input$sim.thres)
  })
  
  output$GOMF.mol <- DT::renderDataTable({
    foo <- gene.ont.mol()[["GO_Molecular_Function_2017"]] %>% select(Term, 
                                                                     Overlap, Adjusted.P.value, Z.score) %>% filter(Adjusted.P.value < 
                                                                                                                      0.05) %>% arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy", 
                                                                  "excel", "pdf", "print", "colvis")), extensions = "Buttons")
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
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy", 
                                                                  "excel", "pdf", "print", "colvis")), extensions = "Buttons")
  }, server = FALSE)
  
  output$PPI.hub.genes <- DT::renderDataTable({
    foo <- gene.ont.mol()[["PPI_Hub_Proteins"]] %>% select(Term, 
                                                                      Overlap, Adjusted.P.value, Z.score) %>% filter(Adjusted.P.value < 
                                                                                                                       0.05) %>% arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy", 
                                                                  "excel", "pdf", "print", "colvis")), extensions = "Buttons")
  }, server = FALSE)
  
  # output$simmat <- renderImage({ getSimMolMatrix(input$smiles,
  # input$sim.thres) })
  
  ## gene tab
  output$genetargets <- renderDataTable({
    getMolsFromGenes(input$inp.gene)
  })
  
  output$genetargetnet <- renderVisNetwork({
    edges <- getMolsFromGeneNetworks.edges(input$inp.gene)
    nodes <- getMolsFromGeneNetworks.nodes(input$inp.gene)
    visNetwork(nodes = nodes, edges = edges, height = "1000px") %>% 
      visEdges(smooth = FALSE) %>% visPhysics(stabilization = FALSE) %>% 
      visLayout(randomSeed = 123) %>% visIgraphLayout()
  })
  
  output$smileslookup <- renderText({
    getSmiles(input$input.name)
  })
  
})