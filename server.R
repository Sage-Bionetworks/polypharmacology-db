loading <- function() {
  Sys.sleep(2)
  shinyjs::hide("loading_page")
  shinyjs::show("main_content")
}

source("helpers.R")
library(DT)
library(png)

shinyServer(function(input, output, session) {

  loading()
  
  simmols <- reactive({
    getSimMols(input$smiles, input$sim.thres)})
  
  output$sims <- renderUI({
    mols <- simmols()
    choices <- mols$Common_Name
    checkboxGroupInput(inputId = "selectdrugs", "Molecules (similarity)", choices = choices, selected = choices)
  })
  
  output$simmoltab <- renderDataTable({
    simmols()
  })
  
  ##molecule (SMILES) input
  
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
  
  
  output$structureimage <- renderImage({
    outfile <- tempfile(fileext='.png')
    img<-getMolImage(input$smiles)
    writePNG(img, target = outfile, dpi = 600)
    list(src = outfile,
         alt = paste("Input molecule structure:", input$smiles))
  }, deleteFile = T)


  output$net <- renderVisNetwork({
    edges <- getNetwork(input$smiles, input$sim.thres, input$selectdrugs)
    nodes <- distinct(data.frame(id = as.character(c("input", edges$to)), 
                                 label = c("INPUT", edges$to)))
    visNetwork(nodes = nodes, edges = edges, height = "2000px") #%>% visIgraphLayout()
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
  
  output$ccle <- renderPlotly({
    ccleoutput <- plotSimCTRPDrugs(input$smiles, input$sim.thres)
    
    validate(
      need(((nrow(as.matrix(ccleoutput[[1]]))>0)+(ncol(as.matrix(ccleoutput[[1]]))>0)==2), "No drugs found!") 
    )
    
    
    if(ncol(as.matrix(ccleoutput[[1]])) == 1){
      
      p <- heatmaply(as.matrix(ccleoutput[[1]]),
                     margins = c(120,100,40,20), 
                     colors = viridis(option = "C",
                                      direction = -1, 
                                      n = 256),
                     showticklabels = c(TRUE,FALSE),
                     Colv = FALSE,
                     key.title = "AUC")
    }else{
      p <- heatmaply(as.matrix(ccleoutput[[1]]),
                     margins = c(120,100,40,20),
                     colors = viridis(option = "C",
                                      direction = -1, 
                                      n = 256),
                     showticklabels = c(TRUE,FALSE),
                     key.title = "AUC")
    }                
  })
  
  output$sang <- renderPlotly({
    sangoutput <- plotSimSangDrugs(input$smiles, input$sim.thres)
    
    validate(
      need(((nrow(as.matrix(sangoutput[[1]]))>0)+(ncol(as.matrix(sangoutput[[1]]))>0)==2), "No drugs found!") 
    )
    
    
    if(ncol(as.matrix(sangoutput[[1]])) == 1){
      
      p <- heatmaply(as.matrix(sangoutput[[1]]),
                     margins = c(120,100,40,20), 
                     colors = viridis(option = "C",
                                      direction = -1, 
                                      n = 256),
                     showticklabels = c(TRUE,FALSE),
                     Colv = FALSE,
                     key.title = "AUC")
    }else{
      p <- heatmaply(as.matrix(sangoutput[[1]]),
                     margins = c(120,100,40,20),
                     colors = viridis(option = "C",
                                      direction = -1, 
                                      n = 256),
                     showticklabels = c(TRUE,FALSE),
                     key.title = "AUC")
    }                
  })
  
####### gene tab
  getMols <- eventReactive(input$genebutton, {
    validate(
      need(input$inp.gene %in% db.genes, "Please enter a valid gene.")
    )
    getMolsFromGenes(input$inp.gene)
  })
  
  getMolNodes <- eventReactive(input$genebutton, {
    validate(
      need(input$inp.gene %in% db.genes, "Please enter a valid gene.")
    )
    getMolsFromGeneNetworks.nodes(input$inp.gene)
  })
  
  getMolEdges <- eventReactive(input$genebutton, {
    validate(
      need(input$inp.gene %in% db.genes, "Please enter a valid gene.")
    )
    getMolsFromGeneNetworks.edges(input$inp.gene)
  })
  

  output$genetargets <- DT::renderDataTable({
    mol<-getMols()
      DT::datatable(mol, options = list(dom = "Bfrtip", 
                                         buttons = c("copy","excel", "pdf", "print", "colvis")), 
                    extensions = "Buttons")
  }, server = FALSE)

  output$genetargetnet <- renderVisNetwork({
    
    nodes <- getMolNodes()
    edges <- getMolEdges()
    
    visNetwork(nodes = nodes, edges = edges, height = "2000px") %>% 
      visEdges(smooth = FALSE) %>% visPhysics(stabilization = FALSE) %>% 
      visLayout(randomSeed = 123) %>% visIgraphLayout()
})


})