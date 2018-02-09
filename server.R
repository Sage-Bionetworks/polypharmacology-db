source("global.R")

shinyServer(function(input, output, session) {
  
  session$sendCustomMessage(type="readCookie",
                           message=list(name="org.sagebionetworks.security.user.login.token'"))
  
# foo <- observeEvent(input$cookie, {
  # synLogin(sessionToken=input$cookie)
  synLogin()
  output$title <- renderText({
    paste0("Welcome, ", synGetUserProfile()$displayName)
  })

  source("helpers.R")
  
  loading()
  
  similarity <- reactive({
    validate(
      need(input$smiles!="", "Please select a molecule.")
    )
    validate(
      need(is.smiles(input$smiles)==TRUE, "")
    )
    similarityFunction(input$smiles, input$snappy)
  })
  
  simmols <- reactive({
    sims<-similarity()
    getSimMols(sims, input$sim.thres) %>% as.data.frame()
    })
  
  output$sims <- renderUI({
    mols <- simmols()
    choices <- mols$common_name
    checkboxGroupInput(inputId = "selectdrugs", "Molecules (similarity)", choices = choices, selected = choices)
  })
  
  output$simmoltab <- renderDataTable({
    DT::datatable(simmols(), colnames = c("Molecule Name", "Tanimoto Similarity"))
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
  

  updateSelectizeInput(session, "drugnames", choices = db.names$common_name, server = TRUE) 

  observeEvent(input$ppdbsearchbutton, {
    pp.smiles<-as.character(convertDrugToSmiles(input$drugnames)[1,1])
    updateTextInput(session, "smiles", value = pp.smiles)
  })
  
  
  ## molecule tab
  output$value <- DT::renderDataTable({
    validate(
      need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
    )
    targ <- getTargetList(input$selectdrugs) %>% as.data.frame()
    DT::datatable(targ, options = list(dom = "Bfrtip", 
                                       buttons = c("copy", 
                                                   "excel", 
                                                   "pdf", 
                                                   "print", 
                                                   "colvis")), extensions = "Buttons",
                  colnames =  c("Molecule Name", "HGNC Symbol", "Mean pChEMBL", "n Quantitative", "n Qualitative"))
  }, server = FALSE)
  
  
  output$structureimage <- renderImage({
    validate(
      need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
    )
    outfile <- tempfile(fileext='.png')
    img<-getMolImage(input$smiles)
    writePNG(img, target = outfile, dpi = 600)
    list(src = outfile,
         alt = paste("Input molecule structure:", input$smiles))
  }, deleteFile = T)


  output$net <- renderVisNetwork({
    drugsfound <- simmols()
    edges <- getNetwork(drugsfound, input$selectdrugs)
    nodes <- distinct(data.frame(id = as.character(c("input", edges$to)), 
                                 label = c("INPUT", edges$to)))
    visNetwork(nodes = nodes, edges = edges, height = "100%") %>% 
      visExport(type = "pdf", name = "exported-network", 
                float = "right", label = "Export PDF", background = "white", style= "")
    
  })
  
  output$targetnet <- renderVisNetwork({
    validate(
      need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
    )
    edges <- getTargetNetwork(input$selectdrugs)
    nodes <- distinct(data.frame(id = c(as.character(edges$from),  as.character(edges$to)), 
                                 label = c(as.character(edges$from),  as.character(edges$to)), color = c(rep("blue", length(edges$from)), 
                                                                            rep("green", length(edges$to)))))
    visNetwork(nodes = nodes, edges = edges, height = "100%") %>% visEdges(smooth = FALSE) %>% 
      visPhysics(stabilization = FALSE) %>% visLayout(randomSeed = 123) %>% 
      visIgraphLayout() %>% 
      visExport(type = "pdf", name = "exported-network", 
                float = "right", label = "Export PDF", background = "white", style= "")
  })

  
  gene.ont.mol <- reactive({
    validate(
      need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
    )
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
  
    ccleoutput <- reactive({
      validate(
        need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
      )
      plotSimCTRPDrugs(input$smiles)
    })
     
    
    output$ccle_mol <- renderText({
      validate(
        need(((nrow(as.matrix(ccleoutput()))>0)+(ncol(as.matrix(ccleoutput()))>0)==2), "No drugs found!")
      )
      
      mol <- ccleoutput() %>% top_n(1,`Tanimoto Similarity`)
      print(paste0("Nearest molecule in CCLE data is ", mol[1,1], ". Correlation calculation relative to ", mol[1,1], "."))
      
    })
    
    output$ccle_1 <- renderPlotly({
      validate(
        need(((nrow(as.matrix(ccleoutput()))>0)+(ncol(as.matrix(ccleoutput()))>0)==2),"") 
      )
    
    plot1<-ggplot(data = ccleoutput()) +
      geom_point(aes(x = Correlation, y = -`BH adj p.val`, text = cpd_name, color = (`BH adj p.val` < 0.05))) +
      theme_bw() +
      scale_color_discrete(name = "p-value < 0.05") +
      labs(x = "Drug Response Correlation", y = "BH adjusted p-value")
    
    ggplotly(p = plot1, tooltip = c("text", "x", "y"))
    
    })
    
    output$ccle_2 <- renderPlotly({
      validate(     
        need(((nrow(as.matrix(ccleoutput()))>0)+(ncol(as.matrix(ccleoutput()))>0)==2),"") 
      )
      
    plot2<-ggplot(data = ccleoutput(), aes(x = Correlation, y = `Tanimoto Similarity`, text = cpd_name, color = (`BH adj p.val` < 0.05))) +
      geom_point() +
      theme_bw() +
      scale_color_discrete(name = "p-value < 0.05") +
      labs(x = "Drug Response Correlation", y = "Chemical Similarity")
    
    ggplotly(p = plot2, tooltip = c("text", "x", "y"))
    
    })
    
    sangoutput <- reactive({
      validate(
        need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
      )
      plotSimSangDrugs(input$smiles)
      })
    
    
    output$sang_mol <- renderText({
    validate(
      need(((nrow(as.matrix(sangoutput()))>0)+(ncol(as.matrix(sangoutput()))>0)==2), "No drugs found!")
    )
    
    mol <- sangoutput() %>% top_n(1,`Tanimoto Similarity`)
    print(paste0("Nearest molecule in Sanger data is ", mol[1,1], ". Correlation calculation relative to ", mol[1,1], "."))
      
    })
    
    output$sang_1 <- renderPlotly({
      validate(
        need(((nrow(as.matrix(sangoutput()))>0)+(ncol(as.matrix(sangoutput()))>0)==2), "No drugs found!")
      )
      
      plot1<-ggplot(data = sangoutput(), aes(x = Correlation, y = -`BH adj p.val`, text = sanger_names, color = (`BH adj p.val` < 0.05))) +
        geom_point() +
        theme_bw() +
        scale_color_discrete(name = "p-value < 0.05") +
        labs(x = "Drug Response Correlation", y = "BH adjusted p-value")
      
      ggplotly(p = plot1, tooltip = c("text", "x", "y"))
      
    })
    
    output$sang_2 <- renderPlotly({
      validate(
        need(((nrow(as.matrix(sangoutput()))>0)+(ncol(as.matrix(sangoutput()))>0)==2), "")
      )
      
      plot2<-ggplot(data = sangoutput(), aes(x = Correlation, y = `Tanimoto Similarity`, text = sanger_names, color = (`BH adj p.val` < 0.05))) +
        geom_point() +
        theme_bw() +
        scale_color_discrete(name = "p-value < 0.05") +
        labs(x = "Drug Response Correlation", y = "Chemical Similarity")
      
      ggplotly(p = plot2, tooltip = c("text", "x", "y"))
      
    })
  

    
  
####### gene tab

  updateSelectizeInput(session, "inp.gene", choices = unique(db$hugo_gene), server = TRUE)
  eventReactive(input$genebutton, {  print(input$inp.gene) })

  getMols <- eventReactive(input$genebutton, {
    genes <- input$inp.gene
    validate(
      need(genes %in% db.genes, "Please enter a valid gene.")
    )
    getMolsFromGenes(input$inp.gene)
  })

  getMolNodes <- eventReactive(input$genebutton, {
    genes <- input$inp.gene
    validate(
      need(genes %in% db.genes, "Please enter a valid gene.")
    )
    getMolsFromGeneNetworks.nodes(input$inp.gene, getMols())
  })

  getMolEdges <- eventReactive(input$genebutton, {
    genes <- input$inp.gene
    validate(
      need(genes %in% db.genes, "Please enter a valid gene.")
    )
    getMolsFromGeneNetworks.edges(input$inp.gene, getMols())
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
      visLayout(randomSeed = 123) %>% visIgraphLayout() %>% 
      visExport(type = "pdf", name = "exported-network", 
                float = "right", label = "Export PDF", background = "white", style= "")
    })
  })
# })
