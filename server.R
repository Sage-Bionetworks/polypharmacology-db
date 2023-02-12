source("global.R")
source("helpers.R")

shinyServer(function(input, output, session) {

  loading()
  
  showNotification("Database last updated on March 29, 2021 (v4).", type = "message", duration = NULL)
  
  similarity <- reactive({
    validate(
      need(input$smiles!="", "Please select a molecule.")
    )
    validate(
      need(is.smiles(input$smiles)==TRUE, "")
    )
    similarityFunction(input$smiles, input$fp.type)
  })
  
  simmols <- reactive({
  	message("simmols: At start")
    sims<-similarity()
  	message("simmols: After similarity")
    result <- getSimMols(sims, input$sim.thres) %>% as.data.frame()
  	message("simmols: After getSimMols")
  	result
  })
  
  output$sims <- renderUI({
    mols <- simmols()
    choice.vals <- mols$inchikey
    choice.names <- mols$pref_name
    checkboxGroupInput(inputId = "selectdrugs", 
                       label = "Molecules (similarity)",  
                       selected = choice.vals,
                       choiceNames = choice.names,
                       choiceValues = choice.vals)
  	message("simmols: At End")
  })
  
  
  output$simmoltab <- renderDataTable({
  	message("At start of renderDataTable")
    dat <- simmols()
  	message("\trenderDataTable: After simmols")
    dat$external_links <- sapply(dat$inchikey, getExternalDrugLinks) 
  	message("\trenderDataTable: After sapply")
    dat <- select(dat, pref_name, external_links, `Tanimoto Similarity`)
  	message("\trenderDataTable: After select")
    DT::datatable(dat, colnames = c("Molecule Name", "External Links", "Tanimoto Similarity"), escape = FALSE)
  	message("At end of renderDataTable")
  })
  
  ##molecule (SMILES) input
  
  observeEvent(input$pugrestbutton, {
    pc.smiles <- convert_id_to_structure_pubchem(input$input.name, id_type = "name", output_type = "IsomericSMILES")
    updateTextInput(session, "smiles", value = pc.smiles)
    if(is.na(pc.smiles)){ 
      output$cirsearchNA <- renderText({
      "No Pubchem Molecule Found"
      })
    }
  })
  

  updateSelectizeInput(session, "drugnames", choices = db.names$synonym, server = TRUE) 

  observeEvent(input$ppdbsearchbutton, {
    pp.smiles<-convertDrugToSmiles(input$drugnames)
    updateTextInput(session, "smiles", value = pp.smiles)
  })
  
  
  ## molecule tab
  output$value <- DT::renderDataTable({
    validate(
      need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
    )
    validate(
      need(nrow(getTargetList(input$selectdrugs) %>% as.data.frame())>=1, "No targets found.")
    )
    targ <- getTargetList(input$selectdrugs) %>% as.data.frame() 
    DT::datatable(targ, options = list(dom = "Bfrtip", 
                                       buttons = c("copy", 
                                                   "excel")), extensions = "Buttons",
                  colnames =  c("InChIKey", "Molecule Name", "HGNC Symbol", "Mean pChEMBL", "n Quantitative", "n Qualitative", "KSI", "Confidence"))
  }, server = FALSE)
  
  
  output$structureimage <- renderImage({
    validate(
      need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
    )
    outfile <- tempfile(fileext='.png')
    getMolImage(input$smiles, outfile)
    list(src = outfile,
         alt = paste("Input molecule structure:", input$smiles))
  }, deleteFile = T)


  output$net <- renderVisNetwork({
    validate(
      need(nrow(getTargetList(input$selectdrugs))>=1, "No molecules found for plotting.")
    )
    drugsfound <- simmols()
    net <- getNetwork(drugsfound, input$selectdrugs) 
    edges <- net %>% select(-title)
    nodes <- distinct(data.frame(id = as.character(c("input", net$to)), 
                                 label = c("INPUT", net$to),
                                 title = c("Input Molecule", net$title)))
    visNetwork(nodes = nodes, edges = edges, height = "100%") %>% 
      visExport(type = "png", name = "exported-network", 
                float = "right", label = "Export PNG", background = "white", style= "")
      
    
  })
  
  output$targetnet <- renderVisNetwork({
    validate(
      need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
    )
    validate(
      need(nrow(getTargetList(input$selectdrugs))>=1, "No targets found for plotting.")
    )
    edges <- getTargetNetwork(input$selectdrugs, input$edge.size)

    druglinks <- sapply(edges$inchikey, function(x){
      druglinks <- getExternalDrugLinks(x)
    })
    genelinks <- sapply(edges$to, function(x){
      getExternalGeneLinks(x)
    })
    
    nodes <- distinct(data.frame(id = c(as.character(edges$from),  as.character(edges$to)), 
                                 label = c(as.character(edges$from),  as.character(edges$to)), 
                                 color = c(rep("blue", length(edges$from)), rep("green", length(edges$to))),
                                 title = c(druglinks, genelinks)
                              ))
    visNetwork(nodes = nodes, edges = edges, height = "100%") %>% visEdges(smooth = FALSE) %>% 
      visPhysics(stabilization = FALSE) %>% visLayout(randomSeed = 123) %>% 
      visIgraphLayout() %>% 
      visExport(type = "png", name = "exported-network", 
                float = "right", label = "Export PNG", background = "white", style= "") %>% 
      visOptions(highlightNearest = T)
    
  })

  
  gene.ont.mol <- reactive({
    validate(
      need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
    )
    validate(
      need(nrow(getTargetList(input$selectdrugs))>=1, "No targets found for enrichment.")
    )
    getGeneOntologyfromTargets(input$selectdrugs)
  })
  
  output$GOMF.mol <- DT::renderDataTable({
    foo <- gene.ont.mol()[["GO_Biological_Process_2018"]] %>% 
      dplyr::select(Term, Overlap, Adjusted.P.value, Combined.Score) %>%
      filter(Adjusted.P.value < 0.05) %>%
      arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy", "excel")), extensions = "Buttons")
  }, server = FALSE)
  
  output$GOCC.mol <- DT::renderDataTable({
    foo <- gene.ont.mol()[["GO_Cellular_Component_2018"]] %>% dplyr::select(Term, 
                                                                     Overlap, Adjusted.P.value, Combined.Score) %>% filter(Adjusted.P.value < 
                                                                                                                      0.05) %>% arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy", 
                                                                  "excel")), extensions = "Buttons")
  }, server = FALSE)
  
  output$GOBP.mol <- DT::renderDataTable({
    foo <- gene.ont.mol()[["GO_Biological_Process_2018"]] %>% dplyr::select(Term, 
                                                                     Overlap, Adjusted.P.value, Combined.Score) %>% filter(Adjusted.P.value < 
                                                                                                                      0.05) %>% arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy","excel")), extensions = "Buttons")
  }, server = FALSE)
  
  output$kegg <- DT::renderDataTable({
    foo <- gene.ont.mol()[["KEGG_2019_Human"]] %>% dplyr::select(Term, Overlap, Adjusted.P.value, Combined.Score) %>% 
      filter(Adjusted.P.value <0.05) %>% 
      arrange(Adjusted.P.value)
    
    DT::datatable(foo, options = list(dom = "Bfrtip", buttons = c("copy", 
                                                                  "excel")), extensions = "Buttons")
  }, server = FALSE)
  
    ccleoutput <- reactive({
      validate(
        need((is.smiles(input$smiles)==TRUE & input$smiles!=""), "Please enter a valid SMILES.")
      )
      plotSimCTRPDrugs(input$smiles, input$fp.type)
    })
     
    
    output$ccle_mol <- renderText({
      validate(
        need(is.smiles(input$smiles)==TRUE, "")
      )
      validate(
        need(((nrow(as.matrix(ccleoutput()))>0)+(ncol(as.matrix(ccleoutput()))>0)==2), "No drugs found!")
      )
      
      mol <- ccleoutput() %>% top_n(1,`Tanimoto Similarity`)
      print(paste0("Nearest molecule in CCLE data is ", mol[1,1], ". Correlation calculation relative to ", mol[1,1], "."))
      
    })
    
    # output$ccle_1 <- renderPlotly({
    #   validate(
    #     need(is.smiles(input$smiles)==TRUE, "Please enter a valid SMILES.")
    #   )
    #   
    #   validate(
    #     need(((nrow(as.matrix(ccleoutput()))>0)+(ncol(as.matrix(ccleoutput()))>0)==2),"") 
    #   )
    # 
    # plot1<-ggplot(data = ccleoutput()) +
    #   geom_point(aes(x = Correlation, y = -`BH adj p.val`, text = cpd_name, color = (`BH adj p.val` < 0.05))) +
    #   theme_bw() +
    #   scale_color_discrete(name = "p-value < 0.05") +
    #   labs(x = "Drug Response Correlation", y = "BH adjusted p-value")
    # 
    # ggplotly(p = plot1, tooltip = c("text", "x", "y"))
    # 
    # })
    # 
    output$ccle_2 <- renderPlotly({
      validate(     
        need(((nrow(as.matrix(ccleoutput()))>0)+(ncol(as.matrix(ccleoutput()))>0)==2),"") 
      )
      
    plot2<-ggplot(data = ccleoutput(), aes(x = Correlation, y = `Tanimoto Similarity`, text = cpd_name, color = (`BH adj p.val` < 0.05))) +
      geom_point() +
      theme_bw() +
      scale_color_manual(name = "p-value < 0.05", values= c("TRUE" = "#30638E", "FALSE" = "#E16036")) +
      labs(x = "Drug Response Correlation", y = "Chemical Similarity")
    
    ggplotly(p = plot2, tooltip = c("text", "x", "y"))
    
    })
    
    sangoutput <- reactive({
      validate(
        need((is.smiles(input$smiles)==TRUE & input$smiles!=""), "Please enter a valid SMILES.")
      )
      plotSimSangDrugs(input$smiles, input$fp.type)
      })
    
    
    output$sang_mol <- renderText({
    validate(
      need(((nrow(as.matrix(sangoutput()))>0)+(ncol(as.matrix(sangoutput()))>0)==2), "No drugs found!")
    )
    
    mol <- sangoutput() %>% top_n(1,`Tanimoto Similarity`)
    print(paste0("Nearest molecule in Sanger data is ", mol[1,1], ". Correlation calculation relative to ", mol[1,1], "."))
      
    })
    
    # output$sang_1 <- renderPlotly({
    #   validate(
    #     need(((nrow(as.matrix(sangoutput()))>0)+(ncol(as.matrix(sangoutput()))>0)==2), "No drugs found!")
    #   )
    #   
    #   plot1<-ggplot(data = sangoutput(), aes(x = Correlation, y = -`BH adj p.val`, text = sanger_names, color = (`BH adj p.val` < 0.05))) +
    #     geom_point() +
    #     theme_bw() +
    #     scale_color_discrete(name = "p-value < 0.05") +
    #     labs(x = "Drug Response Correlation", y = "BH adjusted p-value")
    #   
    #   ggplotly(p = plot1, tooltip = c("text", "x", "y"))
    #   
    # })
    
    output$sang_2 <- renderPlotly({
      validate(
        need(((nrow(as.matrix(sangoutput()))>0)+(ncol(as.matrix(sangoutput()))>0)==2), "")
      )
      
      plot2<-ggplot(data = sangoutput(), aes(x = Correlation, y = `Tanimoto Similarity`, text = sanger_names, color = (`BH adj p.val` < 0.05))) +
        geom_point() +
        theme_bw() +
        scale_color_manual(name = "p-value < 0.05", values= c("TRUE" = "#30638E", "FALSE" = "#E16036")) +
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
    getMolsFromGeneNetworks.nodes(input$inp.gene, getMols(), input$gene.filter.metric)
  })

  getMolEdges <- eventReactive(input$genebutton, {
    genes <- input$inp.gene
    validate(
      need(genes %in% db.genes, "Please enter a valid gene.")
    )
    getMolsFromGeneNetworks.edges(input$inp.gene, getMols(), input$edge.size, input$gene.filter.metric)
  })


  output$genetargets <- DT::renderDataTable({
    mol<-getMols()
    validate(
      need(isTRUE(nrow(mol)>0), "No molecules with this combination found.")
    )
      DT::datatable(mol, options = list(dom = "Bfrtip",
                                         buttons = c("copy","excel")),
                    colnames =  c("InChIKey", "Molecule Name", "HGNC Symbol", "Mean pChEMBL", "n Quantitative", "n Qualitative", "KSI", "Confidence"),
                    extensions = "Buttons")
  }, server = FALSE)

  output$genetargetnet <- renderVisNetwork({

    nodes <- getMolNodes()
    edges <- getMolEdges()

    validate(
      need(isTRUE(nrow(nodes)>0), "No molecules with this combination found.")
    )
    validate(
      need(isTRUE(nrow(edges)>0), "No molecules with this combination found.")
    )
    
      
    visNetwork(nodes = nodes, edges = edges, height = "100%") %>%
      visEdges(smooth = FALSE) %>% visPhysics(stabilization = FALSE) %>%
      visLayout(randomSeed = 123) %>% visIgraphLayout() %>% 
      visExport(type = "png", name = "exported-network", 
                float = "right", label = "Export PNG", background = "white", style= "") %>% 
      visOptions(highlightNearest = T)
    })
  })
#})
