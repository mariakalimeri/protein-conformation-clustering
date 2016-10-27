shinyServer(function(input, output, session) {
    source("clusterFunctionGUI.R")
    source("file.choose2.R")
    source("create.graph.from.membershipVector.R")
    source("convert2kineticmembershipVector.R")
    
    shinyjs::disable("downloadMeM")
    shinyjs::disable("downloadClus")
    
    observe({
        if (input$goButton > 0) {
            Sys.sleep(1)
            # enable the download button
            shinyjs::enable("downloadMeM")
            shinyjs::enable("downloadClus")
        }
    })

    observe({
        if (input$browsePdb == 0){ 
            return()
        } else {
            updateTextInput(session, "pdbpath",  value = file.choose2())
        }
    })
    
    observe({
        if (input$browseTrj == 0){ 
            return()
        } else {
            updateTextInput(session, "trjpath",  value = file.choose2())
        }
    })

    clusteringOutput <- eventReactive(input$goButton, {

        validate(
            need(input$trjpath != "", "Please choose a trajectory file"),
            need(input$pdbpath != "", "Please choose a PDB file"),
            need(input$dcdFrequency != "", paste("Please define the damping",
                                         "frequency of the dcd trajectory in ns")),
            need(input$chainId != "", paste("Please choose a chain identifier",
                                            "as found in the PDB file"))
        )
      
      validate(
        need(substr(input$trjpath, nchar(input$trjpath)-3,
                    nchar(input$trjpath)) == ".dcd" , 
             "Trajectory file needs to have extension '.dcd'"),
        need(substr(input$pdbpath, nchar(input$pdbpath)-3,
                    nchar(input$pdbpath)) == ".pdb" , 
             "Structure file needs to have extension '.pdb'")
      )

        chainId <- strsplit(input$chainId, ",")[[1]]
        atomTypes <- strsplit(input$atomTypes, ",")[[1]]

            withProgress(message = 'Reading input', value = 0, {
                pdb <- read.pdb(input$pdbpath)
                if (!file.exists(input$trjpath)) {
                  stop(paste("No input DCD file found with name:", input$trjpath))
                }
                trJ <- file(input$trjpath, "rb")
                nAtomsIndcd <- dcd.header(trJ)$natom
                validate(
                  need(nrow(pdb$atom)==as.numeric(nAtomsIndcd), 
                       "Trajectory .dcd and structure .pdb seem to have different
                       number of atoms")
                )
                selectedAtomsFit <- atom.select.pdb(pdb, elety=atomTypes, chain=chainId)
                selectedAtomsProt <- atom.select.pdb(pdb, string="protein", chain=chainId)  
            })
            withProgress(message = "Doing the clustering", value = 0, {
                clusteringOutput <- clusterFunctionGUI(input$trjpath, input$cutoff,
                                   fittingInds=selectedAtomsFit$xyz, 
                                   outputProtInds=selectedAtomsProt$xyz)
            })
    }) # clusteringOutput
    
    output$clusterMemberships <- renderPlot({
        totalNoOfFrames <- clusteringOutput()$nframes
        timeVector <- seq(1, totalNoOfFrames)*as.numeric(input$dcdFrequency)
        plot(timeVector, clusteringOutput()$membershipVector, type="p", 
             cex.axis=1.4, cex.lab=1.4, cex.main=1.7, col="black", 
             main="Cluster memberships", xlab="time (ns)", ylab="Cluster ID")
    }) # output$clusterMemberships
    
    output$downloadMeM <- downloadHandler(
          filename = "membershipVector.dat",
           content = function(con) {
            write.table(clusteringOutput()$membershipVector, con, row.names = F, 
                        col.names = F)
           }
    )
    
    output$noOfclusters <- renderText({ 
        max(clusteringOutput()$membershipVector) 
    }) # output$noOfclusters
    
    output$noOfclusters2 <- renderText({ 
        max(clusteringOutput()$membershipVector) 
    }) # output$noOfclusters
    
    output$vectorCluster <- renderPlot( {
        totalNoOfFrames <- clusteringOutput()$nframes
        timeVector <- seq(1, totalNoOfFrames)*as.numeric(input$dcdFrequency)
        plot(timeVector, clusteringOutput()$vectorCluster, pch=".", 
             cex.axis=1.4, cex.lab=1.4, cex.main=1.7, col="black", 
             main="Cluster growth versus time", xlab="time (ns)", 
             ylab="Number of clusters")
        try(lines(timeVector, fitting(), lty = 1, col = "blue"))
    }) # output$vectorCluster
    
    output$downloadClus <- downloadHandler(
        filename = "clusterVector.dat",
        content = function(con) {
            write.table(clusteringOutput()$vectorCluster, con, row.names = F, 
                        col.names = F)
        }
    )
    
    fitting <- eventReactive(input$fittingButton, {
        output$possibleErrorFit <- renderText({""})
        output$charactTime <- renderText({""})
        output$N_inf <- renderText({""})
        totalNoOfFrames <- clusteringOutput()$nframes
        timeVector <- seq(1, totalNoOfFrames)*as.numeric(input$dcdFrequency)
        yValues <- clusteringOutput()$vectorCluster
        mydata <- data.frame(timeVector, clusteringOutput()$vectorCluster)
        tryCatch({
            m.e <- nls(yValues ~ I(a*(1 - exp(1)^(-timeVector/b))), 
                       data = mydata, start = list(a = input$NGuess,
                                                   b = input$charTimeGuess), trace = T)
            charactTime <- round(summary(m.e)$coefficients[2, 1], 2)
            N_inf <- round(summary(m.e)$coefficients[1, 1], 2)
            output$charactTime <- renderText({ charactTime })
            output$charactTime2 <- renderText({ charactTime })
            output$N_inf <- renderText({ N_inf })
            output$N_inf2 <- renderText({ N_inf })
            m.e <- predict(m.e, list(x = timeVector))
        },
        error=function(cond) {
            output$possibleErrorFit <- renderText({"Fitting is not possible. 
                If you are confident it should be, try different initial guesses"})
            return(NA)
        })
    }) # fitting
    
    # ===================== STEP 2 SERVER ======================================
    
    output$trjfileName2Display <- renderText({tail(strsplit(input$trjpath, 
                                                            "/")[[1]], 1)})
    
    output$pdbfileName2Display <- renderText({tail(strsplit(input$pdbpath, 
                                                            "/")[[1]], 1)})
    
    output$clusterAtoms <- renderText({ input$atomTypes })
    
    output$clusteringCutoff <- renderText({ input$cutoff })
    
    writePDBfunc <- function(coN){
        noOfClus <- max(clusteringOutput()$membershipVector)
        pdb <- read.pdb(input$pdbpath)
        selectedAtomsProt <- atom.select.pdb(pdb, string="protein",
                                             chain=input$chainId)
        xyzCoords <- clusteringOutput()$xyzLeaders
        write.pdb(NULL, file=coN, xyz = xyzCoords[1,], 
                  type=pdb$atom[[1]][selectedAtomsProt$atom], 
                  resid=pdb$atom[[5]][selectedAtomsProt$atom], 
                  resno=pdb$atom[[7]][selectedAtomsProt$atom], 
                  eleno=pdb$atom[[2]][selectedAtomsProt$atom], 
                  elety=pdb$atom[[3]][selectedAtomsProt$atom], 
                  chain=pdb$atom[[6]][selectedAtomsProt$atom])
        for (i in 2:noOfClus){
            write.pdb(NULL, file=coN, xyz = xyzCoords[i,], 
                      type=pdb$atom[[1]][selectedAtomsProt$atom], 
                      resid=pdb$atom[[5]][selectedAtomsProt$atom], 
                      resno=pdb$atom[[7]][selectedAtomsProt$atom], 
                      eleno=pdb$atom[[2]][selectedAtomsProt$atom], 
                      elety=pdb$atom[[3]][selectedAtomsProt$atom], 
                      chain=pdb$atom[[6]][selectedAtomsProt$atom], append = T)
        }
    }
    
    output$downloadPDB <- downloadHandler(
        filename = "leaders.pdb",
        content = function(con) {
            withProgress(message = 'Writing PDB file', value = 0,{
                writePDBfunc(con)
            })
        }
    )
    
    observeEvent(input$renderButton, {
        
        withProgress(message = 'Rendering images', value = 0, {
            filename = "leadersTemp.pdb"
            writePDBfunc(filename)
            freqStride <- input$stridingForRendering
            system(paste("./render.sh", freqStride))
            output$renderingComplete <- renderText({"Rendering completed!"})
            output$renderedImage1 <- renderImage({
                list(src = "www/renderedImage1.tga",
                     contentType = 'image/png',
                     width = 350,
                     alt = "")
            }, deleteFile = F)
            output$renderedImage2 <- renderImage({
                list(src = "www/renderedImage2.tga",
                     contentType = 'image/png',
                     width = 350,
                     alt = "")
            }, deleteFile = F)
            output$renderedImage3 <- renderImage({
                list(src = "www/renderedImage3.tga",
                     contentType = 'image/png',
                     width = 350,
                     alt = "")
            }, deleteFile = F)
            system("rm leadersTemp.pdb")
        })
    })
    
    output$downloadImage1 <- downloadHandler(
        filename = c("image1.tga"),
        content = function(con) {
            system(paste("cp www/renderedImage1.tga", con, sep=" "))
        }
    )
    
    output$downloadImage2 <- downloadHandler(
        filename = c("image2.tga"),
        content = function(con) {
            system(paste("cp www/renderedImage2.tga", con, sep=" "))
        }
    )
    
    output$downloadImage3 <- downloadHandler(
        filename = c("image3.tga"),
        content = function(con) {
            system(paste("cp www/renderedImage3.tga", con, sep=" "))
        }
    )
    
    # ===================== STEP 3 SERVER ======================================
    
    # WORKAROUND I had to do the actual plotting outside the eventReactive
    # because the downloadHandler below does not save reactive objects!
    myG <- eventReactive(input$drawGraph, { 
        validate(
            need(input$labelColor != "", "Choose a color for labels"),
            need(input$nodeColor != "", "Choose a color for the nodes"),
            need(input$minNodeSize != "", "Choose minimum node size"),
            need(input$maxNodeSize != "", "Choose maximum node size")
        )
        # For debugging of Step 3
        # fileName <- "sampleFiles_dcd_n_pdb/membershipVector.dat"
        # membershipVectorFromFile <- as.matrix(read.table(fileName))
        # enrichedGraph <- create.graph.from.membershipVector(membershipVectorFromFile)
        enrichedGraph <- 
            create.graph.from.membershipVector(clusteringOutput()$membershipVector)
        myG <- enrichedGraph$myG
        weightedEdgeList <- enrichedGraph$weightedEdgeList
        if (input$drawingAlg=="layout_with_graphopt"){
            co <- eval(parse(text=input$drawingAlg))(myG, 
                                spring.constant=as.numeric(input$springConst))
        } else {
            co <- eval(parse(text=input$drawingAlg))(myG)
        }
        return(list(graph=myG, layout=co, graphAsWeightedEdgeList=weightedEdgeList))
    })

    myGplotted <- function(){
        nodeWeights <- get.vertex.attribute(myG()$graph, "weight") 
        desiredInterval <- c(as.numeric(input$minNodeSize), as.numeric(input$maxNodeSize))
        normalizedNodeSizes <- (desiredInterval[2]-desiredInterval[1])*
                   ((nodeWeights-min(nodeWeights))/(max(nodeWeights)-min(nodeWeights)))+
                                                        desiredInterval[1]
        netPlot <- plot.igraph(myG()$graph, layout=myG()$layout, 
                               vertex.color=input$nodeColor, 
                               vertex.size=normalizedNodeSizes,
                               edge.color="black", edge.arrow.size=0.5, edge.width=1.5,
                               vertex.label.cex=input$labelSize, 
                               vertex.label.color=input$labelColor)
        return(netPlot)
    }
        
    output$networkOfStates <- renderPlot({myGplotted()})
    
    output$downloadNet <- downloadHandler(
      filename = "network.pdf",
      content = function(con) {
        pdf(con)
        myGplotted()
        dev.off()
      })
    
    myKineticG <- eventReactive(input$drawKineticGraph, {
        validate(
            need(input$labelKinColor != "", "Choose a color for labels"),
            need(input$nodeKinColor != "", "Choose a color for the nodes"),
            need(input$minNodeKinSize != "", "Choose minimum node size"),
            need(input$maxNodeKinSize != "", "Choose maximum node size")
        )
        # Write temp file
        write.table(myG()$graphAsWeightedEdgeList, file="mclInput.txt", col.names = F, 
                    row.names = F)
        system(paste("mcl mclInput.txt --abc -I", input$inflation, " -o mclOut.txt"))
        kinMembershipVector <-
            convert2kineticMembershipVector(clusteringOutput()$membershipVector)
        # For debugging of Step 3
        # fileName <- "sampleFiles_dcd_n_pdb/membershipVector.dat"
        # membershipVectorFromFile <- as.matrix(read.table(fileName))
        # kinMembershipVector <-
        #     convert2kineticMembershipVector(membershipVectorFromFile)
        system("rm mclInput.txt mclOut.txt")
        tryCatch({
            myKinG <- create.graph.from.membershipVector(kinMembershipVector)
            myKinG <- myKinG$myG
            coKin <- eval(parse(text=input$drawingKinAlg))(myKinG)
        },
        error=function(cond) {
            stop({"Coarse-graining was not possible. Try increasing the inflation 
                parameter. Otherwise, your initial network may be already too small."})
                return(NA)
        })
        return(list(graph=myKinG, layout=coKin))
    })
    
    myKinGplotted <- function(){
        nodeKinWeights <- get.vertex.attribute(myKineticG()$graph, "weight")
        desiredKinInterval <- c(as.numeric(input$minNodeKinSize),
                             as.numeric(input$maxNodeKinSize))
        normalizedKinNodeSizes <- (desiredKinInterval[2]-desiredKinInterval[1])*
            ((nodeKinWeights-min(nodeKinWeights))/(max(nodeKinWeights)-min(nodeKinWeights)))+
            desiredKinInterval[1]
        myKinPlot <- plot.igraph(myKineticG()$graph,
                               layout=myKineticG()$layout, vertex.color=input$nodeKinColor,
                               vertex.size=normalizedKinNodeSizes,
                               edge.color="black", edge.arrow.size=0.5, edge.width=1.5,
                               vertex.label.cex=input$labelKinSize,
                               vertex.label.color=input$labelKinColor)
        return(myKinPlot)
    }

    output$kineticnetworkOfStates <- renderPlot({myKinGplotted()})

    output$downloadKineticNet <- downloadHandler(
        filename = "kineticallyClusteredNetwork.pdf",
        content = function(con) {
            pdf(con)
            myKinGplotted()
            dev.off()
    })
})