shinyUI(fluidPage(

  # Use shinyjs package to disable download buttons before their due time
  shinyjs::useShinyjs(),

  # Make error-validation messages red
  tags$head(
    tags$style(HTML("
      .shiny-output-error-validation {
        color: red;
        font-family: 'Monaco', cursive;
        font-size: 12px;
      }
    "))
  ),

  titlePanel(
    wellPanel(
      span("Protein Conformational Clustering",
        style = "font-family: 'Helvetica'; font-si16pt"
      )
    ),
    windowTitle = "Protein Conformational Clustering"
  ),

  tabsetPanel(
    # =================== TAB STEP 1 ========================================
    tabPanel(
      p(strong("Step 1:"), "Leader clustering"),
      fluidRow(
        column(
          9,
          wellPanel(
            helpText(
              p("Welcome to protein conformational clustering!"),
              p("This is the first step in a methodology used to 
                               probe protein dynamics and flexibility. Please 
                               visit the 'Read me' tab for an overview of this app.
                               After that, you may click below to read more about
                               the first step.")
            ),
            actionLink("readMoreTab1Link", "Read more about this step")
          )
        )
      ),
      fluidRow(
        # =============== INPUT PANEL ON THE LEFT ==============================
        column(
          4,
          wellPanel(
            h4("Input parameters"),
            wellPanel(
              style = "background-color: #E3E3E3;",
              textInput("trjpath", "Choose a trajectory file (.dcd)"),
              actionButton("browseTrj", "Browse",
                style = "padding:4px; font-size:80%"
              )
            ),
            wellPanel(
              style = "background-color: #E3E3E3;",
              textInput("pdbpath", "Choose a structure file (.pdb)"),
              actionButton("browsePdb", "Browse",
                style = "padding:4px; font-size:80%"
              )
            ),
            textInput(
              inputId = "dcdFrequency", label = "DCD frequency (ns)",
              placeholder = "e.g. 0.1"
            ),
            textInput(
              inputId = "chainId", label = "PDB chain identifier/s",
              placeholder = "e.g. A or A,B for dimers e.t.c. "
            ),
            selectInput(
              inputId = "atomTypes", label = "Atom set for clustering",
              choices = c("Calpha" = "CA", "Backbone" = "C,O,N,CA")
            ),
            sliderInput(
              inputId = "cutoff", label = "RMSD cutoff for clustering (Å)",
              min = 1, max = 5, value = 2.0, step = 0.1
            ),
            actionButton("goButton", "Run clustering")
          )
        ),

        # =============== RESULT PANEL ON THE RIGHT ============================
        column(
          8,
          # ================ UPPER PLOT ===================================
          wellPanel(
            style = "background-color: #ffffff;",
            h4("Output"),
            h5("Number of clusters versus time"),
            plotOutput("vectorCluster"),
            wellPanel(
              strong("Fitting on cluster growth"),
              withMathJax(p("In principle, the number of clusters would
                                     saturate exponentially fast to a limiting 
                                     value according to the model 
                                     $$N = N_\\infty (1-e^{t/\\tau})$$")),
              fluidRow(
                column(
                  6,
                  textInput(
                    inputId = "NGuess", label = "",
                    placeholder = "Initial guess for N_inf "
                  ),
                  span("Number of clusters found:"),
                  span(textOutput("noOfclusters"),
                    style = "color:blue"
                  ),
                  span("Number of clusters at the limit: "),
                  span(textOutput("N_inf"), style = "color:blue"),
                  span("Characteristic time in ns: "),
                  span(textOutput("charactTime"),
                    style = "color:blue"
                  )
                ),
                column(6,
                  align = "right",
                  textInput(
                    inputId = "charTimeGuess", label = "",
                    placeholder = "Initial guess for tau (ns)"
                  ),
                  actionButton("fittingButton", "Try fitting"),
                  code(textOutput("possibleErrorFit")),
                  br(),
                  downloadButton("downloadClus", "Download data")
                )
              )
            )
          ), # wellPanel vectorCluster plot

          # ================ LOWER PLOT  ==================================
          wellPanel(
            style = "background-color: #ffffff;",
            h5("Population of clusters over time"),
            fluidRow(
              plotOutput("clusterMemberships"),
              column(12,
                align = "right",
                downloadButton("downloadMeM", "Download data")
              )
            )
          ) # wellPanel for membershipVector plot
        ) # column result panel on the right
      ), # fluidrow
      bsModal("tab1Expain",
        title = "Leader clustering of MD protein trajectories",
        trigger = "readMoreTab1Link",
        h4("Input"),
        tags$ul(
          tags$li(
            strong("Files:"), "Provide a molecular 
              dynamics trajectory file in a DCD binary format alongside with its 
              corresponding PDB file. 
              The files must include at least one protein on which the clustering 
              will be performed and they may be a raw output of a simulation,
              including water, ions etc; the protein for clustering will be selected,
              automatically, given the input PDB chain identifier (therefore the 
              species that is intented for clustering must have its own chain ID 
              in the PDB).",
            br(),
            span(strong("Note 1:"), "You do not need to worry about memory overheads
              concerning the size of the DCD trajectory: the clustering is done 
              in a read-frame-and-trash manner."),
            br(),
            span(strong("Note 2:"), "Support for xtc trajectory files will be  
              available in the upcoming next version.")
          ),
          tags$li(strong("DCD frequency (ns):"), "The damping frequency of 
              the trajectory file is needed for the fitting and the units of the 
              x-axis in the output plots."),
          tags$li(strong("RMSD cutoff:"), "This may be the hardest thing to predict.
                 In general, the smaller the cutoff, the larger the number of the 
                 resulting clusters and the less likely is to be able to fit the
                 cluster growth with the model shown below. Also the more clusters
                    you get the slower the program will be. A good initial
                    guess may be, minus one the RMSD distance between the aligned first 
                    and last frames of the trajectory.")
        ),
        h4("Output"),
        span("You will get two plots and download buttons that give you the
                 possibility to download the data. Assuming you have a good 
                 understanding of what the output means, the top plot is the cluster growth
                 versus time and you may try to fit it with the provided model. 
                 The second plot contains the same information as in the first,
                 with the additional occupation of each cluster across time:
                 a cluster k was revisited as many times as there are point along 
                 the line y=k.")
      )
    ),
    # ========================== TAB STEP 2 ====================================
    tabPanel(
      p(strong("Step 2:"), "Visualizing the leaders"),
      fluidRow(
        column(9,
          offset = 1,
          wellPanel(
            helpText(
              p("As a second step, you need to visualize 
                                  the most prevalent conformations of the protein. 
                                  This will also help to spot where, exactly, 
                                  the flexibility is localized on the protein matrix."),
              p("To that end, we will focus on the leaders 
                                  (a.k.a. centroids) as found by 
                                  the previous clustering. Those are, hopefully, well
                                  spread out across the accessible conformational space.")
            ),
            actionLink("readMoreTab2Link", "Read more about 
                                           this step")
          )
        )
      ),
      fluidRow(
        column(
          4,
          wellPanel(
            h4("Clustering review"),
            span("Analysed trajectory:"),
            # span("and the PDB file:"),
            span(textOutput("pdbfileName2Display"), style = "color:blue"),
            span(textOutput("trjfileName2Display"), style = "color:blue"),
            span("Atoms used for RMSD:"),
            span(textOutput("clusterAtoms"), style = "color:blue"),
            span("with cutoff distance:"),
            span(textOutput("clusteringCutoff"), style = "color:blue"),

            wellPanel(
              style = "background-color: #E3E3E3;",

              span("Number of clusters found:"),
              # Note for me: Unfortunately shiny doesn't support the
              # same output twice. So if I want to display the follo-
              # wing stuff in tab2 as well, I have to redefine them...
              span(textOutput("noOfclusters2"), style = "color:blue"),
              span("Characteristic time in ns: "),
              span(textOutput("charactTime2"), style = "color:blue"),
              span("Limiting number of clusters: "),
              span(textOutput("N_inf2"), style = "color:blue")
            ),
            downloadButton("downloadPDB", "Download leaders PDB"),
            helpText("(This may take a few seconds)")
          )
        ),
        column(8,
          style = "background-color: #ffffff;",
          h4("Molecular representation of protein conformations"),
          fluidRow(
            column(
              6,
              textInput("stridingForRendering",
                "Choose striding frequency for rendering",
                value = "1", placeholder = "e.g. 2, 5, ..."
              )
            ),
            column(6,
              align = "left",
              br(),
              br(),
              actionLink("stridingExplanation", "What is this?")
            )
          ),
          actionButton("renderButton", "Render"),
          helpText(p("(This may take a few seconds)")),
          span(textOutput("renderingComplete"), style = "color:blue"),
          conditionalPanel(
            "output.renderingComplete =='Rendering completed!'",
            br(),
            column(
              6,
              imageOutput("renderedImage1"),
              downloadButton("downloadImage1", "Save rendering")
            ),
            column(
              6,
              imageOutput("renderedImage2"),
              downloadButton("downloadImage2", "Save rendering")
            ),
            column(12,
              offset = 3,
              imageOutput("renderedImage3"),
              downloadButton("downloadImage3", "Save rendering")
            )
          )
        )
      ),
      bsModal("tab2Expain",
        title = "Molecular visualization",
        trigger = "readMoreTab2Link",
        p("No extra input should be provided here. On the left column
                  you may review the clustering output from the previous tab.
                  At the end of that column, you can download the pdb file that 
                  contains the alined centroid configurations. You can use this file 
                  in VMD or any other molecular visualization software."),
        p(strong("The right panel is not available on the server version of this
                  app.")),
        p("If you are using this app locally and, additionally, you have VMD
                  installed in your machine, you may render the molecular representations
                  by pressing the render button on the right panel. The rendering
                  may take 1 to 2 minutes. You will get the centroid conformations
                  aligned on top of each other. The 3 output images are the same 
                  represenation viewed from 3 different angles.")
      ),
      bsModal("stridingExplanationModal",
        title = "Striding frequency for visualization",
        trigger = "stridingExplanation",
        p("If you found many clusters (i.e. chose a fine cutoff for your system)
                  that is not necessarily bad for step 3, yet it doesn't fascilitate
                  visualization because leaders will be displayed on top of each other."),
        p("Ideally, 3 to 5 overlapping conformations are easy to interpret. 
                  So depending on the total number of clusters
                  you found, choose a striding frequency that would give you no
                  more than 5 leading confomrations.
                  E.g. for 37 number of clusters, a frequency of 9 is good")
      )
    ),
    # =========================== TAB STEP 3 ===================================
    tabPanel(
      p(strong("Step 3:"), "Network of states"),
      fluidRow(
        column(9,
          offset = 2,
          wellPanel(
            helpText(
              p("You will now use the clustering of the first step
                                   to map the protein dynamics on a network of states."),
              p(" That is, you will create a graph, whose nodes 
                                   represent the centroid conformations of the 
                                   clustering in step one.")
            ),
            actionLink("readMoreTab3Link", "Read more about 
                                           this step")
          )
        )
      ),
      fluidRow(
        column(
          5,
          wellPanel(
            h4("Input parameters"),
            selectInput("drawingAlg", "Select graph drawing algorithm",
              choices = c(
                "graphOpt" = "layout_with_graphopt",
                "Fruchterman and Reingold" = "layout_with_fr",
                "GEM" = "layout_with_gem",
                "DrL" = "layout_with_drl"
              )
            ),
            conditionalPanel(
              'input.drawingAlg=="layout_with_graphopt"',
              sliderInput(
                inputId = "springConst", label = "Spring constant",
                min = 0.1, max = 2, value = 0.6, step = 0.2
              )
            ),
            column(
              6,
              textInput(
                inputId = "nodeColor", label = "Node color",
                value = "red", placeholder = "e.g. red"
              ),
              textInput(
                inputId = "minNodeSize", label = "Smallest node",
                value = "20", placeholder = "e.g. 20"
              )
            ),
            column(
              6,
              textInput(
                inputId = "labelColor", label = "Label color",
                value = "white", placeholder = "e.g. white"
              ),
              textInput(
                inputId = "maxNodeSize", label = "Largest node",
                value = "50", placeholder = "e.g. 50"
              )
            ),

            sliderInput(
              inputId = "labelSize", label = "Label size",
              min = 0.5, max = 3, value = 1.5, step = 0.2
            ),
            actionButton("drawGraph", "Draw graph")
          )
        ),

        # =============== RESULT PANEL ON THE RIGHT ============================
        column(
          7,
          wellPanel(
            style = "background-color: #ffffff;",
            h4("Network of protein states"),
            br(),
            plotOutput("networkOfStates", width = 500),
            fluidRow(
              column(12,
                align = "right",
                downloadButton("downloadNet", "Save image")
              )
            )
          )
        )
      ),
      tags$hr(),
      fluidRow(
        column(
          5,
          wellPanel(
            h4("Input parameters"),
            sliderInput("inflation", "Inflation parameter for MCL",
              min = 1.2, max = 5, value = 2, step = 0.1
            ),
            selectInput("drawingKinAlg", "Select graph drawing algorithm",
              choices = c(
                "graphOpt" = "layout_with_graphopt",
                "Fruchterman and Reingold" = "layout_with_fr",
                "GEM" = "layout_with_gem",
                "DrL" = "layout_with_drl"
              )
            ),
            column(
              6,
              textInput(
                inputId = "nodeKinColor", label = "Node color",
                value = "blue", placeholder = "e.g. red"
              ),
              textInput(
                inputId = "minNodeKinSize", label = "Smallest node",
                value = "20", placeholder = "e.g. 20"
              )
            ),
            column(
              6,
              textInput(
                inputId = "labelKinColor", label = "Label color",
                value = "white", placeholder = "e.g. white"
              ),
              textInput(
                inputId = "maxNodeKinSize", label = "Largest node",
                value = "50", placeholder = "e.g. 50"
              )
            ),

            sliderInput(
              inputId = "labelKinSize", label = "Label size",
              min = 0.5, max = 3, value = 1.5, step = 0.2
            ),
            actionButton("drawKineticGraph", "Draw coarse-grained graph")
          )
        ),

        # =============== RESULT PANEL ON THE RIGHT ============================
        column(
          7,
          wellPanel(
            style = "background-color: #ffffff;",
            h4("Kinetically coarse-grained network"),
            br(),
            plotOutput("kineticnetworkOfStates", width = 500),
            fluidRow(
              column(12,
                align = "right",
                downloadButton("downloadKineticNet", "Save image")
              )
            )
          )
        )
      ),
      bsModal("tab3Expain",
        title = "Network of states",
        trigger = "readMoreTab3Link",
        p(strong("Top Panel:"), "To begin with you don't need to modify the 
              default input in the left column. Simply press the 'Draw graph' 
              button. With the graph plotted you may adjust the size and color 
              of the nodes and labels, since this will depend on your system, 
              the total number of nodes e.t.c."),
        p(
          strong("Graph drawing with force based algorithms:"), "At the very 
              beginning of the input column, there is a small colection of force 
              based algorithms to chose from and each of the will draw your network 
              slightly differently ",
          tags$a(
            href = "http://igraph.org/c/doc/igraph-Layout.html",
            "(see igraph manual)", target = "_blank"
          ), ". For more info on 
              force-directed algorithms see ",
          tags$a(
            href = "https://en.wikipedia.org/wiki/Force-directed_graph_drawing",
            "wikipedia.", target = "_blank"
          )
        ),
        p(strong("Lower panel:"), "The lower panel will take your analysis 
              one step further. It will coarse-grain your existing network
              merging the states into kinetically separated substates. It uses 
              and algorithm based on Markov Modeling and the software can be 
              downloaded from", tags$a(
          href = "http://micans.org/mcl/", "here",
          target = "_blank"
        ), "(Stijn van Dongen, 2000)")
      )
    ),
    # ====================== TAB README ========================================
    tabPanel(
      "Read me!",
      column(
        8,
        wellPanel(
          fluidRow(
            column(
              12,
              h3("Introduction"),
              p("This is a user friendly edition of a pipeline of scripts 
                that were assembled to probe protein dynamics and flexibility. 
                The methodology applied here is an analysis tool of molecular 
                dynamics (MD) simulations of proteins.")
            )
          )
        ),
        wellPanel(
          fluidRow(
            column(
              12,
              h3("The methodology in a nutshell"),
              p("Cluster analysis is the process of grouping together 
                          multidimensional data based on a similarity measure.
                          It can be very useful for identifying the pricipal 
                          conformations of a protein in an MD simulation."),
              p(
                actionLink("leaderLink", "The leader algorithm"),
                " (or leader/follower)
                          used here is a partitional clustering algorithm. It is
                          particularly useful for MD because it is memory- and time-
                          efficient as compared to both of its equally popular 
                          candidates,", actionLink("kmeansLink", "k-means"), " and",
                actionLink("gromosLink", "gromos"), ". The leader algorithm
                          takes as an input the cutoff distance between the cluster
                          centroids. So as opposed to, for example k-means, the number 
                          of clusters is one of the outputs of the algorithm."
              ),
              p(
                actionLink("pitFall", "With some care"),
                "the output number of clusters, can be 
                          used as a quantification of a molecule's flexibility, 
                            at least in the timescale explored by the specific 
                            simulation (see step 1). The protein conformations that
                            correspond to the cluster leaders can be visualized to 
                            inspect where the flexibility is localized on the matrix 
                            of a molecule (see step 2). What is more, the initial MD
                            trajectory of the protein can be mapped onto a network of
                            states that effectively projects the topology of the 
                            underlying conformational landscape on 2 dimensions 
                            (see step 3)."
              ),
              p("For more information on the methodology and the possible 
                        pitfalls have a look on pages 42-46 of this", tags$a(
                href = "https://drive.google.com/file/d/0B12DgiTeqEVfTUxMc05QbjRZaGc/view", "thesis.", target = "_blank"
              ))
            )
          )
        ),
        wellPanel(
          fluidRow(
            column(
              12,
              h3("Using this app"),
              p("There are 3 steps in the proposed methodology:"),
              tags$ul(
                tags$li(strong("Step 1:"), "The actual protein clustering that 
                          will give you as output the number of clusters and 
                          the way they increase over the course of the trajectory"),
                tags$li(strong("Step 2:"), "Visualization of  the 
                          principal conformations of the protein, 
                          i.e. the conformations that correspond to the 
                          leaders of the clusters. "),
                tags$li(strong("Step 3:"), "Mapping the protein dynamics 
                          on a network of states.")
              ),
              p("If you feel confortable enough with the idea of clustering
                        (as well as with basic network concepts if you plan to use 
                        step 3) go ahead and visit Tab 1. In each step, there is 
                        a 'Read more about this step' link that should give you 
                        useful additional info."),
              br(),
              h4("Output assesment"),
              p("If you are new to using clustering in your MD analysis,
                        here are a few tips that can help you assess the 
                        goodness of your results."),
              tags$ul(
                tags$li("In the first step, the second plot shows you how 
                          the clusters are populated over time: a cluster n has 
                          been revisited many times if the line y=n has more than
                          one plot points. Cluster occupancy, in general, can be used 
                          as a good measure of 'equilibration' for your system: 
                          The protein is well equilibrated when it has started 
                          revisiting back and forth previous clusters, i.e. as 
                          we say is ergodic. Also, ideally a saturation in both 
                          plots should be visible. If you don't see a saturation,
                           try increasing the cutoff."),
                tags$li("When comparing different systems, the trend, as far as 
                          their maximum number of clusters is concerned, tends to 
                          be robust against different cutoffs. However, it is 
                         always a good tactic to verify this for yourself. Remember 
                           to compare only systems that are homologous, with more 
                           or less the same number of a.a., the same cutoff and 
                           the same trajectory length and dumping frequency."),
                tags$li(
                  strong("For gromacs users:"), "You may use the ",
                  tags$a(
                    href = "http://mdtraj.org/latest/mdconvert.html",
                    em("mdconvert"), target = "_blank"
                  ),
                  "to convert your trajectories from
                           .xtc to .dcd. However, remember to unwrap your proteins
                           beforehand (specially multidomain proteins must not
                           be broken across the box!)"
                )
              ),
              p("The list of tips is far from complete. Contact me for 
                        any questions/suggestions (mkalime@gmail.com)")
            )
          )
        )
      ),
      column(
        4,
        wellPanel(
          h4("Version 1.0b1", h5("(Sep. 2016)")),
          p("Feedback and bug report at mkalime@gmail.com")
        ),
        wellPanel(
          align = "center",
          style = "background-color: #ffffff;",
          h4("References", align = "left"),
          img(src = "jpcb.png", width = 220),
          tags$ul(
            br(),
            tags$li(
              "Methodology on pages 42-46 of this ",
              tags$a(
                href = "https://drive.google.com/file/d/0B12DgiTeqEVfTUxMc05QbjRZaGc/view", 
                "pdf document", 
                target = "_blank"
              )
            ),
            tags$li("Katava M, Kalimeri M, Stirnemann G and Sterpone F, 
                    JPCB (2016)"),
            tags$li("Kalimeri M, Derreumaux P and Sterpone F, 
                    JNCS (2015)"),
            tags$li("Kalimeri M, Girard E, Madern D and Sterpone F, 
                    PLoS One (2015)"),
            tags$li("Kalimeri M, Rahaman O, Melchionna S and Sterpone F, 
                    JPCB (2015)")
          ),
          br(),
          h5("External software used:"), align = "left",
          tags$ul(
            align = "left",
            tags$li(tags$a(
              href = "http://thegrantlab.org/bio3d/",
              "Bio3d, R library", target = "_blank"
            )),
            tags$li(tags$a(
              href = "http://igraph.org/redirect.html",
              "igraph, R library", target = "_blank"
            )),
            tags$li(tags$a(
              href = "http://www.ks.uiuc.edu/Research/vmd/",
              "Visual Molecular Dynamics", target = "_blank"
            )),
            tags$li(tags$a(
              href = "http://micans.org/mcl/",
              "Markov Clustering Algorithm", target = "_blank"
            ))
          )
          # tags$link(href="http://www.lalala.com"),
          #     div(class="header", checked=NA,
          #         p("Ready to take the Shiny tutorial? If so"),
          #         a(href="shiny.rstudio.com/tutorial", "Click Here!")
          #     )
          #
        )
      ),

      # ============ AUXILIARY MODALS ============================================
      bsModal("pitFallFlex",
        title = "Using conformational clustering to quantify flexibility",
        trigger = "pitFall",
        p("Roughly speaking, the folded state of proteins is kept intact as 
              a result of cohesive forces between the constituent amino acids, 
              i.e. van de Waals interactions, electrostatics e.t.c. Thus, folded 
              proteins are soft-matter entities, some more flexible than others. 
              In fact, defining protein flexibility is a hard task. One of the 
              reasons is that their motions span a huge range of different 
              timescales (from thermal fluctuations, in the order of femptoseconds
              or picoseconds, to allosteric movements and folding/unfolding events,
              in the order of microseconds to hours). 
              Another reason is that flexibility may be localized in a particular 
              region of the protein matrix or may span the  
              the whole protein matrix (see intrinsically disordered proteins)."),
        p("Here, we attempt a computational approach in order to quantify protein
              flexibility within the timescales explored by the MD simulation: 
              Given a certain cutoff, the number of clusters found by the leader 
              algorithm quantify the avaliable conformational space of the protein."),
        p(
          strong("There are several possible pitfalls when comparing proteins
             with each other:"),
          tags$li("Since the clustering depends on the chosen cutoff, proteins
              should always be compared for the same cutoff!"),
          tags$li("The proteins under comparison must have the same, or more 
              or less the same, amino acid length. (In principle, for a random 
              polymer the possible number of clusters grows exponentially with 
              the length of the chain. In practice, for folded proteins this 
              dependency is much weaker.)"),
          tags$li("The protein conformational space is multidimensional. Therefore,
              it suffers from the curse of multidimentionality, i.e. it is very 
              sparse. This, in principle, may cause problems
              when comparing across different proteins. However, in practice, we can proceed
              in caution given the fact that protein motions are not random or completely
              decorrelated in time.")
        )
      ),

      bsModal("gromosPseudo",
        title = "Gromos algorithm",
        trigger = "gromosLink",
        tags$ul(
          tags$li("Calculate all pairwise distances between objects"),
          tags$li("For a given cutoff distance find the configuration with the 
                    largest number of neighbors. Make it the centroid of the first 
                    cluster. Eliminate this and its neighbors from pool"),
          tags$li("Repeat for the configuration with the second largest number 
                    of neighbors"),
          tags$li("Repeat till pool is empty"),
          helpText(tags$div(
            HTML(paste("Time complexity O(N", tags$sup(2), ")", sep = ""))
          ))
        )
      ),
      bsModal("leaderPseudo",
        title = "Leader/follower algorithm",
        trigger = "leaderLink",
        tags$ul(
          tags$li("Assign first input as leader of first cluster"),
          tags$li("Read next input and calculate distance from existing leader(s)"),
          tags$li("Find the first leader with distance smaller than cutoff from 
                    input and assign input to that cluster. If none exists, make 
                    input a new leader"),
          tags$li("Repeat for all input"),
          helpText(tags$div(
            HTML(paste("Time complexity << O(N", tags$sup(2), "), depending 
                           on the choice of the cutoff.", sep = ""))
          ))
        )
      ),
      bsModal("kmeansPseudo",
        title = "K-means algorithm",
        trigger = "kmeansLink",
        tags$ul(
          tags$li("Choose number of clusters and initialize their components"),
          tags$li("Select input and alculate sum of squares between input and 
                    all cluster centroids"),
          tags$li("Assign input to the cluster with the 'nearest mean'"),
          tags$li("Update cluster centroids "),
          tags$li("Repeat for all input"),
          helpText(tags$div(
            HTML(paste("Time complexity O(N", tags$sup(2), ")", sep = ""))
          ))
        )
      )
    )
  ) # tabset
))
