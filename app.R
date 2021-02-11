#### LIBRARIES #### 
library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(shinybusy) #install.packages("shinybusy")
library(tidyverse)
library(enrichR) # install.packages("enrichR")

#### extrasensory codes ####


#### LOADING DATA ####
# below line is commented for shinyapp.io deployment temp
stanford <- readRDS(file = "data/final_stanford_extendedlabels_02082021.rds")
# stanford <- readRDS(file = url("https://virginia.box.com/shared/static/oyo1bicpvlxen940zmciqapvg0y3n6gb.rds"))

# enrichR functions
# handcurate db names 
dbs <- c("KEGG_2019_Human",
         "WikiPathways_2019_Human",
         "GO_Biological_Process_2018",
         "ChEA_2016",
         "GWAS_Catalog_2019",
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
         "Gene_Perturbations_from_GEO_down",
         "Gene_Perturbations_from_GEO_up")
enrichRdb <- sort(dbs)

#### SHINY OPTIONS, COLORs ####
# make the graphs match color of the UI
shinyOptions(plot.autocolors = TRUE)


# color definitions
manual_color_list <-
  {c("rosybrown2",
     "cadetblue1",
     "lemonchiffon3",
     "darkseagreen",
     "skyblue3",
     "thistle3",
     "cadetblue3",
     "darkseagreen1",
     "palevioletred3",
     "palevioletred1",
     "darkseagreen2",
     "rosybrown3",
     "thistle2",
     "lightsteelblue3",
     "salmon1",
     "palevioletred4",
     "lemonchiffon4",
     "cadetblue2"
  )}

#### UI ####
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # set theme
  theme = shinytheme("flatly"),
  add_busy_bar(color = "#ff9142"), # THIS IS THE BUSY BAR
  
  
  
  # defining each 'tab' here
  navbarPage("PlaqView", # title 
             
             # PANEL 1: QUERY GENE   ----
             tabPanel("Quick Gene Lookup", 
                      mainPanel(width = 12, # 12/12 is full panel
                                fluidRow(## panel for gene input
                                  column(
                                    width = 5,
                                    wellPanel(                                      
                                      # must add up to 12 all columns
                                      textInput(
                                        "genes",
                                        width = '100%',
                                        h3("Query Gene Expression", h5("please follow HUGO conventions")),
                                        value = "NOX4, CYBB",
                                        placeholder = "try: NOX4, CYBB"
                                      ),
                                      
                                      # choose the type of output graph 
                                      helpText("Preferred Plot Type"),
                                      
                                      selectInput("selectaplot", label = NULL, 
                                                  choices = list(
                                                    "Dot Plot (up to 9 genes)" = "Dot",
                                                    "Feature Plot (up to 4 genes)" = "Feature",
                                                    "Ridge Plot (single gene)" = "Ridge"),
                                                  width = '75%',
                                                  selected = "Dot Plot"),
                                      
                                      helpText("Preferred Cell Labeling Method"),
                                      selectInput("selectlabelmethodforgenequery", label = NULL, 
                                                  choices = list(
                                                    "Wirka et al. (Nature Med. 2019)" = "manually_annotated_labels",
                                                    "SingleR (Individual Cell ID)" = "SingleR.calls",
                                                    "Seurat Clusters (Numbered)" = "seurat_clusters",
                                                    "scCATCH (Heart)" = "scCATCH_Heart",
                                                    "scCATCH (Blood Vessels)" = "scCATCH_BV"),
                                                  width = '75%',
                                                  selected = "Wirka et al. (Nature Med. 2019)"),
                                      # 'go' button
                                      actionButton(
                                        inputId = "runcode",
                                        label = "Start Query",
                                        width = '100%'
                                      ))
                                    
                                  ),
                                  
                                  ## panel for description
                                  column(
                                    width = 7,
                                    wellPanel(includeMarkdown("descriptionfiles/helptext_singlegenepage.Rmd"))
                                  )
                                ),
                                
                                
                                #spacer
                                br(),
                                
                                
                                ## lower panel for graphic outputs
                                wellPanel(width = 12,
                                  fluidRow( # top splite rows
                                    column(width = 6, align="center", plotOutput("umaps", 
                                                                                 width = "auto",
                                                                                 height = '500px')),
                                    column(width = 6, align="center", 
                                           conditionalPanel('input.selectaplot=="Ridge"', plotOutput("Ridge",
                                                                                                     width = "auto",
                                                                                                     height = '500px')),
                                           conditionalPanel('input.selectaplot=="Dot"', plotOutput("Dot",
                                                                                                   width = "auto",
                                                                                                   height = '500px')), # conditonal panels renders only if conditions are met
                                           conditionalPanel('input.selectaplot=="Feature"', plotOutput("Feature",
                                                                                                       width = "auto",
                                                                                                       height = '500px'))
                                           
                                    ) # column 
                                    
                                  ), # fluidrow
                                  br(),
                                  fluidRow( # bottom whole for GO output
                                    column(width = 12, 
                                           selectInput("selectedenrichRdb", label = h5("Top 25 Pathway Enrichment"), 
                                                       choices = enrichRdb, 
                                                       selected = "GO_Biological_Process_2018")
                                    ),
                                    column(width = 12, 
                                           tableOutput("enrichtable"),
                                           helpText("You must restart query if you change database"),
                                           downloadButton("downloadenrichRdata", "Download Complete Pathway Enrichment Data")
                                           
                                    ),
                                  )# another fluidrow 
                                  
                                )# wellpanel
                                
                                
                      )# MAIN PANEL CLOSURE
             ), # TAB PANEL CLOSURE
             
             
             
             # PANEL 2: COMPARE LABELING METHODS  ----  
             tabPanel("Compare Labeling Methods",
                      mainPanel(width = 12, # 12/12 is full panel,
                                wellPanel(includeMarkdown("descriptionfiles/helptext_comparelabels.Rmd")),
                                wellPanel(
                                  fluidRow(
                                    column(width = 6, 
                                           selectInput("leftlabeltype", 
                                                       label = NULL,
                                                       choices = list(
                                                         "Wirka et al. (Nature Med. 2019)" = "manually_annotated_labels",
                                                         "SingleR (Individual Cell ID)" = "SingleR.calls",
                                                         "Seurat Clusters (Numbered)" = "seurat_clusters",
                                                         "scCATCH (Heart)" = "scCATCH_Heart",
                                                         "scCATCH (Blood Vessels)" = "scCATCH_BV"),
                                                       selected = "Wirka et al. (Nature Med. 2019)"),
                                           plotOutput("leftlabelplot",
                                                      height = '500px')),
                                    
                                    column(width = 6,
                                           selectInput("rightlabeltype", 
                                                       label = NULL,
                                                       choices = list(
                                                         "SingleR (Individual Cell ID)" = "SingleR.calls",
                                                         "Seurat Clusters (Numbered)" = "seurat_clusters",
                                                         "scCATCH (Heart)" = "scCATCH_Heart",
                                                         "scCATCH (Blood Vessels)" = "scCATCH_BV"),
                                                       selected = "SingleR (Individual Cell ID)"),
                                           plotOutput("rightlabelplot",
                                                      height = '500px')
                                    ),# column
                                    
                                    br(), 
                                    
                                    column(width = 6, h4("Differential Expression by Cluster"),
                                           downloadButton("diffbyseurat", "Numbered Only (Unlabeled)"),
                                           downloadButton("diffbywirka", "Wirka et al."),
                                           helpText("This will download a .csv of every cluster identified in singleR")
                                    ), # column
                                    
                                    column(width = 6, h4("Differential Expression by Cell Type (SingleR)"),
                                           downloadButton("diffbysingleR", "Download"),
                                           helpText("This will download a .csv of every cluster by cluster number only, intended for manually identifying cell type")
                                    ) # column
                                  ) # fluidrow
                                ) # close wellpanel
                      ), # mainPanel
                      
                      
                      
             ), # tabPanel
             # PANEL 3: COMPARE TRAJECTORY METHODS ----
             tabPanel("Compare Trajectory Methods",
                      mainPanel(width = 12, # 12/12 is full panel,
                                wellPanel(includeMarkdown("descriptionfiles/helptext_comparetrajectories.Rmd")),
                                wellPanel(
                                  fluidRow(
                                    column(width = 6,
                                           selectInput("lefttrajectory",
                                                       label = NULL,
                                                       choices = list(
                                                         "Monocle3  (Trapnell Lab)" = "monocle3",
                                                         "PAGA (Theis Lab)" = "paga",
                                                         "Slingshot (Dudoit Lab)" = "slingshot",
                                                         "SCORPIUS (Ginhoux Lab)" = "scorpius"
                                                       ),
                                                       selected = "Monocle3  (Trapnell Lab)"),
                                           plotOutput("lefttrajectory",
                                                      height = '500px')),

                                    
                                           column(width = 6,
                                                  selectInput("righttrajectory",
                                                              label = NULL,
                                                              choices = list(
                                                                "PAGA (Theis Lab)" = "paga",
                                                                "Slingshot (Dudoit Lab)" = "slingshot",
                                                                "SCORPIUS (Ginhoux Lab)" = "scorpius"
                                                              ),
                                                              selected ="PAGA (Theis Lab)"),
                                                  plotOutput("righttrajectory",
                                                             height = '500px')        
                                           ) # column
                                   

          

                        
                                  ) # fluidrow
                                ) # close wellpanel
                      ), # mainPanel



             ), # tabPanel


             # PANEL 4: EXPLORE YOUR OWN DATA ----  
             tabPanel("Explore Your Own Dataset",
                      mainPanel(width = 12, # 12/12 is full panel,
                                fileInput(inputId = "upload",
                                          label = "Upload .rds or count matrix",
                                          width = "100%",
                                          accept = c(".txt", ".rds")),
                                helpText("This feature is under development! Check back soon for updates.")
                      ), # close mainpanel
                      
                      
             ), # close tabPanel
             
             # PANEL 5: ABOUT PANEL ----
             tabPanel("About & Help",
                      mainPanel(
                        # descriptions
                        includeMarkdown("descriptionfiles/aboutusdescription.Rmd"),
                        br(),
                        img(src = "MSTPlogo.png", width = 233, height = 83),
                        br(),
                        br(),
                        downloadButton("downloadsessioninfo", "Download Session and Package Information")
                        
                        
                      )) # close tab panel
             
  )# close navbarpage
  
)# close fluidpage





#### SERVER ####
# Define server logic required to draw a histogram
server <- function(input, output) {
  #### PANEL #1 FUNCTIONS ####
  # UMAP plot, interactive #
  observeEvent(input$runcode,{ 
    output$umaps <- 
      renderPlot(
        DimPlot(
          stanford,
          reduction = "umap",
          label = TRUE,
          label.size = 4,
          repel = T,
          # repel labels
          pt.size = 1,
          cols = manual_color_list,
          group.by = input$selectlabelmethodforgenequery) + # group.by is important, use this to call metadata separation
          theme(legend.position="bottom", 
                legend.box = "vertical") +
          ggtitle("UMAP by Cell Type") +
          theme(plot.title = element_text(hjust =  0.5)) +
          guides(color = guide_legend(nrow = 5))
      ) # closes renderPlot
  })# closes observe event
  
  # Gene expression plots #
  # dot
  observeEvent(input$runcode,{ # observe event puts a pause until pushed
    output$Dot <- renderPlot({
      # parse string input 
      user_genes <- str_split(input$genes, ", ")[[1]]
      validate(need(input$selectaplot=="Dot", message=FALSE))
      DotPlot(stanford, 
              group.by = input$selectlabelmethodforgenequery,
              features = user_genes) + # a trick to sep long string input
        ggtitle("Expression Dot Plot") +
        theme(plot.title = element_text(hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5)) 
      
      
    })
    # feature
    parsed.genes <- str_split(input$genes, ", ")[[1]]
    output$Feature <- renderPlot({
      user_genes <- str_split(input$genes, ", ")[[1]]
      validate(need(input$selectaplot=="Feature", message=FALSE))
      FeaturePlot(stanford, 
                  features = user_genes[1:4]) + # a trick to sep long string input
        theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
        theme(plot.title = element_text(hjust = 1)) +
        theme(plot.title = element_text(hjust =  0.5)) 
      
    }) 
    #ridge
    output$Ridge <- renderPlot({
      user_genes <- str_split(input$genes, ", ")[[1]]
      validate(need(input$selectaplot=="Ridge", message=FALSE))
      RidgePlot(stanford,
                cols = manual_color_list,
                group.by = input$selectlabelmethodforgenequery,
                features =  user_genes[1:1]
                ) + # a trick to sep long string input
        #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
        theme(plot.title = element_text(hjust =  0.5)) +
        guides(color = guide_legend(nrow = 3))
      
    })
    
    
    # Gene ontology table #
    {
      parsed.genes <- str_split(input$genes, ", ")[[1]]
      enriched <- enrichr(genes = parsed.genes, 
                          database = enrichRdb) # this queries all of them
      cleanedenrichedtable <- select(enriched[[input$selectedenrichRdb]], -Old.Adjusted.P.value, -Old.P.value,)
      cleanedenrichedtable <- top_n(cleanedenrichedtable, 25) # top 25 will be rendered
      output$enrichtable <- renderTable(cleanedenrichedtable,
                                        striped = T,
                                        spacing = "xs",
                                        align = 'l',
                                        width = "98%",
                                        colnames = T,
                                        digits = 3)
      
      # Downloadable csv of selected dataset
      output$downloadenrichRdata <- downloadHandler(
        filename = function() {
          paste(input$genes, "_pathwayenrichment.csv", sep = "")
        },
        content = function(file) {
          write.csv(enriched[[input$selectedenrichRdb]], file, row.names = FALSE)
        } 
      )# close downloadhandler
    }
    
    
  }) # observe event closure
  
  
  
  #### PANEL #2 FUNCTIONS ####
  output$leftlabelplot <-
    renderPlot(
      DimPlot(
        stanford,
        reduction = "umap",
        label = TRUE,
        label.size = 3,
        repel = T,
        # repel labels
        pt.size = 1,
        cols = manual_color_list,
        group.by = input$leftlabeltype ) + # group.by is important, use this to call metadata separation
        theme(legend.position="bottom", 
              legend.box = "vertical") +
        ggtitle("UMAP by Cell Type") +
        theme(plot.title = element_text(hjust =  0.5)) +
        guides(color = guide_legend(nrow = 5))
    )# render plot
  
  output$rightlabelplot <- 
    renderPlot(
      DimPlot(
        stanford,
        reduction = "umap",
        label = TRUE,
        label.size = 3,
        repel = T,
        # repel labels
        pt.size = 1,
        cols = manual_color_list,
        group.by = input$rightlabeltype ) + # group.by is important, use this to call metadata separation
        theme(legend.position="bottom", 
              legend.box = "vertical") +
        ggtitle("UMAP by Cell Type") +
        theme(plot.title = element_text(hjust =  0.5)) +
        guides(color = guide_legend(nrow = 5))
    )# render plot
  
  ## download diffex
  output$diffbyseurat <- downloadHandler(
    filename = "differential_markergenes_by_seurat_clusters.csv",
    content = function(file) {
      
    }  )# close downloadhandler
  
  output$diffbywirka <- downloadHandler(
    filename = "Original_wirka_et_al_2019_NatureMed_supplemental.xlms",
    content = function(file) {
      
    }  )# close downloadhandler
  
  output$diffbysingleR <- downloadHandler(
    filename = "differential_markergenes_by_singleR_labels.csv",
    content = function(file) {
    }  )# close downloadhandler
  
  # #### PANEL #3 FUNCTIONS ####
  # output$lefttrajectory <-
  #   renderPlot(
  #    
  #   )# render plot
  # 
  # output$righttrajectory <- 
  #   renderPlot(
  #    
  #   )# render plot
  # 
  # 
  #### PANEL #4 FUNCTIONS #### 
  
  output$downloadsessioninfo <- downloadHandler(
    filename = paste(date(), "sesssioninfo.txt"),
    content = function(file) {
      write_lines(sessionInfo(), file)
    }  )# close downloadhandler
  
  
  
} # ends server function

# Run the application 
shinyApp(ui = ui, server = server)


#### NOTES ####
# PUT IN P VALUES
# GENE BASE RESULTS/TESTS TO HELP INTERPRET GWAS SIGNAL
# DOWNLOAD MORAN'S I AND OTHER SUMMARY TABLE
# DOWNLOAD DEG STATISTICS ON PG 1
# https://bioconductor.org/packages/release/bioc/html/Qtlizer.html includsion
# https://mrcieu.github.io/gwasglue/ # FUTURE ADDITIONAL
# PATHWAY ENRICHMENT RESULTS? -> CAUSAL GWASGLUE
# PUT IN 'COMING SOON' FEATURE MAP IN ABOUT SECTION TO PROTECT OURSELF IN NEXT VERSION
### E.G. GWAS LOOK UP ETC.