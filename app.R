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
stanford <- readRDS(file = "data/final_stanford_labeled.rds")
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
                                    width = 4,
                                    wellPanel(                                      
                                      # must add up to 12 all columns
                                      textInput(
                                        "genes",
                                        width = '100%',
                                        h3("Query Gene Expression", h5("please follow HUGO conventions")),
                                        placeholder = "try: NOX4, CYBB",
                                        value = "NOX4, CYBB"
                                      ),
                                      
                                      # choose the type of output graph 
                                      selectInput("selectaplot", label = h4("Plot Type"), 
                                                  choices = list(
                                                    "Dot Plot (one or more genes)" = "Dot",
                                                    "Feature Plot (single gene)" = "Feature",
                                                    "Ridge Plot (single gene)" = "Ridge"), 
                                                  selected = "Dot Plot"),                                        
                                      # 'go' button
                                      actionButton(
                                        inputId = "runcode",
                                        label = "Start Query",
                                        width = '100%'
                                      ))
                                    
                                  ),
                                  
                                  ## panel for description
                                  column(
                                    width = 8,
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
                                           selectInput("selectedenrichRdb", label = h5("Pathway Enrichment Database"), 
                                                       choices = enrichRdb, 
                                                       selected = "GO_Biological_Process_2018")
                                    ),
                                    column(width = 12, tableOutput("enrichtable"),
                                           downloadButton("downloadenrichRdata", "Download Pathway Enrichment Data")
                                           
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
                                                         "SingleR (Default)" = "SingleR.calls",
                                                         "Seurat Clusters" = "seurat_clusters",
                                                         "scCATCH (Heart)" = "scheart",
                                                         "scCATCH (Blood Vessels)" = "scvessels"
                                                       ),
                                                       selected = "SingleR (Default)"),
                                           plotOutput("leftlabelplot",
                                                      height = '500px')),
                                    
                                    column(width = 6,
                                           selectInput("rightlabeltype", 
                                                       label = NULL,
                                                       choices = list(
                                                         "Seurat Clusters" = "seurat_clusters",
                                                         "scCATCH (Heart)" = "scheart",
                                                         "scCATCH (Blood Vessels)" = "scvessels"
                                                       ),
                                                       selected = "Seurat Clusters"),
                                           plotOutput("rightlabelplot",
                                                      height = '500px')
                                    ),# column
                                    
                                    br(), 
                                    
                                    column(width = 6, h4("Differential Expression by Cell Type (SingleR)"),
                                           downloadButton("diffbysingleR", "Download"),
                                           helpText("This will download a .csv of every cluster identified in singleR")
                                    ), # column
                                    
                                    column(width = 6, h4("Differential Expression by Seurat (Unlabeled)"),
                                           downloadButton("diffbyseurat", "Download"),
                                           helpText("This will download a .csv of every cluster by cluster number only, intended for manually identifying cell type")
                                    ) # column
                                  ) # fluidrow
                                ) # close wellpanel
                      ), # mainPanel
                      
                      
                      
             ), # tabPanel
             
             # PANEL 3: EXPLORE YOUR OWN DATA ----  
             tabPanel("Explore Your Own Dataset",
                      mainPanel(width = 12, # 12/12 is full panel,
                                fileInput(inputId = "upload",
                                          label = "Upload .rds or count matrix",
                                          width = "100%",
                                          accept = c(".txt", ".rds")),
                                helpText("This feature is still under development.")
                      ), # close mainpanel
                      
                      
             ), # close tabPanel
             
             # PANEL 4: ABOUT PANEL ----
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
          label.size = 3,
          repel = T,
          # repel labels
          pt.size = 1,
          cols = manual_color_list,
          group.by = "SingleR.calls" ) + # group.by is important, use this to call metadata separation
          theme(legend.position="bottom", 
                legend.box = "vertical") +
          ggtitle("UMAP by Cell Type") +
          theme(plot.title = element_text(hjust =  0.5)) +
          guides(color = guide_legend(nrow = 5))
      ) # closes renderPlot
  })# closes observe event
  
  # Gene feature plot, interactive #
  observeEvent(input$runcode,{ # observe event puts a pause until pushed
    output$Dot <- renderPlot({
      validate(need(input$selectaplot=="Dot", message=FALSE))
      DotPlot(stanford, 
              features = (str_split(input$genes, ", "))[[1]])+ # a trick to sep long string input
        ggtitle("Expression Dot Plot") +
        theme(plot.title = element_text(hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5)) 
      
    })
    parsed.genes <- str_split(input$genes, ", ")[[1]]
    output$Feature <- renderPlot({
      validate(need(input$selectaplot=="Feature", message=FALSE))
      FeaturePlot(stanford, 
                  features = (str_split(input$genes, ", "))[[1]])+ # a trick to sep long string input
        theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
        theme(plot.title = element_text(hjust = 1)) +
        ggtitle("Expression Feature Plot") +
        theme(plot.title = element_text(hjust =  0.5)) 
      
    }) 
    
    output$Ridge <- renderPlot({
      validate(need(input$selectaplot=="Ridge", message=FALSE))
      RidgePlot(stanford,
                cols = manual_color_list,
                group.by = "SingleR.calls",
                features = (str_split(input$genes, ", "))[[1]],)+ # a trick to sep long string input
        theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
        ggtitle("Expression Ridge Plot") +
        theme(plot.title = element_text(hjust =  0.5)) +
        guides(color = guide_legend(nrow = 5))
      
    })
    
    
    # Gene ontology check #
    {
      parsed.genes <- str_split(input$genes, ", ")[[1]]
      enriched <- enrichr(genes = parsed.genes, 
                          database = enrichRdb) # this queries all of them
      cleanedenrichedtable <- select(enriched[[input$selectedenrichRdb]], -Old.Adjusted.P.value, -Old.P.value,)
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
  
  #### PANEL #3 FUNCTIONS ####
  
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