###### PlaqView Master Code ######
##### Author: Wei Feng Ma, UVA.
##### wm5wt@virginia.edu

#### DATABASE NAMES AND COLOR SCHEMES ####
# below line is commented for shinyapp.io deployment temp
### set this once in terminal before deploying to shinyapps.io ###
# options(repos = BiocManager::repositories())

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

# color definitions
original_color_list <-
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

metcolors <- MetBrewer::met.brewer("Juarez", n = 6) 
# color_function <- colorRampPalette(original_color_list)
color_function <- colorRampPalette(metcolors)

manual_color_list <- color_function(40) # change this if clusters >40

#### LIBRARIES #### 
# library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(shinybusy) #install.packages("shinybusy")
library(enrichR) # install.packages("enrichR")
library(waiter)
library(DT)
library(readxl)
library(shinyWidgets)
library(shinyjs)
# library(RColorBrewer)
library(rDGIdb) # BiocManager::install("rDGIdb")
library(tidyverse)
# library(rsconnect)
library(monocle3)

#### UI ####
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # AESTHETICS####
  theme = shinytheme("flatly"),
  add_busy_bar(color = "#ff9142", height = "100px"), # THIS IS THE BUSY BAR
  use_waiter(), 
  waiter_show_on_load(html = spin_rotate()),
  useShinyjs(),
  
  # defining each 'tab' here
  navbarPage("PlaqView", id = "inTabset",
             
             #### FRONT PANEL Panel ####
             tabPanel("Select Dataset", 
                      mainPanel(width = 12,
                                fluidRow(
                                  column(width = 12,
                                         wellPanel(
                                           includeMarkdown("descriptionfiles/helptext_welcome.Rmd"),
                                           img(src = "abstract.png", width = '100%'),
                                         )),
                                  
                                )),
                      mainPanel(width = 12,
                                DT::dataTableOutput('availabledatasettable'),
                                br(),
                                fluidRow(
                                  column(width = 6,
                                         actionBttn(
                                           inputId = "loaddatabutton",
                                           label = "Step 1: Click Here to Load Dataset",
                                           style = "unite",
                                           color = "primary",
                                           block = T)
                                  ),
                                  column(width = 6,
                                         hidden(
                                           actionBttn(
                                             inputId = "jumpto1",
                                             label = "Step 2: Start PlaqView",
                                             color = "success",
                                             block = T)
                                         )
                                  )
                                ),
                                br(),
                                br(),
                      )
             ),
             
             #### PANEL 1: QUERY GENE   ----
             tabPanel(title = "Gene Lookup", value = "panel1",
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
                                        value = "TREM2, CYBB",
                                        placeholder = "try: TREM2, CYBB"
                                      ),
                                      
                                      # choose the type of output graph 
                                      pickerInput("selectaplot",
                                                  label = "Select Plot Type", 
                                                  choices = list(
                                                    "Dot Plot (up to 9 genes)" = "Dot",
                                                    "Feature Plot (up to 4 genes)" = "Feature",
                                                    "Ridge Plot (single gene)" = "Ridge"),
                                                  width = '90%',
                                                  selected = "Dot Plot"),
                                      
                                      pickerInput(
                                        inputId = "selectlabelmethodforgenequery",
                                        label = "Select Labeling Method", 
                                        choices = list (
                                          "Seurat_Clusters",
                                          # "scCATCH_Blood",
                                          # "scCATCH_BV",
                                          # "scCATCH_Heart",
                                          "Author_Provided",
                                          "SingleR.calls",
                                          "Seurat_with_Tabula_Ref"  
                                        ), 
                                        selected = "Seurat_with_Tabula_Ref",
                                        width = '90%' #neeed to fit this
                                      ),
                                      
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
                                          textOutput("selecteddatasetID"),  
                                          fluidRow( # top split rows
                                            column(width = 6, align = "center", 
                                                   plotOutput("umaps", width = "auto", height = '500px'),
                                                   br(),
                                                   downloadButton("downloadumapplot", "Download UMAP", width = '100%')
                                            ),
                                            
                                            column(width = 6, align="center", 
                                                   conditionalPanel('input.selectaplot=="Ridge"', 
                                                                    plotOutput("Ridge", width = "auto", height = '500px'),
                                                                    br(),
                                                                    downloadButton("downloadridgeplot", "Download Ridge Plot", width = '100%')
                                                   ),
                                                   
                                                   conditionalPanel('input.selectaplot=="Dot"', 
                                                                    plotOutput("Dot", width = "auto", height = '500px'),
                                                                    br(),
                                                                    downloadButton("downloaddotplot", "Download Dot Plot", width = '100%')
                                                   ), # conditional panels renders only if conditions are met
                                                   
                                                   conditionalPanel('input.selectaplot=="Feature"', 
                                                                    plotOutput("Feature", width = "auto", height = '500px'),
                                                                    br(),
                                                                    downloadButton("downloadfeatureplot", "Download Feature Plot", width = '100%')
                                                                    
                                                   )
                                                   
                                            ) # column 
                                            
                                          ), # fluidrow
                                          br(),
                                          fluidRow( # bottom whole for GO output
                                            column(width = 12, 
                                                   selectInput("selectedenrichRdb", label = h5("Top Significantly Enriched Pathways"), 
                                                               choices = enrichRdb, 
                                                               selected = "GO_Biological_Process_2018")
                                            ),
                                            column(width = 12, 
                                                   DT::dataTableOutput("enrichtable"),
                                                   helpText("You must restart query if you change database"),
                                                   downloadButton("downloadenrichRdata", "Download Complete Pathway Enrichment Data")
                                                   
                                            ),
                                          )# another fluidrow 
                                          
                                )# wellpanel
                                
                                
                      )# MAIN PANEL CLOSURE
             ), # TAB PANEL CLOSURE
             
             
             
             #### PANEL 2: COMPARE LABELING METHODS  ----  
             tabPanel("Cell Labeling",
                      mainPanel(width = 12, # 12/12 is full panel,
                                wellPanel(includeMarkdown("descriptionfiles/helptext_comparelabels.Rmd")),
                                wellPanel(
                                  fluidRow(
                                    column(width = 6, 
                                           
                                           pickerInput("leftlabeltype", 
                                                       label = "Select Labeling Method #1",
                                                       choices = list (
                                                         "Seurat_Clusters",
                                                         # "scCATCH_Blood",
                                                         # "scCATCH_BV",
                                                         # "scCATCH_Heart",
                                                         "Author_Provided",
                                                         "SingleR.calls",
                                                         "Seurat_with_Tabula_Ref"  
                                                       ), 
                                                       selected = "Seurat_with_Tabula_Ref"),
                                                       
                                           plotOutput("leftlabelplot",
                                                      height = '500px')),
                                    
                                    column(width = 6,
                                           pickerInput("rightlabeltype", 
                                                       label = "Select Labeling Method #2",
                                                       choices = list (
                                                         "Seurat_Clusters",
                                                         # "scCATCH_Blood",
                                                         # "scCATCH_BV",
                                                         # "scCATCH_Heart",
                                                         "Author_Provided",
                                                         "SingleR.calls",
                                                         "Seurat_with_Tabula_Ref"  
                                                       ), 
                                                       selected = "SingleR.calls"),
                                                       
                                           plotOutput("rightlabelplot",
                                                      height = '500px')
                                    ),# column
                                    
                                    br(), 
                                    
                                    column(width = 6, h4("Differential Expression by Cluster"),
                                           downloadButton("diffbyseurat", "Numbered Only (Unlabeled)"),
                                           downloadButton("diffbyauthor", "Author Supplied (Manual)"),
                                           
                                           
                                           helpText("This will download a .csv of differentially expressed genes by cluster")
                                    ), # column
                                    
                                    column(width = 6, h4("Differential Expression by Cell Type"),
                                           downloadButton("diffbysingleR", "SingleR"),
                                           downloadButton("diffbyTS", "Seurat + Tabula sapiens"),
                                           
                                           helpText("This will download a .csv of differentially expressed genes as identified by individual cells")
                                    ) # column
                                  ) # fluidrow
                                ) # close wellpanel
                      ), # mainPanel
                      
                      
                      
             ), # tabPanel
             #### PANEL 3: MONOCLE3 ----
             tabPanel("Trajectory",
                      mainPanel(width = 12, # 12/12 is full panel,
                                ### top panel
                                wellPanel(
                                  fluidRow(
                                  column(width = 6,
                                         includeMarkdown("descriptionfiles/helptext_comparetrajectories.Rmd"),
                                         # choose button
                                         actionButton("choose_toggle", "Choose/Unchoose"),
                                         # clear button
                                         actionButton("reset", "Clear Selection"),
                                         # recal button
                                         actionButton("redomonocle3", "Recalculate"),
                                         h4("Instructions:"),
                                         tags$ol(
                                           tags$li("Highlight points by clicking and dragging."),
                                           tags$li("Click the 'Choose/unchoose' button."),
                                           tags$li("Repeat until all of the desired cells are black."),
                                           tags$li("Click 'Done'.")
                                                ),
                                         h4("Details:"),
                                         tags$ul(
                                           tags$li("To start over, click 'Clear'"),
                                           tags$li(paste("You can also choose/unchoose specific cells",
                                                         "by clicking on them directly.")),
                                           ),
                                         ), # column
                                       column(width = 6, 
                                             
                                              plotOutput("plot1", click = "plot1_click", 
                                                         brush = brushOpts(id = "plot1_brush"), height = "600px"),
                                              actionButton("downloadoriginaltrajectory", label = "Download Original Trajectory"),
                                       ), # column
                                ), # fluidrow
                                ), # wellpanel
             
                                ### bottom panel
                                wellPanel(
                                  fluidRow(
                             
                                    plotOutput("originaltrajectory",
                                               width = '500px'),
                                    plotOutput("subsettrajectory",
                                               width = '500px'),

                                   
                                  )# fluidrow
                                ),
                      ), # mainPanel monocle3
                      
             ), # tabPanel
             
             #### PANEL 5: DRUGGABLE GENOME ----  
             tabPanel("Druggable Genome",
                      mainPanel(width = 12,
                                fluidRow(width = 12,
                                         column(width = 6,
                                                includeMarkdown("descriptionfiles/helptext_druggablegenome.Rmd"),
                                                wellPanel(
                                                  textInput(
                                                    inputId = "druggeneinput",
                                                    label = "Gene to Drug",
                                                    value = "EGFR"
                                                  ),
                                                  pickerInput("drugcelllabelmethod", 
                                                              label = "Select Labeling Method",
                                                              choices = list (
                                                                "Seurat_Clusters",
                                                                # "scCATCH_Blood",
                                                                # "scCATCH_BV",
                                                                # "scCATCH_Heart",
                                                                "Author_Provided",
                                                                "SingleR.calls",
                                                                "Seurat_with_Tabula_Ref"  
                                                              ), 
                                                              selected = "Seurat_with_Tabula_Ref")
                                                              
                                                ),
                                                
                                                checkboxGroupInput(
                                                  inputId = "dgidbdatabase",
                                                  label = "Choose a Database", 
                                                  inline = TRUE, 
                                                  selected = sourceDatabases(), # preselects all
                                                  choices = sourceDatabases()
                                                ),
                                                
                                                actionButton(inputId = "rundgidb",
                                                             label = "Start Query",
                                                             width = '100%')
                                         ), 
                                         br(),
                                         column(width = 6, 
                                                wellPanel(plotOutput("featurefordrugs",
                                                                     height = '500px')  )),
                                         br(),
                                         column(width = 12,
                                                DT::dataTableOutput("dgidboutput", width = "100%"),
                                                helpText("You must restart query if you change database. PubMed ID and citations of interactions are available in full download file."),
                                                br(),
                                                downloadButton("downloaddgidboutput", label = "Download Full Gene-Drug Interaction Table")
                                         )
                                )
                      )
                      
             ),
             
             
             ### JS to jump to top of page on click ###
             tags$script(" $(document).ready(function () {
                         $('#jumpto1').on('click', function (e) {
                         window.scrollTo(0, 0)                
                         });                 
                         });")
             
  )# close navbarpage
  
)# close fluidpage








#### SERVER ####
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  #### WELCOME PANEL ####
  df <- read_excel("Available_datasets.xlsx")
  df$DOI <- paste("<a href=",  df$DOI,">", "Link", "</a>") # this converts to clickable format
  # df <- column_to_rownames(df, var = "DataID")
  
  output$availabledatasettable <-
    DT::renderDataTable(df, server = F, # server is for speed/loading
                        selection = list(mode = 'single', selected = c(1)),
                        options=list(columnDefs = list(list(visible=FALSE, targets=c(11)))), # this hides the 8th column, which is datasetID
                        escape = FALSE) # this escapes rendering html (link) literally and makes link clickable
  
  observeEvent(input$loaddatabutton, {
    path <- file.path(paste("data/", df$DataID[input$availabledatasettable_rows_selected], "/",
                            df$DataID[input$availabledatasettable_rows_selected],
                            ".rds", sep=""))
    pathcds <- file.path(paste("data/", df$DataID[input$availabledatasettable_rows_selected], "/",
                            df$DataID[input$availabledatasettable_rows_selected],
                            "_cds.rds", sep=""))
    plaqviewobj <<- readRDS(file = path)
    plaqviewobj.cds <<- readRDS(file = pathcds)
    
    ## these are just for displaying current data name in other tabs##
    output$selecteddatasetID <- renderText({
      paste0("Current dataset: ", df$DataID[input$availabledatasettable_rows_selected])
    }) 
    output$selecteddatasetID2 <- renderText({
      paste0("Current dataset: ", df$DataID[input$availabledatasettable_rows_selected])
    }) 
    output$selecteddatasetID3 <- renderText({
      paste0("Current dataset: ", df$DataID[input$availabledatasettable_rows_selected])
    }) 
   
    show("jumpto1")
    
  })
  
  observeEvent(input$jumpto1, {
    updateTabsetPanel(session = getDefaultReactiveDomain(), "inTabset",
                      selected = "panel1") # this is to switch to tab1
    
  })

  
  #### PANEL #1 FUNCTIONS ####
  #### umap ###
  # UMAP plot#
  observeEvent(input$runcode,{ 
    output$umaps <- 
      renderPlot(
        DimPlot(
          plaqviewobj,
          reduction = "umap",
          label = TRUE,
          label.size = 5,
          repel = T,
          # repel labels
          pt.size = 1,
          cols = color_function(length(unique(plaqviewobj@meta.data[[input$selectlabelmethodforgenequery]]))),
                        # this enables dynamic # of colors based on # of labels given, sorry about the horrible nesting
          group.by = input$selectlabelmethodforgenequery) + # group.by is important, use this to call metadata separation
          theme(legend.position="bottom", 
                legend.box = "vertical") +
          ggtitle("UMAP by Cell Type") +
          theme(plot.title = element_text(hjust =  0.5)) +
          guides(color = guide_legend(nrow = 5)
          )
      ) # closes renderPlot
    
    #### NOMENCLATURE UPDATE ####
    if(df$Species[input$availabledatasettable_rows_selected] == "Human"){
      corrected <- str_to_upper(input$genes)
    } else{
      corrected <- str_to_title(input$genes)
    }
    
    updateTextInput(getDefaultReactiveDomain(),
                    "genes", value = corrected)
    
  })# closes observe event
  
  # this is for the download
  output$downloadumapplot<- downloadHandler(
    filename = function() {
      paste("UMAP.pdf", sep = "")
    },
    content = function(file) {
      pdf(file, paper = "default") # paper = defult is a4 size
      user_genes <- str_split(input$genes, ", ")[[1]]

      validate(need(input$selectaplot=="Dot", message=FALSE))
      temp <- DimPlot(
        plaqviewobj,
        reduction = "umap",
        label = TRUE,
        label.size = 5,
        repel = T,
        # repel labels
        pt.size = 1,
        cols = color_function(length(unique(plaqviewobj@meta.data[[input$selectlabelmethodforgenequery]]))),
        group.by = input$selectlabelmethodforgenequery) + # group.by is important, use this to call metadata separation
        theme(legend.position="bottom", 
              legend.box = "vertical") +
        ggtitle("UMAP by Cell Type") +
        theme(plot.title = element_text(hjust =  0.5)) +
        guides(color = guide_legend(nrow = 5))
      
      plot(temp) #this is all you need
      
      dev.off()
    }
    
  )# close downloadhandler
  
  #### dot plot ####
  observeEvent(input$runcode,{ # observe event puts a pause until pushed
    
    # this is for the display
    output$Dot <- renderPlot({
      # parse string input 
      user_genes <- str_split(input$genes, ", ")[[1]]
      validate(need(input$selectaplot=="Dot", message=FALSE))
      DotPlot(plaqviewobj, 
              group.by = input$selectlabelmethodforgenequery,
              features = user_genes) + # a trick to sep long string input
        ggtitle("Expression Dot Plot") +
        theme(plot.title = element_text(hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })
    
    # this is for the download
    output$downloaddotplot<- downloadHandler(
      filename = function() {
        paste("dotplot.pdf", sep = "")
      },
      content = function(file) {
        pdf(file, paper = "default") # paper = defult is a4 size
        user_genes <- str_split(input$genes, ", ")[[1]]
        validate(need(input$selectaplot=="Dot", message=FALSE))
        temp <- DotPlot(plaqviewobj, 
                        group.by = input$selectlabelmethodforgenequery,
                        features = user_genes) + # a trick to sep long string input
          ggtitle("Expression Dot Plot") +
          theme(plot.title = element_text(hjust = 1)) +
          theme(plot.title = element_text(hjust = 0.5)) 
        
        plot(temp) #this is all you need
        
        dev.off()
      }
      
    )# close downloadhandler
    
    #### feature plot ####
    parsed.genes <- str_split(input$genes, ", ")[[1]]
    output$Feature <- renderPlot({
      user_genes <- str_split(input$genes, ", ")[[1]]
      validate(need(input$selectaplot=="Feature", message=FALSE))
      FeaturePlot(plaqviewobj, 
                  features = user_genes[1:4]) + # a trick to sep long string input
        theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
        theme(plot.title = element_text(hjust = 1)) +
        theme(plot.title = element_text(hjust =  0.5)) 
    }) 
    
    # this is for the download
    output$downloadfeatureplot<- downloadHandler(
      filename = function() {
        paste("featureplot.pdf", sep = "")
      },
      content = function(file) {
        pdf(file, paper = "default") # paper = defult is a4 size
        user_genes <- str_split(input$genes, ", ")[[1]]
        validate(need(input$selectaplot=="Feature", message=FALSE))
        temp <- FeaturePlot(plaqviewobj, 
                            features = user_genes[1:4]) + # a trick to sep long string input
          theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
          theme(plot.title = element_text(hjust = 1)) +
          theme(plot.title = element_text(hjust =  0.5)) 
        
        plot(temp) #this is all you need
        
        dev.off()
      }
      
    )# close downloadhandler
    
    #### ridge plot ####
    output$Ridge <- renderPlot({
      user_genes <- str_split(input$genes, ", ")[[1]]
      validate(need(input$selectaplot=="Ridge", message=FALSE))
      RidgePlot(plaqviewobj,
                cols = color_function(length(unique(plaqviewobj@meta.data[[input$selectlabelmethodforgenequery]]))),
                group.by = input$selectlabelmethodforgenequery,
                features =  user_genes[1:1]
      ) + # a trick to sep long string input
        #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
        theme(plot.title = element_text(hjust =  0.5)) +
        guides(color = guide_legend(nrow = 3))
      
    })
    
    
    output$downloadridgeplot<- downloadHandler(
      filename = function() {
        paste("ridgeplot.pdf", sep = "")
      },
      content = function(file) {
        pdf(file, paper = "default") # paper = defult is a4 size
        user_genes <- str_split(input$genes, ", ")[[1]]
        validate(need(input$selectaplot=="Ridge", message=FALSE))
        temp <- RidgePlot(plaqviewobj,
                          cols = color_function(length(unique(plaqviewobj@meta.data[[input$selectlabelmethodforgenequery]]))),
                          group.by = input$selectlabelmethodforgenequery,
                          features =  user_genes[1:1]
        ) + # a trick to sep long string input
          #theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
          theme(plot.title = element_text(hjust =  0.5)) +
          guides(color = guide_legend(nrow = 3))
        
        plot(temp) #this is all you need
        
        dev.off()
      }
      
    )# close downloadhandler
    
    #### cell population ####
    output$cellpopulation <- renderPlot({
      # this chunk is to calculate population statistics 
      cellcounts <- plaqviewobj@meta.data[[input$selectlabelmethodforgenequery]] # extract cell names
      cellcounts <- as.factor(cellcounts) # convert to factors
      
      cellcounts.summary <- as.data.frame(summary(cellcounts)) # for levels laters
      
      cellcounts <- as.data.frame((cellcounts)) # summarize factors into table
      cellcounts <- dplyr::rename(cellcounts, Celltype = '(cellcounts)' ) # just to rename the columns
      
      # this is to create levels so i can flip the graph to the way i want
      cellcounts.summary <- rownames_to_column(cellcounts.summary)
      cellcounts.summary <- reorder(cellcounts.summary$rowname, cellcounts.summary$`summary(cellcounts)`)
      sortedlevel <- levels(cellcounts.summary)
      
      cellpop <- ggplot(data = cellcounts, aes(y = Celltype)) +
        geom_bar(fill = manual_color_list) +
        xlab("Counts") +
        ylab("Cell Types") +
        xlim(c(-150,1800)) +
        theme_light() +
        scale_y_discrete(limits=sortedlevel) +
        stat_count(geom = "text", # this stat_count function gives percentages
                   aes(label = paste(round((..count..)/sum(..count..)*100,2),"%; n =", round((..count..)))))
      cellpop
    })
    
    # Gene ontology table #
    {
      parsed.genes <- str_split(input$genes, ", ")[[1]]
      enriched <- enrichr(genes = parsed.genes, 
                          database = enrichRdb) # this queries all of them
      cleanedenrichedtable <- select(enriched[[input$selectedenrichRdb]], -Old.Adjusted.P.value, -Old.P.value,)
      cleanedenrichedtable <- top_n(cleanedenrichedtable, 100) # top 100 will be rendered
      output$enrichtable <- DT::renderDataTable(cleanedenrichedtable)
      
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
  
  
  
  #### PANEL #2 LABELS FUNCTIONS ####
  output$leftlabelplot <-
    renderPlot(
      DimPlot(
        plaqviewobj,
        reduction = "umap",
        label = TRUE,
        label.size = 5,
        repel = T,
        # repel labels
        pt.size = 1,
        cols = color_function(length(unique(plaqviewobj@meta.data[[input$leftlabeltype]]))),
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
        plaqviewobj,
        reduction = "umap",
        label = TRUE,
        label.size = 5,
        repel = T,
        # repel labels
        pt.size = 1,
        cols = color_function(length(unique(plaqviewobj@meta.data[[input$rightlabeltype]]))),
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
      file.copy(paste("data/", df$DataID[input$availabledatasettable_rows_selected], "/",
                      "diff_by_seurat.csv", sep = ""), file)
      
    }  )# close downloadhandler
  
  output$diffbyauthor <- downloadHandler(
    filename = "differential_markergenes_by_author_annotated_clusters.csv",
    content = function(file) {
      
      file.copy(paste("data/", df$DataID[input$availabledatasettable_rows_selected], "/",
                      "diff_by_author.csv", sep = ""), file)
      
    }  )# close downloadhandler
  
  output$diffbysingleR <- downloadHandler(
    filename = "differential_markergenes_by_singleR_labels.csv",
    content = function(file) {
      file.copy(paste("data/", df$DataID[input$availabledatasettable_rows_selected], "/",
                      "diff_by_singleR.csv", sep = ""), file)      
    }  )# close downloadhandler
  
  output$diffbyTS <- downloadHandler(
    filename = "differential_markergenes_by_tabulus_sapien_reference.csv",
    content = function(file) {
      file.copy(paste("data/", df$DataID[input$availabledatasettable_rows_selected], "/",
                      "diff_by_predicted.id_tabulus.sapien.csv", sep = ""), file)      
    }  )# close downloadhandler
  
  #### PANEL #3 TRAJ FUNCTIONS ####
  output$originaltrajectory <-
    renderPlot(
      plot_cells(plaqviewobj.cds,
                 color_cells_by = "assigned_cell_type",
                 label_groups_by_cluster=F,
                 show_trajectory_graph = T,
                 trajectory_graph_segment_size = 1,
                 graph_label_size = 1, # size of # in circle
                 group_label_size = 4)
      ) # renderplot
  
  #### selectable plot
  reduction_method <- "UMAP"

  observeEvent(input$loaddatabutton, {
    vals <<- reactiveValues(
      keeprows = rep(FALSE, nrow(colData(plaqviewobj.cds)))
    )
    reduced_dims <- as.data.frame(reducedDims(plaqviewobj.cds)[[reduction_method]])
    names(reduced_dims)[1:2] <- c("V1", "V2")
    
    
    vals <<- shiny::reactiveValues(
      keeprows = rep(FALSE, nrow(colData(plaqviewobj.cds)))
    )
    

    output$plot1 <<- renderPlot({
      # Plot the kept and excluded points as two separate data sets
      colData(plaqviewobj.cds)$keep <- vals$keeprows
      
      suppressMessages(plot_cells(plaqviewobj.cds, reduction_method = reduction_method,
                                  color_cells_by = "assigned_cell_type",
                                  label_groups_by_cluster=F,
                                  show_trajectory_graph = T,
                                  trajectory_graph_segment_size = 1,
                                  graph_label_size = 1, # size of # in circle
                                  group_label_size = 4,
                                  rasterize=FALSE) +
                         geom_point(alpha = colData(plaqviewobj.cds)$keep)) +
        theme(legend.position = "none")
    }, height = function() {
      session$clientData$output_plot1_width
    })
    
    # Toggle points that are clicked
    observeEvent(input$plot1_click, {
      res <- nearPoints(reduced_dims,
                        xvar = "V1", yvar = "V2", input$plot1_click,
                        allRows = TRUE)
      vals$keeprows <<- vals$keeprows | res$selected_
    })
    
    # Toggle points that are brushed, when button is clicked
    observeEvent(input$choose_toggle, {
      res <- brushedPoints(reduced_dims,
                           xvar = "V1", yvar = "V2", input$plot1_brush,
                           allRows = TRUE)
      vals$keeprows <- vals$keeprows | res$selected_
    })
    
    # Reset all points
    observeEvent(input$reset, {
      vals$keeprows <- rep(FALSE, nrow(colData(plaqviewobj.cds)))
    })
    
  })
  
  ### selected/recalulated trajectory
  observeEvent(input$redomonocle3, {
    
    # this is the selected cells
    selectedcells <- vals$keeprows
    
    # get selected cells id
    cds_subsetIDs <- row.names(colData(plaqviewobj.cds)[selectedcells,])
    
    # recreate a cds obj that doesnt contain any prior UMAPS
    expressiondata <- plaqviewobj@assays[["RNA"]]@data
    
    cellmd <- plaqviewobj@meta.data
    
    genemd <- data.frame(gene_short_name = row.names(expressiondata), 
                         row.names = row.names(expressiondata))

    plaqviewobj.cds_NEW <- new_cell_data_set(expression_data = expressiondata,
                                         cell_metadata = cellmd,
                                         gene_metadata = genemd)
    
    subsetted <- plaqviewobj.cds[,cds_subsetIDs]
    
    subsetted <- preprocess_cds(subsetted, num_dim = 30) # we used 30 in earlier seurat scripts
    
    
    # reproject cells now
    subsetted <- reduce_dimension(subsetted, reduction_method = "UMAP")
    subsetted <- cluster_cells(subsetted, reduction_method = "UMAP")
    
    subsetted <- learn_graph(subsetted)
    
    # a helper function to identify the root principal points:
    get_earliest_principal_node <- function(cds, assigned_cell_type= "SMCs"){ # change celltype if desired
      cell_ids <- which(colData(cds)[, "assigned_cell_type"] == assigned_cell_type)
      
      closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
      
      root_pr_nodes
    }
    
    subsetted <<- order_cells(subsetted, root_pr_nodes=get_earliest_principal_node(subsetted),
                                   reduction_method = "UMAP")
    # plot subset plots 
    output$subsettrajectory <-
      renderPlot(
        plot_cells(subsetted,
                   color_cells_by = "assigned_cell_type",
                   label_groups_by_cluster=F,
                   show_trajectory_graph = T,
                   trajectory_graph_segment_size = 1,
                   graph_label_size = 1, # size of # in circle
                   group_label_size = 4,
                   cell_size = 1,
                   alpha = 0.7,
                   scale_to_range = T)       ) # renderplot

  })

  #### PANEL #5 DRUG FUNCTIONS ####
  observeEvent(input$rundgidb, {
    
    #### NOMENCLATURE UPDATE ####
    if(df$Species[input$availabledatasettable_rows_selected] == "Human"){
      corrected <- str_to_upper(input$druggeneinput)
    } else{
      corrected <- str_to_title(input$druggeneinput)
    }
    
    updateTextInput(getDefaultReactiveDomain(),
                    "druggeneinput", # input ID of the textinput
                    value = corrected) # correct to
    
    # updated gene name is now here
    updated_druggeneinput <- input$druggeneinput
    
    druggenes <- str_split(updated_druggeneinput, ", ")[[1]]
    
    # set active. identity
    Idents(plaqviewobj) <- input$drugcelllabelmethod
  
    
    result <- queryDGIdb(updated_druggeneinput,
                         sourceDatabases = input$dgidbdatabase)
    fulltable <- result@data[["interactions"]][[1]]
    
    
    # so if table becomes a list (empty), run the following
    # this is a table to show no drugs available
    if (class(fulltable) == "list" ) {
      nodrug <- matrix(c("No Drugs Found"),ncol=4,byrow=TRUE)
      colnames(nodrug) <- c("drugName","interactionTypes","score","drugConceptId")
      nodrug <- as_data_frame(nodrug)
      
      isolatedtable <- nodrug
      fulltable <- nodrug
      
    }
    
    if (class(fulltable) == "data.frame" ) {
      isolatedtable <-  fulltable %>% # reorder columns
        select("Drug_Name" = drugName, 
               "Interaction_Types" = interactionTypes, # rename columns might as well :)
               "Int_Score" = score , 
               #pmids, 
               #sources, 
               "Drug_ConceptID" = drugConceptId
        )
      
      # these need to be run to 'flatten' the list-type columns
      fulltable <- fulltable %>% mutate(interactionTypes = map_chr(interactionTypes, toString))
      fulltable <- fulltable %>% mutate(sources = map_chr(sources, toString))
      fulltable <- fulltable %>% mutate(pmids = map_chr(pmids, toString))
    }
    
    # plots
    output$featurefordrugs <- renderPlot({
      user_genes <- str_split(corrected, ", ")[[1]]
      FeaturePlot(plaqviewobj, 
                  features = user_genes, label = T, repel = T
      ) + # a trick to sep long string input
        theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
        theme(plot.title = element_text(hjust = 1)) +
        theme(plot.title = element_text(hjust =  0.5)) 
    }) # render plot
    
    
    output$dgidboutput <- DT::renderDataTable(isolatedtable, )
    output$downloaddgidboutput <- downloadHandler(
      filename = function() {
        paste(updated_druggeneinput, "_complete_drug-gene_int.csv", sep = "")
      },
      content = function(file) {
        write.csv(fulltable, file, row.names = FALSE, col.names = T)
      } 
    )# close downloadhandler
    
    # complete druggable genome
    # finan.genome <- read_delim("/data/Druggable_genome_finan.txt", 
    #                            "\t", escape_double = FALSE, trim_ws = TRUE)
  })# observer event
  
  #### PANEL #6 ABOUT FUNCTIONS #### 
  
  output$downloadsessioninfo <- downloadHandler(
    filename = paste(date(), "sesssioninfo.txt"),
    content = function(file) {
      write_lines(sessionInfo(), file)
    }  )# close downloadhandler
  
  #### Waiter ####
  waiter_hide()
  
  

  
} # ends server function

# Run the application 
shinyApp(ui = ui, server = server)


#### NOTES ####
