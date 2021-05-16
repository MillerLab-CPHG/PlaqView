#### LIBRARIES #### 
# library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(shinybusy) #install.packages("shinybusy")
library(tidyverse)
library(enrichR) # install.packages("enrichR")
library(imager)
library(waiter)
library(DT)
library(readxl)
library(shinyWidgets)
library(shinyjs)
library(RColorBrewer)
library(rDGIdb)

#### PreReq Codes ####
# below line is commented for shinyapp.io deployment temp

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

color_function <- colorRampPalette(original_color_list)
manual_color_list <- color_function(40) # change this if clusters >40

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
             
             #### WELCOME Panel ####
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
             
             # PANEL 1: QUERY GENE   ----
             tabPanel(title = "Quick Gene Lookup", value = "panel1",
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
                                                    "Author Supplied" = "manually_annotated_labels",
                                                    "SingleR (Individual Cell ID)" = "SingleR.calls",
                                                    "Seurat + Tabula sapiens Ref" = "predicted.id_tabulus.sapien",
                                                    "Seurat Clusters (Numbered)" = "seurat_clusters",
                                                    "scCATCH (Heart)" = "scCATCH_Heart",
                                                    "scCATCH (Blood Vessels)" = "scCATCH_BV"),
                                                  width = '75%',
                                                  selected = "Author Supplied"),
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
             
             
             
             # PANEL 2: COMPARE LABELING METHODS  ----  
             tabPanel("Labeling Methods",
                      mainPanel(width = 12, # 12/12 is full panel,
                                wellPanel(includeMarkdown("descriptionfiles/helptext_comparelabels.Rmd")),
                                wellPanel(
                                  fluidRow(
                                    column(width = 6, 
                                           selectInput("leftlabeltype", 
                                                       label = NULL,
                                                       choices = list(
                                                         "Author Supplied (Manual)" = "manually_annotated_labels",
                                                         "SingleR (Individual Cell ID)" = "SingleR.calls",
                                                         "Seurat + Tabula sapiens Ref" = "predicted.id_tabulus.sapien",
                                                         
                                                         "Seurat Clusters (Numbered)" = "seurat_clusters",
                                                         "scCATCH (Heart)" = "scCATCH_Heart",
                                                         "scCATCH (Blood Vessels)" = "scCATCH_BV"),
                                                       selected = "Author Supplied (Manual)"),
                                           plotOutput("leftlabelplot",
                                                      height = '500px')),
                                    
                                    column(width = 6,
                                           selectInput("rightlabeltype", 
                                                       label = NULL,
                                                       choices = list(
                                                         "SingleR (Individual Cell ID)" = "SingleR.calls",
                                                         "Seurat Clusters (Numbered)" = "seurat_clusters",
                                                         "Seurat + Tabula sapiens Ref" = "predicted.id_tabulus.sapien",
                                                         
                                                         "scCATCH (Heart)" = "scCATCH_Heart",
                                                          "scCATCH (Blood Vessels)" = "scCATCH_BV"),
                                                       selected = "SingleR (Individual Cell ID)"),
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
                                           helpText("This will download a .csv of differentially expressed genes as identified by singleR")
                                    ) # column
                                  ) # fluidrow
                                ) # close wellpanel
                      ), # mainPanel
                      
                      
                      
             ), # tabPanel
             # PANEL 3: COMPARE TRAJECTORY METHODS ----
             tabPanel("Trajectory Methods",
                      mainPanel(width = 12, # 12/12 is full panel,
                                wellPanel(includeMarkdown("descriptionfiles/helptext_comparetrajectories.Rmd")),
                                wellPanel(
                                  fluidRow(
                                    column(width = 6,
                                           selectInput("lefttrajectorymethod",
                                                       label = NULL,
                                                       choices = list(
                                                         "Monocle3  (Trapnell Lab)" = "monocle3",
                                                         "PAGA (Theis Lab)" = "paga",
                                                         "Slingshot (Dudoit Lab)" = "slingshot",
                                                         "SCORPIUS (Saeys Lab)" = "scorpius"
                                                       ),
                                                       selected = "Monocle3  (Trapnell Lab)"),
                                           plotOutput("lefttrajectorygraph",
                                                      height = '500px')
                                           
                                    ), # column
                                    
                                    column(width = 6,
                                           selectInput("righttrajectorymethod",
                                                       label = NULL,
                                                       choices = list(
                                                         "PAGA (Theis Lab)" = "paga",
                                                         "Slingshot (Dudoit Lab)" = "slingshot",
                                                         "SCORPIUS (Saeys Lab)" = "scorpius"
                                                       ),
                                                       selected ="PAGA (Theis Lab)"),
                                           plotOutput("righttrajectorygraph",
                                                      height = '500px')        
                                    ) # column
                                    
                                    
                                    
                                    
                                    
                                  ) # fluidrow
                                ) # close wellpanel
                      ), # mainPanel
                      
                      
                      
             ), # tabPanel
             
             # PANEL 4: DRUGGABLE GENOME ----  
             tabPanel("Drug-Gene Interactions",
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
                                 selectInput("drugcelllabelmethod", 
                                             label = NULL,
                                             choices = list(
                                               "SingleR (Individual Cell ID)" = "SingleR.calls",
                                               "Author Supplied (Manual)" = "manually_annotated_labels",
                                               "Seurat + Tabula sapiens Ref" = "predicted.id_tabulus.sapien",
                                               
                                               "Seurat Clusters (Numbered)" = "seurat_clusters",
                                               "scCATCH (Heart)" = "scCATCH_Heart",
                                               "scCATCH (Blood Vessels)" = "scCATCH_BV"),
                                             selected = "SingleR (Individual Cell ID)"),
                                 ),
                                
                                 checkboxGroupInput(
                                   inputId = "dgidbdatabase",
                                   label = "Choose a Database", 
                                   inline = TRUE, 
                                   selected = c("FDA", "DrugBank"),
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
             
             # PANEL 6: ABOUT PANEL ----
             tabPanel("About & Help",
                      mainPanel(
                        # descriptions
                        includeMarkdown("descriptionfiles/aboutusdescription.Rmd"),
                        br(),
                        img(src = "MSTPlogo.png", width = 233, height = 83),
                        img(src = "PlaqOmics.png", width = 233, height = 83),
                        img(src = "Leducq.png", width = 233, height = 83),
                        
                        br(),
                        br(),
                        downloadButton("downloadsessioninfo", "Download Session and Package Information")
                        
                        
                      )) # close tab panel

             
  )# close navbarpage
  
)# close fluidpage





            
#### SERVER ####
# Define server logic required to draw a histogram
server <- function(input, output) {
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
    plaqviewobj <<- readRDS(file = path)
    
    output$selecteddatasetID <- renderText({
      paste0("Sucessfully loaded the ", rownames(df)[input$availabledatasettable_rows_selected], 
             " dataset!.",
             collapse = ", ")
    }) 
    show("jumpto1")
    
  })
  
  observeEvent(input$jumpto1, {
    updateTabsetPanel(session = getDefaultReactiveDomain(), "inTabset",
                      selected = "panel1") # this is to switch to tab1
    
  })
  
  
  #### PANEL #1 FUNCTIONS ####
  #### umap ####
  # UMAP plot, interactive #
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
          cols = manual_color_list,
          group.by = input$selectlabelmethodforgenequery) + # group.by is important, use this to call metadata separation
          theme(legend.position="bottom", 
                legend.box = "vertical") +
          ggtitle("UMAP by Cell Type") +
          theme(plot.title = element_text(hjust =  0.5)) +
          guides(color = guide_legend(nrow = 5)
                 )
      ) # closes renderPlot
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
        cols = manual_color_list,
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
                cols = manual_color_list,
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
                          cols = manual_color_list,
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
  
  
  
  #### PANEL #2 FUNCTIONS ####
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
        plaqviewobj,
        reduction = "umap",
        label = TRUE,
        label.size = 5,
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
  
  #### PANEL #3 FUNCTIONS ####
  output$lefttrajectorygraph <- renderPlot(
    readRDS(file =  normalizePath(file.path('data/', # this will dynamically read the file containing the right name
                                            df$DataID[input$availabledatasettable_rows_selected], "/",
                                            paste(input$lefttrajectorymethod, '.rds', sep=''))))
    
  ) # renderplot 
  
  output$righttrajectorygraph <- renderPlot(
    readRDS(file =  normalizePath(file.path('data/', # this will dynamically read the file containing the right name
                                            df$DataID[input$availabledatasettable_rows_selected], "/",
                                            paste(input$righttrajectorymethod, '.rds', sep='')))) 
  ) # renderplot 
  
  
  #### PANEL #4 FUNCTIONS ####
  observeEvent(input$rundgidb, {
    druggenes <- str_split(input$druggeneinput, ", ")[[1]]
    Idents(plaqviewobj) <- input$drugcelllabelmethod
    output$featurefordrugs <- renderPlot({
      user_genes <- str_split(input$genes, ", ")[[1]]
      FeaturePlot(plaqviewobj, 
                  features = druggenes[1:1], label = T, repel = T
                  ) + # a trick to sep long string input
        theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
        theme(plot.title = element_text(hjust = 1)) +
        theme(plot.title = element_text(hjust =  0.5)) 
    }) # render plot
    
    result <- queryDGIdb(input$druggeneinput,
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
    
    
    output$dgidboutput <- DT::renderDataTable(isolatedtable, )
    output$downloaddgidboutput <- downloadHandler(
      filename = function() {
        paste(input$druggeneinput, "_complete_drug-gene_int.csv", sep = "")
      },
      content = function(file) {
        write.csv(fulltable, file, row.names = FALSE, col.names = T)
      } 
    )# close downloadhandler
    
    # complete druggable genome
    finan.genome <- read_delim("Druggable_genome_finan.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
  })# observer event
  
  #### PANEL #6 FUNCTIONS #### 
  
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
# PUT IN P VALUES
# GENE BASE RESULTS/TESTS TO HELP INTERPRET GWAS SIGNAL
# DOWNLOAD MORAN'S I AND OTHER SUMMARY TABLE
# DOWNLOAD DEG STATISTICS ON PG 1
# https://bioconductor.org/packages/release/bioc/html/Qtlizer.html includsion
# https://mrcieu.github.io/gwasglue/ # FUTURE ADDITIONAL
# PATHWAY ENRICHMENT RESULTS? -> CAUSAL GWASGLUE
# PUT IN 'COMING SOON' FEATURE MAP IN ABOUT SECTION TO PROTECT OURSELF IN NEXT VERSION
### E.G. GWAS LOOK UP ETC.