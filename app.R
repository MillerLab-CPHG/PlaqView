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

color_function <- colorRampPalette(original_color_list)
# color_function <- colorRampPalette(metcolors)

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
library(ggpubr)
library(gtools)
library(CIPR)
# library(reactlog)
# library(future)
# 
# # tell shiny to log all reactivity
# reactlog_enable()
# 
# # tell shiny to try to paralle compute
# future::plan("multisession")

#### READ GOOGLE SHEET ####
googlesheets4::gs4_deauth() # this tells google sheet to read-only
df <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1hLyjPFA2ZRpBLHnTgUnmDz7kimMZWFbz_ZGTml3-hRA/edit#gid=0")
# df <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vTF5Gw4Dbshlh3wVB8UAMswUEiOn4NEzXaEp8x73NtbWY3n4oIrWEVNMIwNYyInJM7k70G1lUcr7x9g/pub?output=csv")

df$DOI <- paste("<a href=",  df$DOI,">", "Link", "</a>") # this converts to clickable format

# subset data rows that are marked 'deployed = Yes"
df <- filter(df, `Deployed` == "Yes")
df <- df %>% 
  select('DataID', Year, Journal, DOI, Species, Tissue, Notes, Population, Cells = Cell.Number, `Article.Title` ) 
df$`Article.Title` <- str_to_title(df$`Article.Title`) # autocaps

#### UI ####
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Aesthetics ####
  theme = shinytheme("flatly"),
  add_busy_bar(color = "#ff9142", height = "100px"), # THIS IS THE BUSY BAR
  use_waiter(), 
  # waiter_show_on_load(html = spin_rotate()),
  useShinyjs(),
  # div(DT::dataTableOutput("table"), style = "font-size: 75%; width: 75%"), # DT font sinzes
  
  # Pages ####
  navbarPage("PlaqView", id = "inTabset",
             
             #### UI: Data ####
             tabPanel("Select Dataset", 
                      mainPanel(width = 12,
                                fluidRow(
                                  column(width = 5,
                                         wellPanel(
                                           img(src = "logo.png", width = '100%'),
                                           h3("Instructions:"),
                                           tags$ol(
                                             tags$li("Select a dataset from the BLUE drop-down Menu."),
                                             tags$li("Click the blue 'Load Dataset' button."), # change server code to toggle
                                             tags$li("The 'Start Exploring' button will appear when data is loaded."),
                                             tags$li("(Optional) come back to this page to load another dataset.")
                                             
                                           ),
                                           br(),
                                           fluidRow(
                                             column(width = 12,
                                                    # load data button
                                                    h4("Select a Dataset"),
                                                    pickerInput(
                                                      inputId = "dataselector",
                                                      #label = "Select a Dataset", 
                                                      choices = df$DataID,
                                                      choicesOpt = list(
                                                        subtext = paste(df$Species, ": ",
                                                                        df$Cells,
                                                                        " Cells",
                                                                        sep = "")),
                                                      options = list(
                                                        style = "btn-primary")
                                                    ),
                                                    
                                                    actionBttn(
                                                      inputId = "loaddatabutton",
                                                      label = "Load Dataset",
                                                      style = "unite",
                                                      color = "primary",
                                                      block = T),
                                                    
                                                  
                                                    br(),
                                                    # jump to page 1 button
                                                    hidden(
                                                      actionBttn(
                                                        inputId = "jumpto1",
                                                        label = "Start Exploring",
                                                        style = "unite",
                                                        color = "success",
                                                        block = T)
                                                    ),
                                                    br(),
                                                    helpText(textOutput("loadeddatasetID")),
                                                    
                                             ),
                                           
                                                    
                                             
                                           ),
                                           )
                                  ),
                                  column(width = 7,
                                         wellPanel(
                                           includeMarkdown("descriptionfiles/helptext_welcome.Rmd"),
                                           img(src = "abstract.png", width = '100%'),
                                         )
                                         ), # column
                                  
                                 
                                  
                                ) # fluid row 
                                ), # mainpanel
                      mainPanel(width = 12,
                                wellPanel(
                                  h4("Details of Single- Cell Dataset and IDs"),
                                  actionButton(inputId = "refreshtable", "Fetch Latest Dataset Details"),
                                  br(),
                                  br(),
                                  DT::dataTableOutput('availabledatasettable'),
                                  br(),
                                  inlineCSS(list("table" = "font-size: 12px")),
                                          ),
                              
                      )
             ),
             
             #### UI: Genes   ----
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
                                        value = "APOE, COL1A1, FBLN1, FBLN2",
                                        placeholder = "try: TREM2, CYBB"
                                      ),
                                      
                                      # choose the type of output graph 
                                      pickerInput("selectaplot",
                                                  label = "Select Plot Type", 
                                                  choices = list(
                                                    "Dot Plot (up to 9 genes)" = "Dot",
                                                    "Feature Plot (up to 4 genes)" = "Feature",
                                                    "Ridge Plot (single gene)" = "Ridge"),
                                                  width = '95%',
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
                                          "SingleR_calls" = "SingleR.calls",
                                          "Seurat_with_Tabula_Ref"  
                                        ), 
                                        selected = "Seurat_with_Tabula_Ref",
                                        width = '95%' #neeed to fit this
                                      ),
                                      br(),
                                      # 'go' button
                                      actionBttn(
                                        inputId = "runcode",
                                        label = "Start Query",
                                        style = "unite",
                                        color = "success",
                                        block = T)
                                    
                                      )
                                    
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
                                          # textOutput("selecteddatasetID"),  
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
             
             
             
             #### UI: Labels/CIPR  ----  
             tabPanel("Cell Labeling/CIPR",
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
                                                         "SingleR.calls" = "SingleR_calls",
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
                                                         "SingleR_calls" = "SingleR.calls",
                                                         "Seurat_with_Tabula_Ref"  
                                                       ), 
                                                       selected = "SingleR_calls"),
                                                       
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
                                ), # close wellpanel
                                
                                wellPanel(
                                  fluidRow(
                                    # CIPR reference selectors
                                    column(width = 5, 
                                           tags$h3(tags$b("CIPR: Cluster Identity Predictor")),
                                           h5("Use this module to compare our annotation against CIPR references."),
                                           
                                           h4("Instructions:"),
                                           tags$ol(
                                             tags$li("Choose the PlaqView Annotation Method."),
                                             tags$li("Choose the CIPR Reference Database."), 
                                             tags$li("Choose the CIPR Method (Calculates Cell Identity Score)."), 
                                             tags$li("Click Run CIPR."),
                                             tags$li("Optional: Drag a Box on CIPR Graph and Select Clusters for More Details.")
                                           ),
                                           
                                           pickerInput("CIPRoriginal", 
                                                       label = "Select PlaqView Labeling Method to Benchmark",
                                                       choices = list (
                                                         "Seurat_Clusters" = "diff_by_seurat.csv",
                                                         "Author_Provided" = "diff_by_author.csv",
                                                         "SingleR_calls" = "diff_by_singleR.csv",
                                                         "Seurat_with_Tabula_Ref" = "diff_by_Seurat_with_Tabula_Ref.csv"
                                                       ), 
                                                       selected = "SingleR_calls"),
                                           pickerInput("CIPRreference", 
                                                       label = "Select CIPR Reference", 
                                                       choices = list("ImmGen (mouse)" = "immgen",
                                                                   "Presorted RNAseq (mouse)" = "mmrnaseq",
                                                                   "Blueprint-Encode (human)" = "blueprint",
                                                                   "Primary Cell Atlas (human)" = "hpca",
                                                                   "DICE (human)" = "dice",
                                                                   "Hematopoietic diff (human)" = "hema",
                                                                   "Presorted RNAseq (human)" = "hsrnaseq"), 
                                                       selected = "Primary Cell Atlas (human)"), 
                                           
                                           # pickerInput("CIPRmethod", 
                                           #              label = "Select CIPR Method", 
                                           #              choices = list("logFC Dot Product" = "logfc_dot_product", 
                                           #                             "logFC Spearman" = "logfc_spearman", 
                                           #                             "logFC Pearson" = "logfc_pearson",
                                           #                          "Spearman (All Genes)" = "all_genes_spearman", 
                                           #                          "Pearson (All Genes)" = "all_genes_pearson"
                                           #                          ), 
                                           #              selected = "logFC Dot Product"), 
                                           br(),
                                           actionBttn(inputId = "runCIPR", label="Run CIPR", style = "unite", size = 'lg', color = "success"),
                                       
                                           
                                    ), # column
                                    
                                    br(),
                                    
                                    column(width = 7,
                                           includeMarkdown("descriptionfiles/helptext_CIPR.Rmd"),
                                           h4("Definitions:"),
                                           tags$ul( # ul is unordered and ol is ordered
                                             tags$li("Cluster: Cluster annotation provided by PlaqView."),
                                             tags$li("Reference: Broad classification of the reference cell type."), # change server code to toggle
                                             tags$li("Ref.ID: Shortened unique identifier for reference cell type."),
                                             tags$li("Full.Name: Human-readable long name of the reference cell type."),
                                             tags$li("Description: Additional details about this cell type."),
                                             tags$li("Identity.Score: Identity score for the given reference cell type calculated via logFC dot product or correlation methods (see CIPR documents for more details)."),
                                             tags$li("Percent.Pos.Cor.: The percentage of the differentially expressed genes in PlaqView clusters is also differentially expressed in a similar fashion in the reference cell subsets (e.g. both upregulated and downregulated)."),

                                           ),
                                    ), # column
                                    
                                    
                                    column(width = 12, 
                                           br(),
                                           plotOutput("CIPRplot", click = "brushtop5",
                                                      brush = brushOpts(id = "brushtop5")),
                                           # DT::dataTableOutput('CIPRtable'),
                                           disabled(downloadButton("download_top5", "Download CIPR Result Table")), 
                                           disabled(downloadButton("download_CIPRplot", "Download CIPR Plot")), 
                                           br(),
                                           DT::dataTableOutput("brushedtop5"),
                                           br(), br(), 
                                        
                                           
                                    ),
                                  ) # fluid row
                                ),
                      ), # mainPanel
                      
                      
                      
             ), # tabPanel
             #### UI: MetaData   ----
             tabPanel(title = "Metadata Explorer", value = "panelx",
                      mainPanel(width = 12, # 12/12 is full panel
                                fluidRow(
                                  ## panel for description
                                  column(width = 12,
                                    wellPanel(includeMarkdown("descriptionfiles/helptext_metadata.Rmd"))
                                  ), # column
                                  
                                  column(width = 12,
                                         wellPanel(
                                           fluidRow(
                                             column(width = 6, align = "center",
                                                    pickerInput("select.factor.variables", 
                                                                label = "Explore Factor-Type Variables",
                                                                choices = list (
                                                                  "Looks like you forgot to load a dataset"
                                                                ), 
                                                                selected = ""),
                                                    
                                                    plotOutput("plot.factor.variables",
                                                               height = '500px')),
                                            
                                             
                                             column(width = 6, align = "center",
                                                    pickerInput("select.continuous.variables", 
                                                                label = "Explore Continuous-Type Variables",
                                                                choices = list (
                                                                  "Looks like you forgot to load a dataset"
                                                                ), 
                                                                selected = ""),
                                                    
                                                    plotOutput("plot.continuous.variables",
                                                               height = '500px'),
                                                    
                                                    pickerInput("select.continuous.variables.dependency", 
                                                                label = "Change X Axis for Continuous Variable Plot",
                                                                choices = list (
                                                                  "Looks like you forgot to load a dataset"
                                                                ), 
                                                                selected = ""),
                                                    ),
                                             
                                             
                                             
                                           ), # fluidrow
                                           
                                         )# well panel
                                         ), # column
                                  
                                  
                                  column(width = 12,
                                    wellPanel(
                                      # must add up to 12 all columns
                                      fluidRow(
                                        column(width = 5,
                                               textInput(
                                                 "genes.metadata",
                                                 width = '100%',
                                                 h3("Query Gene Expression Split by Metadata", h5("please follow HUGO conventions")),
                                                 value = "APOE, COL1A1, FBLN1, FBLN2",
                                                 placeholder = "try: TREM2, CYBB"
                                               ),
                                               
                                               # choose the type of output graph
                                               pickerInput("selectaplot_metadata",
                                                           label = "Select Plot Type",
                                                           choices = list(
                                                             "Dot Plot (well-suited for multiple genes)" = "Dot",
                                                             "Ridge Plot (better visual for one gene)" = "Ridge"),
                                                           width = '95%',
                                                           selected = "Dot Plot"),
                                               
                                               pickerInput(
                                                 inputId = "selectlabelmethodforgenequery.metadata",
                                                 label = "Select Metadata Column (may be limited depending on dataset)",
                                                 choices = list (
                                                   # this is a reactive choice list based on whats available in the dataset
                                                   "Waiting for User to Load Data"
                                                 ),
                                                 # selected = "nCounts",
                                                 width = '95%' #neeed to fit this
                                               ),
                                               
                                               # 'go' button
                                               actionBttn(
                                                 inputId = "run.metadata",
                                                 label = "Generate Plot",
                                                 style = "unite",
                                                 color = "success",
                                                 size = "md",
                                                 block = F)
                                        ), # column (5)
                                        
                                        column(width = 7, align = "center",
                                               conditionalPanel(condition = "input.selectaplot_metadata == 'Dot'",
                                                                plotOutput("Dot.metadata", width = "auto", height = '500px'),
                                                                br(),
                                                                downloadButton("download.Dot.metadata", "Download This Plot", width = '100%') 
                                               ), # con. panel
                                               conditionalPanel(condition = "input.selectaplot_metadata == 'Ridge'",
                                                                plotOutput("Ridge.metadata", width = "auto", height = '500px'),
                                                                br(),
                                                                downloadButton("download.Ridge.metadata", "Download This Plot", width = '100%') 
                                               ), # con. panel
                                        ) # column (7)
                                        
                                      ),
                                      

                                    ) # wellpanel
                                  ), # column (12)
                                ), # fluid row


                                
                                fluidRow( # top split rows
                      
                                ), # fluidrow


                      )# MAIN PANEL CLOSURE
             ), # TAB PANEL CLOSURE




             #### UI: RNATraject ----
             tabPanel("Trajectory",
                      mainPanel(width = 12, # 12/12 is full panel,
                                ### top panel
                                wellPanel(
                                  fluidRow(
                                  column(width = 6,
                                         includeMarkdown("descriptionfiles/helptext_comparetrajectories.Rmd"),
                                         # choose button
                                         actionBttn("choose_toggle", "Select", style = "unite", color = "primary", block = F, size = "sm"),
                                         # clear button
                                         (actionBttn("reset", "Clear", style = "unite", color = "primary", block = F, size = "sm")),
                                         # recal button
                                         (actionBttn("redomonocle3", "Calculate", style = "unite", color = "success", block = F, size = "sm")),
                                         h4("Instructions:"),
                                         tags$ol(
                                           tags$li("Highlight points by clicking and dragging."),
                                           tags$li("Click the 'Select' button."), # change server code to toggle
                                           tags$li("Repeat until all of the desired cells are black."),
                                           tags$li("Click 'Calculate'.")
                                                ),
                                         h4("Details:"),
                                         tags$ul(
                                           tags$li("To start over, click 'Clear'"),
                                           tags$li(paste("You can also choose/unchoose specific cells",
                                                         "by clicking on them directly.")),
                                           tags$li(paste("You can un-select a cluster of cells",
                                                         "by clicking 'Choose' again.")),
                                           ),
                                         h5("Calculation may take up to 10mins. Results will appear below."),
                                         
                                         ), # column
                                       column(width = 6, 
                                              plotOutput("plot1", click = "plot1_click", 
                                                         brush = brushOpts(id = "plot1_brush")),
                                              br(),
                                              downloadButton("downloadoriginaltrajectory", 
                                                             label = "Original Trajectory", 
                                                             width = '50%'),
                                              disabled(downloadButton("downloadsubsettrajectory", 
                                                                      label = "Subset Trajectory",
                                                                      width = '50%')),

                                       ), # column
                                ), # fluidrow
                                ), # wellpanel
             
                                ### bottom panel
                                hidden(
                                  wellPanel(id = "subset.trajectory", 
                                            fluidRow(width = 12,
                                                     column(width = 6,
                                                            h4("Original Full Trajecotry"),
                                                            plotOutput("originaltrajectory"),
                                                     ),
                                                     column(width = 6,
                                                            h4("Trajectory of Selected Cells"),
                                                            plotOutput("subsettrajectory"),
                                                     ),
                                                     
                                                     
                                            ) # fluidrow
                                  ) # wellpanel
                                ),# hidden
                           
                      ), # mainPanel monocle3
                      
             ), # tabPanel
             
             #### UI: Drugs ----  
             tabPanel("Druggable Genome",
                      mainPanel(width = 12,
                                fluidRow(width = 12,
                                         column(width = 6,
                                                wellPanel(
                                                  includeMarkdown("descriptionfiles/helptext_druggablegenome.Rmd"),
                                                  textInput(
                                                    inputId = "druggeneinput",
                                                    label = "Gene to Drug",
                                                    value = "EGFR"
                                                  ),
                                                  
                                                  actionBttn(
                                                    inputId = "rundgidb",
                                                    label = "Start Query",
                                                    style = "unite",
                                                    color = "success",
                                                    block = T,
                                                    size = "lg"),
                                                  br(),
                                                  pickerInput("drugcelllabelmethod", 
                                                              label = "Change Labeling Method",
                                                              choices = list (
                                                                "Seurat_Clusters",
                                                                "Author_Provided",
                                                                "SingleR_calls" = "SingleR.calls",
                                                                "Seurat_with_Tabula_Ref"  
                                                              ), 
                                                              selected = "Seurat_with_Tabula_Ref"),
                                                  pickerInput(
                                                    inputId = "dgidbdatabase",
                                                    label = "Choose Database(s)", 
                                                    inline = TRUE, 
                                                    selected = c("COSMIC", "DrugBank", "FDA"), # preselect
                                                    choices = sourceDatabases(),
                                                    options = list(
                                                      `actions-box` = TRUE), 
                                                    multiple = TRUE, width = "100%"
                                                  ),
                                                  
                                                  helpText("You must restart query if you change database. PubMed ID and citations of interactions are available in full download file."),
                                                  
                                                ), # wellpanel
                                                
                                         ), 
                                         column(width = 6, 
                                                wellPanel(plotOutput("featurefordrugs",
                                                                     height = '500px'),
                                                          br(),
                                                          disabled(
                                                            downloadButton("downloadfeaturefordrugumap", 
                                                                           label = "Download this UMAP")) #disable
                                                )),
                                         br(),
                                         column(width = 12,
                                                DT::dataTableOutput("dgidboutput", width = "100%"),
                                                br(),
                                                disabled(downloadButton("downloadfulltable", label = "Download Full Gene-Drug Interaction Table")
                                                ),
                                                
                                         )
                                )
                      )
                      
             ),
             
             
             ### JS to jump to top of page on click ###
             tags$script(" $(document).ready(function () {
                         $('#jumpto1').on('click', function (e) {
                         window.scrollTo(0, 0)                
                         });
                         });"),
             tags$script(("$.fn.dataTable.ext.errMode = 'none';")),
             
             
             
  )# close navbarpage
  
)# close fluidpage



#### SERVER FUNCTIONS ####
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  #### SER: Data ####
  output$availabledatasettable <-
    DT::renderDataTable(df, server = F, # server is for speed/loading
                        selection = list(mode = 'single'),
                        # options=list(columnDefs = list(list(visible=FALSE, targets=c(10)))), # this hides the #8 col (datasetID)
                        options = list(pageLength = 20),
                        escape = FALSE) # this escapes rendering html (link) literally and makes link clickable
  
  # second of the same code.. may help resolve datatable not loading error
  output$availabledatasettable <-
    DT::renderDataTable(df, server = F, # server is for speed/loading
                        selection = list(mode = 'single'),
                        # options=list(columnDefs = list(list(visible=FALSE, targets=c(10)))), # this hides the #8 col (datasetID)
                        options = list(pageLength = 20),
                        escape = FALSE) # this escapes rendering html (link) literally and makes link clickable
  
  # refresh button
  observeEvent(input$refreshtable, {
    session$reload()
    
    output$availabledatasettable <-
      DT::renderDataTable(df, server = F, # server is for speed/loading
                          selection = list(mode = 'single'),
                          # options=list(columnDefs = list(list(visible=FALSE, targets=c(10)))), # this hides the #8 col (datasetID)
                          options = list(pageLength = 20),
                          escape = FALSE) # this escapes rendering html (link) literally and makes link clickable
    
  })
  
    # start the page with load data disabled until dataset is clicked
  # disable("loaddatabutton")
  # observeEvent(input$availabledatasettable_rows_selected,{
  #   enable("loaddatabutton")
  # })
  

  
   observeEvent(input$loaddatabutton, {
    # create path for loading data
    path <- file.path(paste("data/", input$dataselector, "/",
                            input$dataselector,
                            ".rds", sep=""))
    pathcds <- file.path(paste("data/", input$dataselector, "/",
                               input$dataselector,
                            "_cds.rds", sep=""))
  
    plaqviewobj <<- readRDS(file = path)
    plaqviewobj.cds <<- readRDS(file = pathcds)
    
    # show which data is read
    loadeddatasetID <<- paste("Dataset Loaded Sucessfully: ", print(input$dataselector))
    output$loadeddatasetID <- renderText(loadeddatasetID)
    
    ## these are just for displaying current data name in other tabs##
    output$selecteddatasetID <- renderText({
      paste0("Current dataset: ", input$dataselector)
    }) 
    output$selecteddatasetID2 <- renderText({
      paste0("Current dataset: ", input$dataselector)
    }) 
    output$selecteddatasetID3 <- renderText({
      paste0("Current dataset: ", input$dataselector)
    }) 
   

  })

  # server code to show jumpto1
  observeEvent(input$loaddatabutton, {
    
    shinyjs::show(id = "jumpto1")  
  }) 
  
  # server code to jump to page 1
  observeEvent(input$jumpto1, {
    updateTabsetPanel(session = getDefaultReactiveDomain(), "inTabset",
                      selected = "panel1") # this is to switch to tab1
    
  })

  
  #### SER: Genes ####
  # UMAP 
  observeEvent(input$runcode, {
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
  })
  observeEvent(input$runcode,{ 
    
    #### NOMENCLATURE UPDATE ####
    if(df$Species[df$DataID == input$dataselector] == "Human"){
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
      user_genes <<- str_split(input$genes, ", ")[[1]]

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
  
  #### dot plot ###
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
                  order = T,
                  pt.size = 1,
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
                            order = T,
                            pt.size = 1,
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
    
    ### GSEA #####
    
      parsed.genes <- str_split(input$genes, ", ")[[1]]
      enriched <- enrichr(genes = parsed.genes, 
                          database = enrichRdb) # this queries all of them
      cleanedenrichedtable <- select(enriched[[input$selectedenrichRdb]], -Old.Adjusted.P.value, -Old.P.value,)
      cleanedenrichedtable <- top_n(cleanedenrichedtable, 100) # top 100 will be rendered
      
      #select columns to display
      cleanedenrichedtable <- cleanedenrichedtable %>% select(Term, Overlap, Adjusted.P.value, Combined.Score, Genes)
      
      # force as.numeric to remove a bug in DT pkg
      #cleanedenrichedtable$Adjusted.P.value <- as.numeric(cleanedenrichedtable$Adjusted.P.value)
      #cleanedenrichedtable$Adjusted.P.value <- as.numeric(cleanedenrichedtable$Combined.Score)
      
      output$enrichtable <- DT::renderDataTable(cleanedenrichedtable, server = F)
      
      # Downloadable csv of selected dataset
      output$downloadenrichRdata <- downloadHandler(
        filename = function() {
          paste(input$genes, "_pathwayenrichment.csv", sep = "")
        },
        content = function(file) {
          write.csv(enriched[[input$selectedenrichRdb]], file, row.names = FALSE)
        } 
      )# close downloadhandler
    
    
    
  }) # observe event closure
  
  
  
  #### SER: Labels ####
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
      file.copy(paste("data/", input$dataselector, "/",
                      "diff_by_seurat.csv", sep = ""), file)
      
    }  )# close downloadhandler
  
  output$diffbyauthor <- downloadHandler(
    filename = "differential_markergenes_by_author_annotated_clusters.csv",
    content = function(file) {
      
      file.copy(paste("data/", input$dataselector, "/",
                      "diff_by_author.csv", sep = ""), file)
      
    }  )# close downloadhandler
  
  output$diffbysingleR <- downloadHandler(
    filename = "differential_markergenes_by_singleR_labels.csv",
    content = function(file) {
      file.copy(paste("data/", input$dataselector, "/",
                      "diff_by_singleR.csv", sep = ""), file)      
    }  )# close downloadhandler
  
  output$diffbyTS <- downloadHandler(
    filename = "differential_markergenes_by_tabulus_sapien_reference.csv",
    content = function(file) {
      file.copy(paste("data/", input$dataselector, "/",
                      "diff_by_Seurat_with_Tabula_Ref.csv", sep = ""), file)      
    }  )# close downloadhandler
  
  #### SER: CIPR ####
  observeEvent(input$runCIPR, {
    pathforCIPR <- file.path(paste("data/", input$dataselector, "/",
                                   input$CIPRoriginal, sep=""))
    ciprinput <- tryCatch(read.csv(file = pathforCIPR), error=function(e) NULL)
    
    # this line is for wei's trouble shooting only
    # ciprinput <- tryCatch(read.csv(file = "../DataProcessing/data/Wirka_2019/diff_by_singleR.csv"), error=function(e) NULL)
    
    if(is.null(ciprinput) == TRUE){ # if author did not provide annotation
      output$CIPRplot <- 
        renderPlot({
          plot.new()
          par(bg = "#f7f7f7")
          text(x = 0.5, y = 0.5, paste("This Dataset Does NOT Have Author-Provided Annotations \n",
                                       "Please Select Another Labeling Method"), 
               cex = 1.6, col = "black")
          
        })
      
    } else { # if its not empty
      
      ciprinput$gene <- str_to_lower(ciprinput$gene)
      
      # actual CIPR calculation
      CIPR <- CIPR(input_dat = ciprinput,
                   comp_method = "logfc_dot_product", 
                   reference = input$CIPRreference, 
                   plot_ind = F,
                   plot_top = T,
                   global_plot_obj = T,
                   global_results_obj = T,
                   top_num = 3)  
      
      
      # plot
      output$CIPRplot <- renderPlot({
        CIPR
      }) # renderplot
      
      # table for selected groups
      output$brushedtop5 <- renderDataTable({
        validate(
          need(input$brushtop5, "Select data points for detailed information")
        )
        brushedtable <- brushedPoints(CIPR_top_results, input$brushtop5)
        
        brushedtable.cleaned <<- brushedtable %>% 
          select(Cluster = cluster,
                 Reference = reference_cell_type, 
                 Ref.ID = reference_id,
                 Full.Name = long_name,
                 Description = description,
                 Identity.Score = identity_score,
                 Percent.Pos.Cor. = percent_pos_correlation
          )
        
        brushedtable.cleaned$Percent.Pos.Cor. <<- paste(round(brushedtable.cleaned$Percent.Pos.Cor., 1), "%")
        brushedtable.cleaned$Identity.Score <<- paste(round(brushedtable.cleaned$Identity.Score, 1), "")
        
        brushedtable.cleaned
      })
      
      output$download_top5<- downloadHandler(
        filename = function() {
          paste(input$dataselector, "_",
                input$CIPRoriginal, "_", input$CIPRreference, "_", input$CIPRmethod,
                "_result_tables.csv", sep = "")
        },
        content = function(file) {
          write.csv(brushedtable.cleaned, file, row.names = FALSE, col.names = T)
        }
        
      )# close downloadhandler
      
      output$download_CIPRplot<- downloadHandler(
        filename = function() {
          paste(input$dataselector, "_",
                input$CIPRoriginal, "_", input$CIPRreference, "_", input$CIPRmethod,
                "_CIPRplot.pdf", sep = "")
        },
        content = function(file) {
          pdf(file, width = 18, height = 5) # paper = defult is a4 size, a4r = a4 rotated
          plot(CIPR)
          dev.off()
        }
        
      )# close downloadhandler
      
      enable("download_top5")
      enable("download_CIPRplot")
    }
    
      
  }
  ) # obsevent runcipr
  


  
  #### SER: Metadata ####
  observeEvent(input$loaddatabutton, {
    
    # pull out names of the different classes within meta.data and combine character.factors
    # Cf is the factor type
    plaqview.metadata.character.type <<- names(plaqviewobj@meta.data %>% select_if(is.character))
    plaqview.metadata.factor.type <<- names(plaqviewobj@meta.data %>% select_if(is.factor))
    plaqview.metadata.cf <<- append(plaqview.metadata.character.type, plaqview.metadata.factor.type)
    
    # numeric is the continuous type
    plaqview.metadata.numeric.type <<- names(plaqviewobj@meta.data %>% select_if(is.numeric))
    
    updatePickerInput(session, "selectlabelmethodforgenequery.metadata",
                      label = "Select Unabridged Metadata Column (may be limited depending on dataset)",
                      choices = plaqview.metadata.cf,
                      selected = plaqview.metadata.cf[1])
    
    updatePickerInput(session, "select.factor.variables",
                      label = "Select Unabridged Factor-Type Variables (may be limited)",
                      choices = plaqview.metadata.cf,
                      selected = plaqview.metadata.cf[1])
    
    updatePickerInput(session, "select.continuous.variables",
                      label = "Select Unabridged Continuous-Type Variables (may be limited)",
                      choices = plaqview.metadata.numeric.type,
                      selected = plaqview.metadata.numeric.type[1])
    
    updatePickerInput(session, "select.continuous.variables.dependency",
                      label = "Change X Axis for Continuous-Type Variables",
                      choices = plaqview.metadata.cf,
                      selected = plaqview.metadata.cf[1])
  })
  #### factor variable plots ####
  output$plot.factor.variables <-
    renderPlot(
      DimPlot(
        plaqviewobj,
        reduction = "umap",
        label = TRUE,
        label.size = 5,
        repel = T,
        # repel labels
        pt.size = 1,
        cols = color_function(length(unique(plaqviewobj@meta.data[[input$select.factor.variables]]))),
        group.by = input$select.factor.variables ) + # group.by is important, use this to call metadata separation
        theme(legend.position="bottom", 
              legend.box = "vertical") +
        ggtitle("Metadata UMAP (Factor Variable)") +
        theme(plot.title = element_text(hjust =  0.5)) +
        guides(color = guide_legend(nrow = 5))
    )# render plot
  
  #### cont. variable plots ####
  output$plot.continuous.variables <-
    renderPlot(
      VlnPlot(plaqviewobj, features = input$select.continuous.variables,
              group.by = input$select.continuous.variables.dependency)
      
    )# render plot
  
  #### dot plot metadata ####
  observeEvent(input$run.metadata,{ # observe event puts a pause until pushed
    #### NOMENCLATURE UPDATE ###
    if(df$Species[df$DataID == input$dataselector] == "Human"){
      corrected <- str_to_upper(input$genes.metadata)
    } else{
      corrected <- str_to_title(input$genes.metadata)
    }
    
    updateTextInput(getDefaultReactiveDomain(),
                    "genes.metadata", value = corrected)
    
    # this is for the display
    output$Dot.metadata <- renderPlot({
      # parse string input 
      user_genes.metadata <- str_split(input$genes.metadata, ", ")[[1]]
      validate(need(input$selectaplot_metadata=="Dot", message=FALSE))
      DotPlot(plaqviewobj, 
              group.by = input$selectlabelmethodforgenequery.metadata,
              features = user_genes.metadata) + # a trick to sep long string input
        ggtitle("Expression Dot Plot") +
        theme(plot.title = element_text(hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })
    
    # this is for the download
    output$download.Dot.metadata<- downloadHandler(
      filename = function() {
        paste("dotplotmetadata.pdf", sep = "")
      },
      content = function(file) {
        pdf(file, width = 10, height = 6 ) # paper = defult is a4 size
        
        # parse string input 
        user_genes.metadata <- str_split(input$genes.metadata, ", ")[[1]]
        validate(need(input$selectaplot_metadata=="Dot", message=FALSE))
        temp <- DotPlot(plaqviewobj, 
                group.by = input$selectlabelmethodforgenequery.metadata,
                features = user_genes.metadata) + # a trick to sep long string input
          ggtitle("Expression Dot Plot") +
          theme(plot.title = element_text(hjust = 1)) +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        plot(temp) #this is all you need
        
        dev.off()
      }
      
    )# close downloadhandler
    

  }) # observe event closure
  
  
  
  
  #### ridge plot metadata ####
  observeEvent(input$run.metadata,{ # observe event puts a pause until pushed
    #### NOMENCLATURE UPDATE ###
    if(df$Species[df$DataID == input$dataselector] == "Human"){
      corrected <- str_to_upper(input$genes.metadata)
    } else{
      corrected <- str_to_title(input$genes.metadata)
    }
    
    updateTextInput(getDefaultReactiveDomain(),
                    "genes.metadata", value = corrected)
    
    # this is for the display
    output$Ridge.metadata <- renderPlot({
      # parse string input 
      user_genes.metadata <- str_split(input$genes.metadata, ", ")[[1]]
      validate(need(input$selectaplot_metadata=="Ridge", message=FALSE))
      RidgePlot(plaqviewobj, 
              group.by = input$selectlabelmethodforgenequery.metadata,
              features = user_genes.metadata) + # a trick to sep long string input
        ggtitle("Expression Ridge Plot") +
        theme(plot.title = element_text(hjust = 1)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    })
    
    # this is for the download
    output$download.Ridge.metadata<- downloadHandler(
      filename = function() {
        paste("ridgeplotmetadata.pdf", sep = "")
      },
      content = function(file) {
        
        pdf(file, width = 10, height = 6 )
        
        
        user_genes.metadata <- str_split(input$genes.metadata, ", ")[[1]]
        validate(need(input$selectaplot_metadata=="Ridge", message=FALSE))
        temp <- RidgePlot(plaqviewobj, 
                        group.by = input$selectlabelmethodforgenequery.metadata,
                        features = user_genes.metadata) + # a trick to sep long string input
          ggtitle("Expression Ridge Plot") +
          theme(plot.title = element_text(hjust = 1)) +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        plot(temp) #this is all you need
        
        dev.off()
      }
      
    )# close downloadhandler
    
    
  }) # observe event closure
  
  
  
  
  
  #### SER: RNATraject ####
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
  
  # this is for the download
  output$downloadoriginaltrajectory<- downloadHandler(
    filename = function() {
      paste(input$dataselector, "_original_trajectory.pdf", sep = "")
    },
    content = function(file) {
      pdf(file, paper = "default") # paper = defult is a4 size
    
      temp <- plot_cells(plaqviewobj.cds,
                 color_cells_by = "assigned_cell_type",
                 label_groups_by_cluster=F,
                 show_trajectory_graph = T,
                 trajectory_graph_segment_size = 1,
                 graph_label_size = 1, # size of # in circle
                 group_label_size = 4)
      plot(temp)
      dev.off()
    }
    
  )# close downloadhandler
  
  
  #### subset plot
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
    }, 
    
    )
    
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
    selectedcells <<- vals$keeprows
    
    
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
    
    subsetted <- preprocess_cds(subsetted, num_dim = 25) # we used 30 in earlier seurat scripts
    
    
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
    }) # observe event
  
  output$downloadsubsettrajectory<- downloadHandler(
    filename = function() {
      paste(input$dataselector, "_subset_trajectory.pdf", sep = "")
    },
    content = function(file) {
      pdf(file, paper = "default") # paper = defult is a4 size
      
      temp <- plot_cells(subsetted,
                         color_cells_by = "assigned_cell_type",
                         label_groups_by_cluster=F,
                         show_trajectory_graph = T,
                         trajectory_graph_segment_size = 1,
                         graph_label_size = 1, # size of # in circle
                         group_label_size = 4,
                         cell_size = 1,
                         alpha = 0.7,
                         scale_to_range = T)
      plot(temp)
      dev.off()
    }
    
  )# close downloadhandler
  
  
  #### show subset.trajectory 
    observeEvent(input$redomonocle3, {
      shinyjs::showElement(id= "subset.trajectory")
    })
    observeEvent(input$redomonocle3, {
      enable("downloadsubsettrajectory")
    })

  #### SER: Drugs ####
    observeEvent(input$rundgidb, {
      
      
      #### NOMENCLATURE UPDATE ###
      if(df$Species[df$DataID == input$dataselector] == "Human"){
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
      
      # show download buttons
      enable("downloadfeaturefordrugumap")
      enable("downloadfulltable")
      
      #fulltable$score <- as.numeric(fulltable$score) # bypass DT error
      
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
        fulltable <<- fulltable %>% mutate(pmids = map_chr(pmids, toString))
      }
      
      # plots
      output$featurefordrugs <- renderPlot({
        user_genes <<- str_split(corrected, ", ")[[1]]
        FeaturePlot(plaqviewobj, 
                    features = user_genes, label = T, repel = T
        ) + # a trick to sep long string input
          theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
          theme(plot.title = element_text(hjust = 1)) +
          theme(plot.title = element_text(hjust =  0.5)) 
      }) # render plot
      
      output$downloadfeaturefordrugumap<- downloadHandler(
        filename = function() {
          paste(input$dataselector, "-querydrug_umap.pdf", sep = "")
        },
        content = function(file) {
          pdf(file, paper = "default") # paper = defult is a4 size
          
          temp <- FeaturePlot(plaqviewobj,
                              order = T,
                              pt.size = 1,
                              features = user_genes, label = T, repel = T) + # a trick to sep long string input
            theme(legend.position="bottom", legend.box = "vertical") + # group.by is important, use this to call metadata separation
            theme(plot.title = element_text(hjust = 1)) +
            theme(plot.title = element_text(hjust =  0.5)) 
          plot(temp)
          dev.off()
        }
        
      )# close downloadhandler
      
      # render table
      output$dgidboutput <- DT::renderDataTable(isolatedtable,  server = F)
      
    
      output$downloadfulltable <- downloadHandler(
        filename = paste("dgidb_full_output.csv"),
        content = function(file) {
          
          temp <- as.data.frame(fulltable)
          
          write_csv(temp, file )
         # write_csv(temp, file = "Documents/test.csv" )
          
         # simplified table worked
         # write_excel_csv(isolatedtable, file)
          
          
        }  )# close downloadhandler
      
      
    })# observer event
    
  #### SER: About Functions #### 
  
  output$downloadsessioninfo <- downloadHandler(
    filename = paste(date(), "sesssioninfo.txt"),
    content = function(file) {
      write_lines(sessionInfo(), file)
    }  )# close downloadhandler
  
  #### SER: Waiter ####
  waiter_hide()
  
  

  
} # ends server function

# Run the application 
shinyApp(ui = ui, server = server)

