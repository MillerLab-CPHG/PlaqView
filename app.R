#### LIBRARIES #### 
library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(plotly)

#### LOADING DATA ####
stanford <- readRDS(file = "data/final_stanford_labeled.rds")

#### SHINY OPTIONS, COLORs ####
# make the graphs match color of the UI
shinyOptions(plot.autocolors = TRUE)

# color definitions
manual_color_list <-
    c("rosybrown2",
      "palevioletred1",
      "cadetblue1",
      "lemonchiffon3",
      "darkseagreen",
      "skyblue3",
      "cadetblue3",
      "lemonchiffon4",
      "darkseagreen1",
      "darkseagreen2",
      "rosybrown3",
      "thistle2",
      "salmon1",
      "palevioletred3",
      "palevioletred4",
      "lightsteelblue3",
      "cadetblue2",
      "thistle3"
    )

#### UI ####
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # set theme
  theme = shinytheme("flatly"),
  
  # defining each 'tab' here
    navbarPage("sCADView", # title 
               
            # PANEL 1    
             tabPanel("Query Single Gene",
                      # Sidebar with a slider input for number of bins
                      sidebarLayout(
                        
                        # text box input
                        sidebarPanel(
                          textInput("genes", 
                                    h3("Search for a Gene", 
                                       h5("please follow HUGO conventions")),
                                    value = "ACTA2"), #default value
                          
                          # 'go' button
                          actionButton(
                            inputId = "runcode",
                            label = "Run",
                            width = '100'
                          ),
                          
                          helpText("Please be patient, it might take a while to query."), 
                          helpText("To compare multiple genes' expression pattern, use the 'Compare Genes' tab."),
                          helpText("Relevant graphs will appear below the original UMAP. 
                                   Hover over graphs for details or to zoom in/out.")
                        ),
                        
                        # Show a plot of the generated distribution
                        mainPanel(
                          
                          # plot of cell type UMAPs by gene expression
                          h2("UMAP by Cell Type"),
                          plotlyOutput("umaps", width = '95%', height = '40%'),
                          br(), # blank space
                          h2("Gene Expression Query"),
                          # plot of feature UMAPs
                          plotlyOutput("feature", width = '95%', height = '40%')
                        )
                      )
             
                      
                      ),
            
            # PANEL 2    
            tabPanel("Compare Genes",
                         helpText("Support to compare genes will be coming soon")

                     ),
                     
            # PANEL 3    
            tabPanel("Compare Trajectory Methods",
                     helpText("Support to compare different trajectory methods will be coming soon")
                     
            ),
            
            # ABOUT PANEL
             tabPanel("About",
                  mainPanel(
                    # descriptions
                    includeMarkdown("description.Rmd")
                  )
                                )
                      
            )
  )
  


#### SERVER ####
# Define server logic required to draw a histogram
server <- function(input, output) {
    # UMAP plot, interactive #
    output$umaps <- renderPlotly(print(
        DimPlot(
            stanford,
            reduction = "umap",
            label = TRUE,
            label.size = 6,
            repel = T,
            # repel labels
            pt.size = 1,
            cols = manual_color_list,
            group.by = "SingleR.calls" # group.by is important, use this to call metadata separation
        )
    ))
    
   # Gene feature plot, interactive #
    observeEvent(input$runcode,{ # observe event puts a pause until pushed
        output$feature <- renderPlotly(print(FeaturePlot(stanford, 
                                                    features = input$genes)
                                                    ))
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
