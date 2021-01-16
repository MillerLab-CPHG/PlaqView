#### LIBRARIES #### 
library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(plotly)
library(shinybusy) #install.packages("shinybusy")
library()

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
  add_busy_bar(color = "#ff9142"), # THIS IS THE BUSY BAR
  
  
  # defining each 'tab' here
    navbarPage("sCADView", # title 
               
            # PANEL 1    ----
             tabPanel("Query Single Gene", 
                      mainPanel(width = 12, # 12/12 is full panel
                                fluidRow(## panel for gene input
                                    column(
                                      width = 4,
                                      wellPanel(                                      
                                        # must add up to 12 all columns
                                        textInput(
                                          "genes",
                                          width = '100%',
                                          h3("Search for a Gene", h5("please follow HUGO conventions")),
                                          placeholder = "try: MYH11"
                                        ),

                                        # choose the type of output graph 
                                        selectInput("select", label = h3("Dataset"), 
                                                    choices = list("Wirka scRNA (2019)" = "BLANK"), 
                                                    selected = "Wirka 2019"),                                        
                                        # 'go' button
                                        actionButton(
                                          inputId = "runcode",
                                          label = "Query Gene",
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
                          wellPanel(
                          fluidRow(
                                      column(width = 6, plotlyOutput("umaps")),
                                      column(width = 6, plotlyOutput("FeaturePlot"))
                                      ) # fluidrow
                          )# wellpanel
            
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
                    includeMarkdown("descriptionfiles/aboutusdescription.Rmd"),
                    br(),
                    img(src = "MSTPlogo.png", width = 233, height = 83)
                    
                  ))
                                )
                      
            )

  


#### SERVER ####
# Define server logic required to draw a histogram
server <- function(input, output) {
    # UMAP plot, interactive #
  observeEvent(input$runcode,{ 
    output$umaps <- renderPlotly(print(
        DimPlot(
            stanford,
            reduction = "umap",
            label = TRUE,
            label.size = 3,
            repel = T,
            # repel labels
            pt.size = 1,
            cols = manual_color_list,
            group.by = "SingleR.calls" # group.by is important, use this to call metadata separation
        ) + theme(legend.position="bottom", legend.box = "horizontal") +
          ggtitle("UMAP by Cell Type") +
          theme(plot.title = element_text(hjust = 0.5))
        
    ))
  })
    
   # Gene feature plot, interactive #
    observeEvent(input$runcode,{ # observe event puts a pause until pushed
      
        # feature plot
        output$FeaturePlot <- renderPlotly(print(FeaturePlot(stanford, 
                                                    features = input$genes) 
                                                   ))
        output$RidgePlot <- renderPlotly(print(RidgePlot(stanford, 
                                                              features = input$genes) 
        ))
        output$VlnPlot <- renderPlotly(print(VlnPlot(stanford, 
                                                              features = input$genes) 
        ))
        output$DotPlot <- renderPlotly(print(DotPlot(stanford, 
                                                              features = input$genes) 
        ))
        
        
 
        
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
