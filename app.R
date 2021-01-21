#### LIBRARIES #### 
library(BiocManager)
library(shiny)
library(shinythemes)
library(Seurat)
library(shinybusy) #install.packages("shinybusy")
library(tidyverse)


#### LOADING DATA ####
# below line is commented for shinyapp.io deployment temp
stanford <- readRDS(file = "data/final_stanford_labeled.rds")
# stanford <- readRDS(file = url("https://virginia.box.com/shared/static/oyo1bicpvlxen940zmciqapvg0y3n6gb.rds"))

#### SHINY OPTIONS, COLORs ####
# make the graphs match color of the UI
shinyOptions(plot.autocolors = TRUE)


# color definitions
manual_color_list <-
    c("rosybrown2",
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
    )

#### UI ####
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # set theme
  theme = shinytheme("flatly"),
  add_busy_bar(color = "#ff9142"), # THIS IS THE BUSY BAR
  
  
  # defining each 'tab' here
    navbarPage("sCADView", # title 
               
            # PANEL 1: QUERY 1 GENE   ----
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
                                          placeholder = "try: CDH2"
                                        ),

                                        # choose the type of output graph 
                                        selectInput("selectaplot", label = h5("Plot Type"), 
                                                    choices = list("Feature Plot" = "Feature",
                                                                   "Dot Plot" = "Dot",
                                                                   "Ridge Plot" = "Ridge"), 
                                                    selected = "Feature Plot"),                                        
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
                                      column(width = 6, plotOutput("umaps")),
                                      column(width = 6, 
                                             conditionalPanel('input.selectaplot=="Ridge"', plotOutput("Ridge")),
                                             conditionalPanel('input.selectaplot=="Dot"', plotOutput("Dot")), # conditonal panels renders only if conditions are met
                                             conditionalPanel('input.selectaplot=="Feature"', plotOutput("Feature")),
                                             conditionalPanel('input.selectaplot=="test"', plotOutput("test"))
                                             
                                      ) # column 
                                             
                                      ) # fluidrow
                          )# wellpanel
            
                        )
                      ),
            
                      
            
            # PANEL 2: MULTIPLE GENES ----   
            tabPanel("Query Multiple Genes",
                         helpText("Support coming soon!")

                     ),
                     
            # PANEL 3: TRAJECTORY MODELS  ----  
            tabPanel("Compare Trajectory Methods",
                     helpText("Support coming soon!")
                     
            ),
            
            # PANEL 4: EXPLORE YOUR OWN DATA ----  
            tabPanel("Upload a dataset",
                     helpText("Support coming soon!")
                     
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
            theme(plot.title = element_text(hjust = 1)) +
            guides(color = guide_legend(nrow = 5))
      ) # closes renderPlot
  })# closes observe event
    
   # Gene feature plot, interactive #
  observeEvent(input$runcode,{ # observe event puts a pause until pushed
      
      output$Dot <- renderPlot({
        validate(need(input$selectaplot=="Dot", message=FALSE))
        print(DotPlot(stanford, 
                          features = input$genes) + # group.by is important, use this to call metadata separation
                ggtitle(paste(input$genes, "Expression Dot Plot")) +
                theme(plot.title = element_text(hjust = 1)) 
        )
      })
      
      output$Feature <- renderPlot({
        validate(need(input$selectaplot=="Feature", message=FALSE))
        print(FeaturePlot(stanford, 
                      features = input$genes)  + # group.by is important, use this to call metadata separation
                theme(legend.position="bottom", legend.box = "horizontal") + # group.by is important, use this to call metadata separation
                ggtitle(paste(input$genes, "Expression Feature Plot")) +
                theme(plot.title = element_text(hjust = 1)) 
        )
      }) 
    
      output$Ridge <- renderPlot({
        validate(need(input$selectaplot=="Ridge", message=FALSE))
        RidgePlot(stanford, 
                      features = input$genes)  + # group.by is important, use this to call metadata separation
          theme(legend.position="bottom", legend.box = "horizontal") + # group.by is important, use this to call metadata separation
          ggtitle(paste(input$genes, "Expression Ridge Plot")) +
          theme(plot.title = element_text(hjust = 1)) 
      })
      
      ### download datasets 
      
    })

    #### PANEL #2 FUNCTIONS ####
    #### PANEL #3 FUNCTIONS ####

}

# Run the application 
shinyApp(ui = ui, server = server)
