library(shiny)
library(DT)

shinyApp(
  ui <- fluidPage(
    DT::dataTableOutput("data"),
    textOutput('myText')
  ),
  
  server <- function(input, output) {
    
    myValue <- reactiveValues(employee = '')
    
    shinyInput <- function(FUN, len, id, ...) {
      inputs <- character(len)
      for (i in seq_len(len)) {
        inputs[i] <- as.character(FUN(paste0(id, i), ...))
      }
      inputs
    }
    
    df <- reactiveValues(data = data.frame(
      
      Name = c('Dilbert', 'Alice', 'Wally', 'Ashok', 'Dogbert'),
      Motivation = c(62, 73, 3, 99, 52),
      Actions = shinyInput(actionButton, 5, 'button_', label = "Fire", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)' ),
      stringsAsFactors = FALSE,
      row.names = 1:5
    ))
    
    
    output$data <- DT::renderDataTable(
      df$data, server = FALSE, escape = FALSE, selection = 'none'
    )
    
    observeEvent(input$select_button, {
      selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
      myValue$employee <<- paste('click on ',df$data[selectedRow,1])
    })
    
    
    output$myText <- renderText({
      
      myValue$employee
      
    })
    
  }
)