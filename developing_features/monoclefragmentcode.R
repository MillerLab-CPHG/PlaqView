#' Choose cells interactively to subset a cds
#'
#' @param cds CDS object to subset
#' @param reduction_method The reduction method to plot while choosing cells.
#' @param return_list Logical, return a list of cells instead of a subsetted
#'   CDS object.
#' @param clear_cds Logical, clear CDS slots before returning.
#'   After clearing the cds, re-run processing from preprocess_cds(), ...
#'   Default is FALSE.
#'
#' @return A subset CDS object. If return_list = FALSE, a list of cell names.
#' @export
#'
choose_cells <- function(cds,
                         reduction_method = c("UMAP", "tSNE", "PCA", "Aligned"),
                         clear_cds = FALSE,
                         return_list = FALSE) 
  
  
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste0("No dimensionality reduction for ",
                                       reduction_method, " calculated. ",
                                       "Please run reduce_dimension with ",
                                       "reduction_method = ", reduction_method,
                                       ", cluster_cells, and learn_graph ",
                                       "before running choose_cells"))
  assertthat::assert_that(is.logical(return_list))
  assertthat::assert_that(interactive(),
                          msg = paste("choose_cells only works in",
                                      "interactive mode."))
  
  
  reduced_dims <- as.data.frame(reducedDims(cds)[[reduction_method]])
  names(reduced_dims)[1:2] <- c("V1", "V2")
  
  ui <- shiny::fluidPage(
    shiny::titlePanel("Choose cells for a subset"),
    
    # Sidebar layout with input and output definitions ----
    shiny::sidebarLayout(
      
      # Sidebar panel for inputs ----
      shiny::sidebarPanel(
        # done button
        shiny::actionButton("choose_toggle", "Choose/unchoose"),
        # clear button
        shiny::actionButton("reset", "Clear"),
        # done button
        shiny::actionButton("done", "Done"),
        shiny::h3("Instructions:"),
        shiny::tags$ol(
          shiny::tags$li("Highlight points by clicking and dragging."),
          shiny::tags$li("Click the 'Choose/unchoose' button."),
          shiny::tags$li("Repeat until all of the desired cells are black."),
          shiny::tags$li("Click 'Done'.")
        ),
        shiny::h4("Details:"),
        shiny::tags$ul(
          shiny::tags$li("To start over, click 'Clear'"),
          shiny::tags$li(paste("You can also choose/unchoose specific cells",
                               "by clicking on them directly"))
        )
      ),
      
      # Main panel for displaying outputs ----
      shiny::mainPanel(
        shiny::plotOutput("plot1", height="auto",
                          click = "plot1_click",
                          brush = shiny::brushOpts(id = "plot1_brush"))
      )
    )
  )
  
  server <- function(input, output, session) {
    
    vals <- shiny::reactiveValues(
      keeprows = rep(FALSE, nrow(colData(cds)))
    )
    
    output$plot1 <- shiny::renderPlot({
      # Plot the kept and excluded points as two separate data sets
      colData(cds)$keep <- vals$keeprows
      
      suppressMessages(plot_cells(cds, reduction_method = reduction_method,
                                  cell_size = 1, label_cell_groups = FALSE,
                                  rasterize=FALSE) +
                         geom_point(alpha = colData(cds)$keep)) +
        theme(legend.position = "none")
    }, height = function() {
      session$clientData$output_plot1_width
    })
    
    # Toggle points that are clicked
    shiny::observeEvent(input$plot1_click, {
      res <- shiny::nearPoints(reduced_dims,
                               xvar = "V1", yvar = "V2", input$plot1_click,
                               allRows = TRUE)
      vals$keeprows <- vals$keeprows | res$selected_
    })
    
    # Toggle points that are brushed, when button is clicked
    shiny::observeEvent(input$choose_toggle, {
      res <- shiny::brushedPoints(reduced_dims,
                                  xvar = "V1", yvar = "V2", input$plot1_brush,
                                  allRows = TRUE)
      vals$keeprows <- vals$keeprows | res$selected_
    })
    
    # Reset all points
    shiny::observeEvent(input$reset, {
      vals$keeprows <- rep(FALSE, nrow(colData(cds)))
    })
    
    shiny::observeEvent(input$done, {
      shiny::stopApp(vals$keeprows)
    })
    
  }
  
  
  sel <- suppressMessages(shiny::runApp(shiny::shinyApp(ui, server)))
  if(return_list) {
    return(row.names(colData(cds)[sel,]))
  } else {
    if( clear_cds )
      return(clear_cds_slots(cds[,sel]))
    return(cds[,sel])
  }

