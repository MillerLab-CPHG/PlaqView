reactable(df, compact = T, searchable = T, defaultPageSize = 20,
          defaultColDef = colDef(
            header = function(value) gsub(".", " ", value, fixed = TRUE),
            align = "center",
            minWidth  = 50,
            headerStyle = list(background = "#f7f7f8")
          ),
          columns = list(
            'Article.Title' = colDef(minWidth = 200,
                                     ),
            # DOI = colDef(html = TRUE, cell = function(value, index) {
            #   # this is raw html
            #   sprintf('<a href="%s" target="_blank">LINK</a>', df$DOI[index], value)
            #   }),
            Cells = colDef(format = colFormat(separators = TRUE),
                           footer = paste0("Total Cells: ", sum(as.numeric(df$Cells)))
                           ),
   
            Journal = colDef(html = TRUE, cell = function(value, index) {
              sprintf('<a href="%s" target="_blank">%s</a>', df$DOI[index], value)
            }),       
            Year = colDef(show = T),
            DOI = colDef(show = F),
  
            Tissue = colDef(show = F),
            Notes = colDef(show = F),
            Species = colDef(show = T),
            DataID = colDef(
              minWidth = 100,
              name = "Data ID",
              html = TRUE,
              cell = function(value, index){
                species <- df$Species[index]
                year <- df$Year[index]
                tissue <- df$Tissue[index]
                Notes <- df$Notes[index]
                DOI <- df$DOI[index]
                div(
                  div(style = list(fontWeight = 600), value),
                  #div(style = list(fontSize = 12), species),
                  #div(style = list(fontSize = 12), year),
                  div(style = list(fontSize = 12), tissue),
                  div(style = list(fontSize = 12), Notes)
                )# div
              }
            )
)
)

