parsed.genes <- str_split(input$genes, ", ")[[1]]
enriched <- enrichr(genes = parsed.genes, 
                    database = enrichRdb) # this queries all of them

cleanedenrichedtable <- select(enriched[["ChEA_2016"]], -Old.Adjusted.P.value, -Old.P.value,)
cleanedenrichedtable <- top_n(cleanedenrichedtable, 100) # top 100 will be rendered

#select columns to display
cleanedenrichedtable <- cleanedenrichedtable %>% select(Term, Overlap, Adjusted.P.value, Combined.Score, Genes)

# force as.numeric to remove a bug in DT pkg
cleanedenrichedtable$Adjusted.P.value <- as.numeric(cleanedenrichedtable$Adjusted.P.value)
cleanedenrichedtable$Adjusted.P.value <- as.numeric(cleanedenrichedtable$Combined.Score)

output$enrichtable <- DT::renderDataTable(cleanedenrichedtable)
