# Error in queryAPIget(url) : Oops, badly formatted query. 


updated_druggeneinput <- "TP53"


result <- queryDGIdb("TP53",
                     sourceDatabases = "KEGG_2019_Human")
fulltable <- result@data[["interactions"]][[1]]


library(httr)
library(jsonlite)

result = GET("https://dgidb.org/api/v2/interactions.json?genes=FLT1,MM1,FAKE&nteraction_sources=TALC")


