test <- queryDGIdb("TP53", sourceDatabases = sourceDatabases())
fulltable <- test@data[["interactions"]][[1]]
