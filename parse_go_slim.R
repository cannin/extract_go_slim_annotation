# This file parses GO SLIMs from: https://geneontology.github.io/docs/go-subset-guide/
library(jsonlite)

# Parameters 
drop_go_ids <- c("response to stress"="GO:0006950", 
                 "mitotic cell cycle"="GO:0000278", 
                 "growth"="GO:0040007", 
                 "enzyme binding"="GO:0019899", 
                 "biosynthetic process"="GO:0009058", 
                 "ion binding"="GO:0043167", 
                 "cellular nitrogen compound metabolic process"="GO:0034641", 
                 "aging"="GO:0007568",
                 "anatomical structure development"="GO:0048856")

add_go_ids <- data.frame(
  go_child_term=c("GO:0032502"),
  go_child_term_name=c("developmental process"),
  ontology=c("BP"),
  stringsAsFactors = FALSE
)

# Extraction 
json <- fromJSON("goslim_generic_20200528.json")
t1 <- json$graphs$nodes[[1]]

idx <- grepl('GO_', t1$id)
t2 <- t1[idx,]

ontologies <- lapply(t2$meta$basicPropertyValues, function(x) {
  idx <- which(grepl('hasOBONamespace', x$pred))
  namespace <- x$val[idx]
})
ontologies <- unlist(ontologies)
ontologies <- sub('biological_process', 'BP', ontologies)
ontologies <- sub('cellular_component', 'CC', ontologies)
ontologies <- sub('molecular_function', 'MF', ontologies)

terms <- sub('http://purl.obolibrary.org/obo/GO_', 'GO:', t2$id)

# Write results
results <- data.frame(go_child_term=terms, go_child_term_name=t2$lbl, ontology=ontologies, stringsAsFactors=FALSE)
results <- results[order(results$ontology),]
write.table(results, "all_go_slim_generic_20200528.txt", sep="\t", row.names = FALSE, quote=FALSE)

## Drop/add any terms 
filtered_results <- results[!(results$go_child_term %in% drop_go_ids),]
filtered_results <- rbind(filtered_results, add_go_ids)

write.table(filtered_results, "filtered_go_slim_generic_20200528.txt", sep="\t", row.names = FALSE, quote=FALSE)

writeLines(filtered_results$go_child_term, "filtered_go_slim_generic_ids_only_20200528.txt")

