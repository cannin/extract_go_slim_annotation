library(readxl)
library(readr)
library(org.Hs.eg.db)
library(GO.db)
library(magrittr)
library(stringr)

source("read_gaf.R")

# PARAMETERS ---- 
max_go <- 3

## GAF file 
go_dir <- "."
gaf_file <- file.path(go_dir, "goa_human.goslim_generic_20200528.gaf")

## Gene Symbols
# NOTE: Replace "symbols" with project specific gene symbol vector
gene_list <- read_excel("Kurnit CCR_Supplemental Table 3_Meric-Bernstam_with title.xlsx", skip=2)
symbols <- gene_list$Gene

### Remove NAs caused by merged cells in Excel 
symbols <- symbols[!is.na(symbols)]

### Replace to symbols 
symbols <- str_remove(symbols, " \\(.*\\)$")

## Evidence 
keep_evidence_codes <- c("EXP", "IDA", "IEP", "IGI", "IMP", "IPI", "NAS", "TAS")
#drop_evidence_codes <- c("IBA", "IBD", "IEA", "ISS", "IC", "HTP", "HDA", "HMP", "HGI", "HEP" ,"ND")
drop_evidence_codes <- c()

### Strongest evidence: https://www.oncomedicmexico.com/wp-content/uploads/2017/04/The-druggable-genome-and-support-for-target-identification-and-validation-in-drug-development.pdf.pdf
tier1_evidence <- c("EXP", "IDA", "TAS") 
### Secondary evidence:  https://academic.oup.com/bioinformatics/article/34/17/i857/5093216
tier2_evidence <- c("IPI", "IMP", "IGI", "IEP", "IC") 
### Others: Largely ordered as on GO page: http://geneontology.org/docs/guide-go-evidence-codes/
tier3_evidence <- c("HTP", "HDA", "HMP", "HGI", "HEP") 
tier4_evidence <- c("IBA", "IBD", "IKR", "IRD", "ISS", "ISO", "ISA", "ISM", "IGC", "RCA") 
tier5_evidence <- c("NAS")
tier6_evidence <- c("ND", "IEA") 
evidence_order <- c(tier1_evidence, tier2_evidence, tier3_evidence, tier4_evidence, tier5_evidence, tier6_evidence)

## Specific GO IDs 
drop_go_ids <- c("molecular_function"="GO:0003674", 
                 "biological_process"="GO:0008150", 
                 "cellular_component"="GO:0005575", 
                 "cellular protein modification process"="GO:0006464", 
                 "response to stress"="GO:0006950",
                 "growth"="GO:0040007",
                 "enzyme binding"="GO:0019899", 
                 "biosynthetic process"="GO:0009058", 
                 "ion binding"="GO:0043167",
                 "cellular nitrogen compound metabolic process"="GO:0034641")

## Aspects (GO Ontologies)
drop_aspects <- c("C")
keep_aspects <- c("biological_process"="P")

keep_gaf_cols <- c("DB_Object_Symbol", "GO_ID", "Aspect", "Evidence_Code")

# READ DATA ----
gaf <- read_gaf(gaf_file, database = org.Hs.eg.db, 
                accession = "SYMBOL", filter.evidence = c())
gaf$Evidence_Code <- factor(gaf$Evidence_Code, levels=evidence_order)

## GO Data for Terms
ids <- keys(GO.db, keytype="GOID")
dat_go <- select(GO.db, keys=ids, columns=c("TERM","ONTOLOGY"), keytype="GOID")

## Filter GAF: Remove cellular compartment, root GO IDs in GO Slim, remove indirect annotations with qualifiers
#filtered_gaf <- gaf[!(gaf$Aspect %in% drop_aspects),]
filtered_gaf <- gaf[(gaf$Aspect %in% keep_aspects),]
filtered_gaf <- filtered_gaf[is.na(filtered_gaf$Qualifier),]
filtered_gaf <- filtered_gaf[!(filtered_gaf$GO_ID %in% drop_go_ids),]

## For debugging; see available annotations
tmp_gaf <- filtered_gaf[filtered_gaf$DB_Object_Symbol == "CCNE1", ]
tmp_gaf_org <- filtered_gaf[filtered_gaf$DB_Object_Symbol == "VHL", ]
table(tmp_gaf_org$GO_ID) %>% sort(., decreasing = TRUE)

## Limit columns
tmp_gaf <- filtered_gaf[, (colnames(filtered_gaf) %in% keep_gaf_cols)]

## Loop over symbols 
go_ids_df <- NULL 
go_ids_df_sm <- NULL 
for(i in 1:length(symbols)) {
  #i <- 1
  symbol <- symbols[i]
  
  tmp <- tmp_gaf[tmp_gaf$DB_Object_Symbol == symbol,]
  named_go <- merge(tmp, dat_go, by.x="GO_ID", by.y="GOID", all.x=TRUE)
  
  tmp <- table(named_go$GO_ID) %>% sort(., decreasing = TRUE)
  go_cnts <- data.frame(GO_ID=names(tmp), Evidence_Code_Count=as.vector(tmp))
  
  extracted_go <- merge(named_go, go_cnts, all.x=TRUE)
  if(nrow(extracted_go) > 0) {
    # Re-order
    extracted_go <- extracted_go[with(extracted_go, order(-Evidence_Code_Count, Evidence_Code)), ]
    
    # Minimal data.frame for taking top N
    min_df <- extracted_go[, !(colnames(extracted_go) %in% c("Evidence_Code", "Evidence_Code_Count"))] %>% unique
    
    # Keep the top N GO IDs by evidence
    if(nrow(min_df) >= max_go) {
      cur_go_ids <- min_df$GO_ID[1:max_go] 
    } else {
      cur_go_ids <- min_df$GO_ID[1:length(min_df)] 
    }
    
    tmp <- min_df[min_df$GO_ID %in% cur_go_ids, c("GO_ID", "DB_Object_Symbol", "TERM", "ONTOLOGY")]
    
    # Get a string 
    tmp_sm <- data.frame(symbol=symbol, annot=paste(tmp$TERM, collapse="|"))
  } else {
    tmp <- data.frame(GO_ID=NA, DB_Object_Symbol=symbol, TERM=NA, ONTOLOGY=NA)
    tmp_sm <- data.frame(symbol=symbol, annot=NA)
  }
  
  if(is.null(go_ids_df)) {
    go_ids_df <- tmp
  } else {
    go_ids_df <- rbind(go_ids_df, tmp)
  }
  
  if(is.null(go_ids_df_sm)) {
    go_ids_df_sm <- tmp_sm
  } else {
    go_ids_df_sm <- rbind(go_ids_df_sm, tmp_sm)
  }
}

# WRITE RESULTS ----
filename <- paste0(paste(names(keep_aspects), collapse = "_"), "_go_terms.txt")
write_tsv(go_ids_df_sm, filename)

# IGNORE ----
symbols %>% length
(symbols %in% filtered_gaf$DB_Object_Symbol) %>% which %>% length 

