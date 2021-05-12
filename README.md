# Introduction 

Extract GO annotations for a set of gene symbols using evidence types and counts from a GO Slim subset 

# Needed Input Files 

Several files need to be download from GO: 

## Annotations
From: http://geneontology.org/docs/download-go-annotations/

* goa_human.gaf: http://current.geneontology.org/products/pages/downloads.html

## Ontology

From: https://geneontology.github.io/docs/download-ontology/

* go-basic.obo: http://purl.obolibrary.org/obo/go/go-basic.obo (used by OWLTools for mapping to GOSlim)
* goslim_generic.json: http://current.geneontology.org/ontology/subsets/goslim_generic.json 

# Code 
## Main Code
* parse_go_slim.R: Parse GO Slim terms from JSON (e.g., goslim_generic.json)
* extract_genes_from_go_gaf.R: Extract GO annotations for a set of gene symbols using evidence types and counts from a GO Slim subset; mapping of GO to GO Slim subset is done with OWLTools 

## Secondary Code
* read_gaf.R: Reads GAF file to data.frame; from: https://github.com/skinnider/flavin (used by extract_genes_from_go_gaf.R)

## OWLTools (For Generating GO Annotation File (GAF) Format GO Subsets)

The Gene Ontology (GO) can be mapped to specific terms using OWLTools https://github.com/owlcollab/owltools/wiki/Map2Slim

```
./owltools go-basic_20200528.obo --gaf goa_human_20200528.gaf --map2slim --idfile filtered_go_slim_generic_ids_only_20200528.txt --write-gaf goa_human.goslim_generic_20200528.gaf > /dev/null 2>&1
```

