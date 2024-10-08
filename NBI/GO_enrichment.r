#!/usr/bin/Rscript
# Written making use of the informative guides at:
# http://avrilomics.blogspot.co.uk/2015/07/using-topgo-to-test-for-go-term.html
# http://rstudio-pubs-static.s3.amazonaws.com/6942_3214de070a324533a02c646bbbbfdb78.html


#get config options
library(optparse)
library(Rgraphviz)
library(topGO)
library(GO.db)
opt_list = list(
  make_option("--all_genes", type="character", help="directory containing fisher tables"),
  make_option("--GO_annotations", type="character", help="output results tab separated file"),
  make_option("--out_dir", type="character", help="output directory for the analysis")
  )

opt = parse_args(OptionParser(option_list=opt_list))
f1 = opt$all_genes
f2 = opt$GO_annotations
o = opt$out_dir

All_genes <- data.frame()
All_genes <- read.table(f1, header = FALSE, sep = "\t", stringsAsFactors = TRUE)

gene_enrichment3 <- t(All_genes[,2])
names(gene_enrichment3) <- All_genes[,1]

GO_relationships <- readMappings(file = f2)


# ---
# Test enrichment of Molecular function ontologies
# ----
GOdata <- new("topGOdata", ontology = "MF", allGenes = gene_enrichment3, geneSel = function(p) p < 0.05, description = "Test", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)
resultFisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
resultFisher
allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
out_prefix <- paste(o, "/", "TopGO_MF", sep="")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = out_prefix, useInfo = "all", pdfSW = TRUE)
length(usedGO(GOdata))

# ---
# Test enrichment of Biological process ontologies
# ----
GOdata <- new("topGOdata", ontology = "BP", allGenes = gene_enrichment3, geneSel = function(p) p < 0.05, description = "Test", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)
resultFisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
resultFisher
allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
out_prefix <- paste(o, "/", "TopGO_BP", sep="")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = out_prefix, useInfo = "all", pdfSW = TRUE)
length(usedGO(GOdata))

# ---
# Test enrichment of Cellular component ontologies
# ----
GOdata <- new("topGOdata", ontology = "CC", allGenes = gene_enrichment3, geneSel = function(p) p < 0.05, description = "Test", annot = annFUN.gene2GO, gene2GO = GO_relationships)
GOdata
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)
resultFisher <- runTest(GOdata, algorithm="weight01", statistic="fisher")
resultFisher
allRes <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo ='all')
out_prefix <- paste(o, "/", "TopGO_CC", sep="")
printGraph(GOdata, resultFisher, firstSigNodes = 10, fn.prefix = out_prefix, useInfo = "all", pdfSW = TRUE)
length(usedGO(GOdata))
