loadPackages <- function() {
    if(!try(require("DESeq"))) stop("package affy is missing")
    if(!try(require("affy"))) stop("package affy is missing")
    if(!try(require("affyPLM"))) stop("package affyPLM is missing")
    if(!try(require("plier"))) stop("package plier is missing")
    if(!try(require("limma"))) stop("package limma is missing")
    if(!try(require("biomaRt"))) stop("package biomaRt is missing")
    if(!try(require("org.Dm.eg.db"))) stop("package org.Dm.eg.db is missing")
    if(!try(require("igraph"))) stop("igraph")
    if(!try(require("marray"))) stop("marray")
    if(!try(require("AnnotationDbi"))) stop("AnnotationDbi")
    if(!try(require("GenomicFeatures"))) stop("GenomicFeatures")


    source("./R code/loadData.R")
    source("./R code/extractFactors.R")
    source("./R code/runQC.R")
    source("./R code/diffExp.R")
    source("./R code/loadGSC.R")
    source("./R code/checkLoadArg.R")
    source("./R code/GSCstatBatch.R")
    source("./R code/runGSA.R")
    source("./R code/runGSAhyper.R")
    source("./R code/networkPlot.R")
    
    options(device = "windows")
}