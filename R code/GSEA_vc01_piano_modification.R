##############
# unfortunatelly, the original Piano is not fully suitable for my case
# so I re-used Piano scripts making them less generic yet fitting my data better
rm(list=ls())

#############
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
    if(!try(require("gplots"))) stop("gplots")

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

#############
setwd("C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertibrates/GSEA/GSEA_vc01")
loadPackages()

#############
l_da <- loadData(datadir = getwd(), countfile = "/Data/all.txt", setup="/Data/setup.txt", verbose=TRUE)


extractFactors(l_da)

#############
runQC(l_da, hist=FALSE, boxplot=FALSE, pca=TRUE, verbose=TRUE)

#############

# rm(diffExp)
# source("./R code/diffExp.R")

# dexp <- diffExp(l_da, 
#                 contrasts=c("Irradiation_fm - Control_fm", 
#                             "Dioxin_fm - Control_fm",
#                             "toluene_fm - Control_fm"),
#                 fitMethod="robust", 
#                 adjustMethod="BY", significance=0.05,
#                 plot=TRUE, heatmapCutoff=1e-10, volcanoFC=2
#                 )


rm("diffExp")
source("./R code/diffExp.R")
dexp <- diffExp(l_da, 
                 contrasts=c("Irradiation_fm - Control_fm", 
                             "toluene_fm - Control_fm",
                             "Dioxin_fm - Control_fm",
                             "formaldehyde_fm - Control_fm",
                             "Irradiation_m - Control_m", 
                             "toluene_m - Control_m",
                             "Dioxin_m - Control_m",
                             "formaldehyde_m - Control_m"),
                 fitMethod="nbinom", 
                 adjustMethod="BY", significance=0.05,
                 plot=FALSE, heatmapCutoff=1e-10, volcanoFC=2
)


#############
rm("loadGSC")
source("./R code/loadGSC.R")
# "biological_process" "molecular_function" "cellular_component" ""  
goSubgroup <- "cellular_component"
myGsc <- loadGSC(dexp, group=goSubgroup)




#############
# source("./R code/runGSAhyper.R")

contrasts=c("Irradiation_fm - Control_fm", 
            "toluene_fm - Control_fm",
            "Dioxin_fm - Control_fm",
            "formaldehyde_fm - Control_fm",
            "Irradiation_m - Control_m", 
            "toluene_m - Control_m",
            "Dioxin_m - Control_m",
            "formaldehyde_m - Control_m")

GSAsign <- 0.05

for (contrast in contrasts) {
    if (file.exists(paste("./Data/GSA_hyper_",contrast,"_",goSubgroup,"_",GSAsign,".rds",sep=""))) {
        gsah<- readRDS(paste("./Data/GSA_hyper_",contrast,"_",goSubgroup,"_",GSAsign,".rds",sep=""))
    }
    else {
        
        dt <- dexp$resTable[[contrast]][complete.cases(dexp$resTable[[contrast]]),]
        gsah <- runGSAhyper(genes = dt[,1], 
                            pvalues = dt[,5], 
                            pcutoff = GSAsign, 
                            gsc = myGsc, 
                            gsSizeLim=c(1,Inf), 
                            adjMethod="fdr")
        saveRDS(gsah,paste("./Data/GSA_hyper_",contrast,"_",goSubgroup,"_",GSAsign,".rds",sep=""))
    
    }
}


#############
networkPlot(gsah,class="non",significance=0.05)
