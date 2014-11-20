setInternet2(use=T)
library(httr)
use_proxy("172.19.1.1", port = 3128, username = NULL,
          password = NULL)

citation("DSS")

setwd("C:/Users/Sony/Documents/")
#suppressPackageStartupMessages(require(clusterProfiler))
#Gff2GeneTable("dmel-all-no-analysis-r5.51.gff")


all_names<-read.delim("C:/Users/Sony/Google Äèñê/ÔÖÏ/RNASeq/Sample_CMR2.reads", header = F)
names(all_names)
eg <- all_names[,1]
#eg <- unique(eg[!is.na(eg)])


rm(all2, IMR_up, IMR_down, IFR_up, IFR_down, DMR_up, DMR_down, DFR_up, DFR_down, FMR_up, FMR_down, FFR_up, FFR_down, TMR_up, TMR_down, TFR_up, TFR_down)
all2<-read.table("C:/Users/Sony/Google Äèñê/ÔÖÏ/ñòàòüÿ/dss_new_tr.csv", header = T, sep=";", stringsAsFactors=F, fill = F, na.strings = 'NA', blank.lines.skip = T)
cat(names(all2))
attach(all2)
IMR_up 
IMR_down 
IFR_up
IFR_down
DMR_up 
DMR_down 
DFR_up 
DFR_down 
FMR_up 
FMR_down 
FFR_up 
FFR_down
TMR_up
TMR_down
TFR_up
TFR_down

IMR_up<-IMR_up[!is.na(IMR_up)]
IMR_down<-IMR_down[!is.na(IMR_down)]
IFR_up<-IFR_up[!is.na(IFR_up)]
IFR_down<-IFR_down[!is.na(IFR_down)]
DMR_up<-DMR_up[!is.na(DMR_up)]
DMR_down<-DMR_down[!is.na(DMR_down)]
DFR_up<-DFR_up[!is.na(DFR_up)]
DFR_down<-DFR_down[!is.na(DFR_down)]
FMR_up<-FMR_up[!is.na(FMR_up)]
FMR_down<-FMR_down[!is.na(FMR_down)]
FFR_up<-FFR_up[!is.na(FFR_up)]
FFR_down<-FFR_down[!is.na(FFR_down)]
TMR_up<-TMR_up[!is.na(TMR_up)]
TMR_down<-TMR_down[!is.na(TMR_down)]
TFR_up<-TFR_up[!is.na(TFR_up)]
TFR_down<-TFR_down[!is.na(TFR_down)]

IMR_up<-unique(IMR_up[IMR_up != ""])
IMR_down<-unique(IMR_down[IMR_down != ""])
IFR_up<-unique(IFR_up[IFR_up != ""])
IFR_down<-unique(IFR_down[IFR_down != ""])
DMR_up<-unique(DMR_up[DMR_up != ""])
DMR_down<-unique(DMR_down[DMR_down != ""])
DFR_up<-unique(DFR_up[DFR_up != ""])
DFR_down<-unique(DFR_down[DFR_down != ""])
FMR_up<-unique(FMR_up[FMR_up != ""])
FMR_down<-unique(FMR_down[FMR_down != ""])
FFR_up<-unique(FFR_up[FFR_up != ""])
FFR_down<-unique(FFR_down[FFR_down != ""])
TMR_up<-unique(TMR_up[TMR_up != ""])
TMR_down<-unique(TMR_down[TMR_down != ""])
TFR_up<-unique(TFR_up[TFR_up != ""])
TFR_down<-unique(TFR_down[TFR_down != ""])


#require(biomaRt)
#drosophila = useMart("ensembl")
#dme = useDataset("dmelanogaster_gene_ensembl", mart = drosophila)
#gomap <- getBM(attributes = c("entrezgene", "go_id"), filters = "flybase_transcript_id",
#               values = as.character(eg), mart = dme)
#dim(gomap)
#head(gomap)
#buildGOmap(gomap)

gs_f<-c("IMR_up", "IMR_down", "IFR_up", "IFR_down", "DMR_up", "DMR_down", "DFR_up", "DFR_down", "FMR_up", "FMR_down", "FFR_up", "FFR_down", "TMR_up", "TMR_down", "TFR_up", "TFR_down")
gs_id<-c("IMR_up_id", "IMR_down_id", "IFR_up_id", "IFR_down_id", "DMR_up_id", "DMR_down_id", "DFR_up_id", "DFR_down_id", "FMR_up_id", "FMR_down_id", "FFR_up_id", "FFR_down_id", "TMR_up_id", "TMR_down_id", "TFR_up_id", "TFR_down_id")

gs<-IMR_up
gs<-IMR_down 
gs<-IFR_up
gs<-IFR_down
gs<-DMR_up 
gs<-DMR_down 
gs<-DFR_up 
gs<-DFR_down 
gs<-FMR_up 
gs<-FMR_down 
gs<-FFR_up 
gs<-FFR_down
gs<-TMR_up
gs<-TMR_down
gs<-TFR_up
gs<-TFR_down

require(biomaRt)
dros = useMart("ensembl")
dm = useDataset("dmelanogaster_gene_ensembl", mart = dros)


map_gs <- getBM(attributes = c("entrezgene", "go_id"), filters = "flybase_transcript_id",
               values = as.character(gs), mart = dm)
dim(map_gs)
head(map_gs)
entrezgene<-map_gs$entrezgene
#names<-gomap$ensembl_gene_id
#names_go<-gomap$go_id
#names_go_name<-gomap$definition_1006

IMR_up_id<-entrezgene
IMR_down_id<-entrezgene
IFR_up_id<-entrezgene
IFR_down_id<-entrezgene
DMR_up_id<-entrezgene 
DMR_down_id<-entrezgene 
DFR_up_id<-entrezgene 
DFR_down_id<-entrezgene 
FMR_up_id<-entrezgene 
FMR_down_id<-entrezgene 
FFR_up_id<-entrezgene 
FFR_down_id<-entrezgene
TMR_up_id<-entrezgene
TMR_down_id<-entrezgene
TFR_up_id<-entrezgene
TFR_down_id<-entrezgene


Reduce(intersect,  list(v1=IMR_up_id,v2=DMR_up_id,v3=FMR_up_id,v4=TMR_up_id))
cat(Reduce(intersect,  list(v1=IMR_up_id,v2=DMR_up_id,v4=TMR_up_id)))
Reduce(intersect,  list(v2=DMR_up_id,v3=FMR_up_id,v4=TMR_up_id))

Reduce(intersect,  list(v1=IFR_up_id,v2=DFR_up_id,v3=FFR_up_id,v4=TFR_up_id))
cat(Reduce(intersect,  list(v1=IFR_up_id,v2=DFR_up_id,v4=TFR_up_id)))

Reduce(intersect,  list(v1=IMR_down_id,v2=DMR_down_id,v3=FMR_down_id,v4=TMR_down_id))
Reduce(intersect,  list(v1=IMR_down_id,v3=FMR_down_id,v4=TMR_down_id))
Reduce(intersect,  list(v1=IMR_down_id,v2=DMR_down_id,v3=FMR_down_id))
cat(Reduce(intersect,  list(v1=IMR_down_id,v2=DMR_down_id,v4=TMR_down_id)))
Reduce(intersect,  list(v2=DMR_down_id,v3=FMR_down_id,v4=TMR_down_id))

Reduce(intersect,  list(v1=IFR_down_id,v2=DFR_down_id,v3=FFR_down_id,v4=TFR_down_id))
cat(Reduce(intersect,  list(v1=IFR_down_id,v2=DFR_down_id,v4=TFR_down_id)))
Reduce(intersect,  list(v2=DFR_down_id,v3=FFR_down_id,v4=TFR_down_id))


###ClusterProfiler example
install.packages(c("GOSemSim", "DOSE", "clusterProfiler"),repos="http://www.bioconductor.org/packages/devel/bioc")
library("clusterProfiler")
require(org.Dm.eg.db)
gene=as.list(org.Dm.egGO2EG)[[5]]
gene
x=enrichGO(gene, ont="BP", organism="fly", readable=TRUE)
head(summary(x))

entrezg<-IMR_up_id
x_Rad_MF <- groupGO(gene=as.character(entrezg),
                 organism="fly",
                 ont="MF",
                 level=2,
                 readable=TRUE)
par(ps = 8, cex = 1, cex.main = 1)
plot(x_Rad_MF)

x_Rad_BP <- groupGO(gene=as.character(entrezg),
                 organism="fly",
                 ont="BP",
                 level=2,
                 readable=TRUE)
par(ps = 8, cex = 1, cex.main = 1)
plot(x_Rad_BP)


y_Rad <- enrichGO(gene=as.character(entrezg),
                  organism="fly",
                  ont="MF",
                  pvalueCutoff=0.01,
                  qvalueCutoff=0.05,
                  readable=TRUE)
par(ps = 8, cex = 1, cex.main = 1)
plot(y_Rad)

z_Rad <- enrichKEGG(gene=as.character(entrezg),
                    organism="fly",
                    pvalueCutoff=0.05,
                    qvalueCutoff=0.05,
                    readable=TRUE)
par(ps = 8, cex = 1, cex.main = 1)
plot(z_Rad)
plot(y_Rad, type="cnet", categorySize="geneNum", output="interactive")

Sample<-list(IMR_up=IMR_up_id, DMR_up=DMR_up_id, FMR_up=FMR_up_id, TMR_up=TMR_up_id)
Sample<-list(IMR_up=IMR_up_id, DMR_up=DMR_up_id, FMR_up=FMR_up_id, TMR_up=TMR_up_id)

Sample<-list(IMR_up=IMR_up_id, DMR_up=DMR_up_id, FMR_up=FMR_up_id, TMR_up=TMR_up_id, IFR_up=IFR_up_id, DFR_up=DFR_up_id, FFR_up=FFR_up_id, TFR_up=TFR_up_id)
Sample<-list(IMR_down=IMR_down_id, DMR_down=DMR_down_id, FMR_down=FMR_down_id, TMR_down=TMR_down_id, IFR_down=IFR_down_id, DFR_down=DFR_down_id, FFR_down=FFR_down_id, TFR_down=TFR_down_id)

IMR_up_id<-IMR_up_id[!is.na(IMR_up_id)]
DMR_up_id<-DMR_up_id[!is.na(DMR_up_id)]
FMR_up_id<-FMR_up_id[!is.na(FMR_up_id)]
TMR_up_id<-TMR_up_id[!is.na(TMR_up_id)]
IFR_up_id<-IFR_up_id[!is.na(IFR_up_id)]
DFR_up_id<-DFR_up_id[!is.na(DFR_up_id)]
FFR_up_id<-FFR_up_id[!is.na(FFR_up_id)]
TFR_up_id<-TFR_up_id[!is.na(TFR_up_id)]

IMR_down_id<-IMR_down_id[!is.na(IMR_down_id)]
DMR_down_id<-DMR_down_id[!is.na(DMR_down_id)]
FMR_down_id<-FMR_down_id[!is.na(FMR_down_id)]
TMR_down_id<-TMR_down_id[!is.na(TMR_down_id)]
IFR_down_id<-IFR_down_id[!is.na(IFR_down_id)]
DFR_down_id<-DFR_down_id[!is.na(DFR_down_id)]
FFR_down_id<-FFR_down_id[!is.na(FFR_down_id)]
TFR_down_id<-TFR_down_id[!is.na(TFR_down_id)]


Sample<-list(IMR_up=as.character(IMR_up_id), DMR_up=as.character(DMR_up_id), FMR_up=as.character(FMR_up_id), TMR_up=as.character(TMR_up_id), IFR_up=as.character(IFR_up_id), DFR_up=as.character(DFR_up_id), FFR_up=as.character(FFR_up_id), TFR_up=as.character(TFR_up_id))
Sample<-list(IMR_down=as.character(IMR_down_id), DMR_down=as.character(DMR_down_id), FMR_down=as.character(FMR_down_id), TMR_down=as.character(TMR_down_id), IFR_down=as.character(IFR_down_id), DFR_down=as.character(DFR_down_id), FFR_down=as.character(FFR_down_id), TFR_down=as.character(TFR_down_id))


names(Sample)

xGGo <- compareCluster(Sample,
                       fun="groupGO",
                       ont="BP",
                       organism="fly"
)
plot(xGGo)
plot(xGGo, type="bar", by="percentage")
plot(xGGo, type="bar", by="count")

xeGO <- compareCluster(Sample,
                       fun="enrichGO",
                       organism="fly",
                       pvalueCutoff=0.05,
                       qvalueCutoff=0.05)
plot(xeGO)
plot(xeGO, type="bar", by="percentage")
plot(xeGO, type="bar", by="count")

xP <- compareCluster(Sample,
                     fun="enrichPathway",
                     organism="fly",
                     pvalueCutoff=0.05,
                     qvalueCutoff=0.05)
plot(xP)
plot(xP, type="bar", by="percentage")
plot(xP, type="bar", by="count")


xGo <- compareCluster(Sample,
                      fun="enrichGO",
                      ont="BP",
                      organism="fly",
                      pvalueCutoff=0.05,
                      qvalueCutoff=0.05)
plot(xGo)
plot(xGo, type="bar", by="percentage")
plot(xGo, type="bar", by="count")

xKG <- compareCluster(Sample,
                      fun="enrichKEGG",
                      organism="fly",
                      pvalueCutoff=0.05)

plot(xKG)
plot(xKG, type="bar", by="percentage")
plot(xKG, type="bar", by="count")






#############################
library("GeneAnswers")
library("org.Dm.eg.db")
library("GO.db")
fb.entrez <- unlist(as.list(org.Dm.egALIAS2EG))
genes1<-IMR_up
iv <- match(genes1, names(fb.entrez))
iv <- iv[!is.na(iv)]

x <- cbind(x,rep(NA, nrow(x)))
x[,3] <- fb.entrez[iv]
topset <- x[myTopHits,3]
topset <- topset[!is.na(topset)]
foo <- geneAnswersBuilder(topset, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG')
go.bp <- foo@enrichmentInfo
###########################
#genes<-c("33156")
#mf1 <- enrichGO(genes, ont = "MF", organism = "fly")
#mf2 <- enrichGO(as.character(entrezgene), ont = "MF", organism = "fly", pvalueCutoff = 0.05, qvalueCutoff = 0.1,
#readable = TRUE)

#summary(mf1)
#summary(mf2)

#library("GeneAnswers")
#library("org.Dm.eg.db")
#library("GO.db")
#fb.entrez <- unlist(as.list(org.Dm.egALIAS2EG))
#genes1<-IMR_up
#cat(genes1)
#iv <- match(genes1, names(fb.entrez))
#iv <- iv[!is.na(iv)]
#iv
#topset <- entrezgene[!is.na(entrezgene)]
#foo <- geneAnswersBuilder(topset, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG')
#go.bp <- foo@enrichmentInfo
#################################################
source("http://bioconductor.org/biocLite.R")
biocLite("GeneAnswers")
citation("GeneAnswers")
library(GeneAnswers)
biocLite("graph")
library(graph)


#data(sampleGroupsData)
#data('humanGeneInput')
#data('humanExpr')
## build a GeneAnswers instance with statistical test based on biological process of GO and saved example data.
#x <- geneAnswersBuilder(humanGeneInput, 'org.Hs.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, geneExpressionProfile=humanExpr)
#x <- geneAnswersBuilder(IMR_up_id, 'org.Dm.eg.db', categoryType='CABIO.PATH', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
#y <- geneAnswersBuilder(DMR_up_id, 'org.Dm.eg.db', categoryType='CABIO.PATH', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
#z <- geneAnswersBuilder(FMR_up_id, 'org.Dm.eg.db', categoryType='CABIO.PATH', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
#a <- geneAnswersBuilder(TMR_up_id, 'org.Dm.eg.db', categoryType='CABIO.PATH', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")

#foo <- geneAnswersBuilder(topset, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG')
#go.bp <- foo@enrichmentInfo
x <- geneAnswersBuilder(IMR_up_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=2, sortBy="pvalue")
y <- geneAnswersBuilder(DMR_up_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
z <- geneAnswersBuilder(FMR_up_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
a <- geneAnswersBuilder(TMR_up_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")

x_f_up <- geneAnswersBuilder(IFR_up_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=2, sortBy="pvalue")
y_f_up <- geneAnswersBuilder(DFR_up_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
z_f_up <- geneAnswersBuilder(FFR_up_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
a_f_up <- geneAnswersBuilder(TFR_up_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")

x_m_down <- geneAnswersBuilder(IMR_down_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=2, sortBy="pvalue")
y_m_down <- geneAnswersBuilder(DMR_down_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
z_m_down <- geneAnswersBuilder(FMR_down_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
a_m_down <- geneAnswersBuilder(TMR_down_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")

x_fm_down <- geneAnswersBuilder(IFR_down_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=2, sortBy="pvalue")
y_fm_down <- geneAnswersBuilder(DFR_down_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
z_fm_down <- geneAnswersBuilder(FFR_down_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")
a_fm_down <- geneAnswersBuilder(TFR_down_id, 'org.Dm.eg.db', categoryType='GO.BP', testType='hyperG',  pvalueT=0.1, FDR.correction=TRUE, level=1, sortBy="pvalue")

go.bp <- x@enrichmentInfo
#'GO', 'KEGG', 'DOLITE', 'REACTOME.PATH', 'CABIO.PATH'

gene_GO <- geneAnswersReadable(x, catTerm=T)
gene_GO_x<-geneAnswersSort(gene_GO, sortBy="pvalue")

gene_GO_y <- geneAnswersReadable(y, catTerm=T)
gene_GO_yy<-geneAnswersSort(gene_GO_y, sortBy="pvalue")

gene_GO_z <- geneAnswersReadable(z, catTerm=T)
gene_GO_zz<-geneAnswersSort(gene_GO_z, sortBy="pvalue")

gene_GO_a <- geneAnswersReadable(a, catTerm=T)
gene_GO_aa<-geneAnswersSort(gene_GO_a, sortBy="pvalue")

gene_GO_f_up <- geneAnswersReadable(x_f_up, catTerm=T)
gene_GO_x_f_up<-geneAnswersSort(gene_GO_f_up, sortBy="pvalue")

gene_GO_y_f_up <- geneAnswersReadable(y_f_up, catTerm=T)
gene_GO_yy_f_up<-geneAnswersSort(gene_GO_y_f_up, sortBy="pvalue")

gene_GO_z_f_up <- geneAnswersReadable(z_f_up, catTerm=T)
gene_GO_zz_f_up<-geneAnswersSort(gene_GO_z_f_up, sortBy="pvalue")

gene_GO_a_f_up <- geneAnswersReadable(a_f_up, catTerm=T)
gene_GO_aa_f_up<-geneAnswersSort(gene_GO_a_f_up, sortBy="pvalue")

gene_GO_m_down <- geneAnswersReadable(x_m_down, catTerm=T)
gene_GO_x_m_down<-geneAnswersSort(gene_GO_m_down, sortBy="pvalue")

gene_GO_y_m_down <- geneAnswersReadable(y_m_down, catTerm=T)
gene_GO_yy_m_down<-geneAnswersSort(gene_GO_y_m_down, sortBy="pvalue")

gene_GO_z_m_down <- geneAnswersReadable(z_m_down, catTerm=T)
gene_GO_zz_m_down<-geneAnswersSort(gene_GO_z_m_down, sortBy="pvalue")

gene_GO_a_m_down <- geneAnswersReadable(a_m_down, catTerm=T)
gene_GO_aa_m_down<-geneAnswersSort(gene_GO_a_m_down, sortBy="pvalue")

gene_GO_fm_down <- geneAnswersReadable(x_fm_down, catTerm=T)
gene_GO_x_fm_down<-geneAnswersSort(gene_GO_fm_down, sortBy="pvalue")

gene_GO_y_fm_down <- geneAnswersReadable(y_fm_down, catTerm=T)
gene_GO_yy_fm_down<-geneAnswersSort(gene_GO_y_fm_down, sortBy="pvalue")

gene_GO_z_fm_down <- geneAnswersReadable(z_fm_down, catTerm=T)
gene_GO_zz_fm_down<-geneAnswersSort(gene_GO_z_fm_down, sortBy="pvalue")

gene_GO_a_fm_down <- geneAnswersReadable(a_fm_down, catTerm=T)
gene_GO_aa_fm_down<-geneAnswersSort(gene_GO_a_fm_down, sortBy="pvalue")

trace(tkplot, edit=TRUE)

topcaBIO.PATH(x, catTerm = TRUE, keepID=TRUE,orderby="pvalue",top="ALL", file=T, fileName="Rad_bioCart.txt")
geneAnswersChartPlots(gene_GO, top = 25, las=2, cex.names=0.5, chartType='pieChart', sortBy="pvalue")
geneAnswersChartPlots(gene_GO, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_x, centroidSize="pvalue")
par(cex.lab=0.1)
geneAnswersConceptNet(gene_GO_x, centroidSize="pvalue", showCats=c(5,17,15,18,19,28,29,30,31), output = "interactive")
geneAnswersConceptNet(gene_GO_x, centroidSize="pvalue", showCats=c(5,17,15,18,19,28,29,30,31))

geneAnswersChartPlots(gene_GO_f_up, top = 35, las=2, cex=0.2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_x_f_up, centroidSize="pvalue", showCats=c(1,9,15), output = "interactive")

geneAnswersChartPlots(gene_GO_m_down, top = 35, las=2, cex=0.2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_x_m_down, centroidSize="pvalue", showCats=c(23,27,28,34), output = "interactive")

geneAnswersChartPlots(gene_GO_fm_down, top = 35, las=2, cex=0.2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_x_fm_down, centroidSize="pvalue", showCats=c(1,5,17), output = "interactive")


topcaBIO.PATH(y, catTerm = TRUE, keepID=TRUE,orderby="pvalue",top="ALL", file=T, fileName="For_bioCart.txt")
geneAnswersChartPlots(gene_GO_y, top = 5, las=2, cex.names=0.5, chartType='pieChart', sortBy="pvalue")
gr<-geneAnswersChartPlots(gene_GO_y, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_yy, centroidSize="pvalue")
geneAnswersConceptNet(gene_GO_yy, centroidSize="pvalue", showCats=c(3,9,10,29,21,23,24,19), output = "interactive")
geneAnswersConceptNet(gene_GO_yy, centroidSize="pvalue", showCats=c(3,9,10,29,21,23,24,19))

gr<-geneAnswersChartPlots(gene_GO_y_f_up, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_yy_f_up, centroidSize="pvalue", showCats=c(4,10,12,20,27), output = "interactive")

gr<-geneAnswersChartPlots(gene_GO_y_m_down, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_yy_m_down, centroidSize="pvalue", showCats=c(1,6,15), output = "interactive")

gr<-geneAnswersChartPlots(gene_GO_y_fm_down, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_yy_fm_down, centroidSize="pvalue", showCats=c(1,13,25,34,35), output = "interactive")


topcaBIO.PATH(z, catTerm = TRUE, keepID=TRUE,orderby="pvalue",top="ALL", file=T, fileName="Dio_bioCart.txt")
geneAnswersChartPlots(gene_GO_z, top = 7, las=2, cex.names=0.5, chartType='pieChart', sortBy="pvalue")
geneAnswersChartPlots(gene_GO_z, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_zz, centroidSize="pvalue")
geneAnswersConceptNet(gene_GO_zz, centroidSize="pvalue", showCats=c(4,6,8,9,10,11,13,32,34), output = "interactive")
geneAnswersConceptNet(gene_GO_zz, centroidSize="pvalue", showCats=c(4,6,8,9,10,11,13,31))

geneAnswersChartPlots(gene_GO_z_f_up, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_zz_f_up, centroidSize="pvalue", showCats=c(4,6,8,9,10,11,13,32,34), output = "interactive")

geneAnswersChartPlots(gene_GO_z_m_down, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_zz_m_down, centroidSize="pvalue", showCats=c(1), output = "interactive")

geneAnswersChartPlots(gene_GO_z_fm_down, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_zz_fm_down, centroidSize="pvalue", showCats=c(1:16), output = "interactive")


topcaBIO.PATH(a, catTerm = TRUE, keepID=TRUE,orderby="pvalue",top="ALL", file=T, fileName="Dio_bioCart.txt")
geneAnswersChartPlots(gene_GO_a, top = 7, las=2, cex.names=0.5, chartType='pieChart', sortBy="pvalue")
geneAnswersChartPlots(gene_GO_a, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_aa, centroidSize="pvalue")
geneAnswersConceptNet(gene_GO_aa, centroidSize="pvalue", showCats=c(3,6,8,9,10,11,15,17), output = "interactive")
geneAnswersConceptNet(gene_GO_aa, centroidSize="pvalue", showCats=c(3,9,15,17,23))

geneAnswersChartPlots(gene_GO_a_f_up, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_aa_f_up, centroidSize="pvalue", showCats=c(3,4,7,13,15), output = "interactive")

trace(tkplot, edit=TRUE)
geneAnswersChartPlots(gene_GO_a_m_down, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_aa_m_down, centroidSize="pvalue", showCats=c(3,4,7,13,15), output = "interactive")

geneAnswersChartPlots(gene_GO_a_fm_down, top = 35, las=2, cex.names=0.6, chartType='barPlot', sortBy="pvalue")
geneAnswersConceptNet(gene_GO_aa_fm_down, centroidSize="pvalue", showCats=c(4,13,23), output = "interactive")

#geneAnswersConceptNet(x, centroidSize="oddsRatio", output='interactive', catTerm=TRUE, geneSymbol=TRUE)
#geneAnswersConceptNet(xx, colorValueColumn='foldChange', centroidSize='pvalue', output='interactive')

geneAnswersHeatmap(x, showCats = c(1:50), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.1)
geneAnswersHeatmap(y, showCats = c(1:25), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.1)
geneAnswersHeatmap(z, showCats = c(1:8), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.1)
geneAnswersHeatmap(a, showCats = c(1:16), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.1)

geneAnswersHeatmap(x, showCats = c(5,17,15,18,19,28,29,30,31), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.8)
#geneAnswersHeatmap(y, showCats = c(3,10,29,21,23,24,19), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.1)
geneAnswersHeatmap(y, showCats = c(3,10,24,21,23,19), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.8)
#geneAnswersHeatmap(z, showCats = c(4,6,8,9,10,11,13,31), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.1)
geneAnswersHeatmap(z, showCats = c(4,6,8,10,13,31), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.8)
geneAnswersHeatmap(a, showCats = c(3,9,15,17,23), catTerm=T, catID=F, geneSymbol=TRUE, showAllGenes=F, cex=1.1)
                      par(mar=c(0,0,0,0))
                      geneAnnotationHeatmap()                      
gAList <- lapply(Sample, geneAnswersBuilder, 'org.Dm.eg.db', categoryType='GO.BP', pvalueT=0.05, verbose=FALSE)
#gAList <- lapply(Sample, geneAnswersBuilder, 'org.Dm.eg.db', categoryType='KEGG', pvalueT=0.01, verbose=FALSE)

#output<- getConceptTable(gAList, items='geneNum')
#output<- getConceptTable(gAList, items='both')
output<- getConceptTable(gAList, catTerm = FALSE, items='geneNum')
output<- getConceptTable(gAList, catTerm = FALSE, items='both')

output
groupReport(output[[1]], gAList,  matrixOfHeatmap=output[[2]], clusterTable=NULL, fileName='GO_BP.html',  catType='GO', colorValueColumn=colnames(Sample[[1]])[-1]) 
groupReport(output[[1]], gAList,  matrixOfHeatmap=output[[2]], clusterTable=NULL, fileName='KEGGtest.html',  catType='KEGG', colorValueColumn=colnames(Sample[[1]])[-1]) 
groupReport(output[[1]], gAList,  matrixOfHeatmap=output[[2]], clusterTable=NULL, fileName='CABIO.html',  catType='CABIO.PATH', colorValueColumn=colnames(Sample[[1]])[-1]) 

drawTable(output[[1]], topCat=35, mar=c(5,25,4,2), heatMap=TRUE, matrixOfHeatmap=NULL, clusterTable=NULL, addRowLabel=TRUE, cex.axis=c(0.9, 0.7), reverseOfCluster=FALSE, xGridLine=FALSE, colorBar=TRUE, newWindow=T)

drawTable(output[[1]], topCat=35, mar=c(5,25,4,2), heatMap=TRUE, matrixOfHeatmap=output, clusterTable=NULL, addRowLabel=TRUE, cex.axis=c(0.9, 0.7), reverseOfCluster=FALSE, xGridLine=FALSE, colorBar=TRUE, newWindow=T)

#############################################
library("limma")
library("GOsummaries")
query_up_m<-list(IMR_up, DMR_up, FMR_up, TMR_up)
gs_up_m = gosummaries(query_up_m, organism = "dmelanogaster",go_branches = c("BP","ke","re"), max_p_value="0.05", min_set_size=1,max_set_size=1000,max_signif=1000)
gs_up_m = add_to_slot.gosummaries(gs_up_m, "Title", list("IMR_up", "DMR_up", "FMR_up", "TMR_up"))
plot(gs_up_m, fontsize = 12, cex=1.2)

query_up_m<-list(IMR_up, DMR_up)
gs_up_m = gosummaries(query_up_m, organism = "dmelanogaster",go_branches = c("BP","ke","re"), max_p_value="0.05", min_set_size=1,max_set_size=1000,max_signif=1000)
gs_up_m = add_to_slot.gosummaries(gs_up_m, "Title", list("IMR_up", "DMR_up"))
plot(gs_up_m, fontsize = 12, cex=1.2)


query_up_fm<-list(IFR_up, DFR_up, FFR_up, TFR_up)
gs_up_fm = gosummaries(query_up_fm, organism = "dmelanogaster",go_branches = c("BP","ke","re"), max_p_value="0.05", min_set_size=1,max_set_size=1000,max_signif=1000)
gs_up_fm = add_to_slot.gosummaries(gs_up_fm, "Title", list("IFR_up", "DFR_up", "FFR_up", "TFR_up"))
plot(gs_up_fm, fontsize = 16, cex=1.6)



query_up<-list(IMR_up, IFR_up, DMR_up, DFR_up, FMR_up, FFR_up, TMR_up, TFR_up)
gs_up = gosummaries(query_up, organism = "dmelanogaster")
gs_up = add_to_slot.gosummaries(gs_up, "Title", list("IMR_up", "IFR_up", "DMR_up", "DFR_up", "FMR_up", "FFR_up", "TMR_up", "TFR_up"))
plot(gs_up, fontsize = 12, cex=1.2)


query_down<-list(IMR_down, IFR_down, DMR_down, DFR_down, FMR_down, FFR_down, TMR_down, TFR_down)
gs_down = gosummaries(query_down, organism = "dmelanogaster")
gs_down = add_to_slot.gosummaries(gs_down, "Title", list("IMR_down", "IFR_down", "DMR_down", "DFR_down", "FMR_down", "FFR_down", "TMR_down", "TFR_down"))
plot(gs_down, fontsize = 12, cex=1.2)

query_down_m<-list(IMR_down, DMR_down, FMR_down, TMR_down)
gs_down_m = gosummaries(query_down_m, organism = "dmelanogaster")
gs_down_m = add_to_slot.gosummaries(gs_down_m, "Title", list("IMR_down", "DMR_down", "FMR_down", "TMR_down"))
plot(gs_down_m, fontsize = 12, cex=1.2)

query_down_fm<-list(IFR_down, DFR_down, FFR_down, TFR_down)
gs_down_fm = gosummaries(query_down_fm, organism = "dmelanogaster")
gs_down_fm = add_to_slot.gosummaries(gs_down_fm, "Title", list("IFR_down", "DFR_down", "FFR_down", "TFR_down"))
plot(gs_down_fm, fontsize = 12, cex=1.2)