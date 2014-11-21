#setwd("C:/Users/Sony/Documents")
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq")

source("http://bioconductor.org/biocLite.R")
biocLite("goseq")

source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

source("http://bioconductor.org/biocLite.R")
biocLite("gene2pathway")

source("http://bioconductor.org/biocLite.R")
biocLite("RCytoscape")

source("http://bioconductor.org/biocLite.R")
biocLite("limma")


install.packages("gplots")
install.packages("VennDiagram")

#################################################################################################

library(DESeq)
library(goseq)
library(biomaRt)
library(gene2pathway)
#library(RCytoscape)
#library(limma)
citation("DESeq")

file<-"C:/Users/Sony/Google Диск/ФЦП/RNASeq/reads/all.txt"
countsTable <- read.delim(file=file, header=TRUE, stringsAsFactors=TRUE )
colnames(countsTable)
rownames(countsTable) <- countsTable$Gene
gene_names<-rownames(countsTable)
gene_names
countsTable <- countsTable[ , -1 ]
conds <- factor(c("Control_m","Control_m","Control_fm","Control_fm","Irradiation_m","Irradiation_m","Irradiation_fm","Irradiation_fm","Dioxin_m","Dioxin_m","Dioxin_fm","Dioxin_fm","formaldehyde_m","formaldehyde_m","formaldehyde_fm","formaldehyde_fm","toluene_m","toluene_m","toluene_fm","toluene_fm"))
cds <- newCountDataSet(countsTable, conds)
#head( counts(cds) )
cds <- estimateSizeFactors( cds )
sizeFactors( cds )

head( counts( cds, normalized=TRUE ) )

#cds <- estimateDispersions( cds, method="blind", sharingMode="fit-only" )
cds <- estimateDispersions( cds )

#самки

res_fm <- nbinomTest(cds, "Control_fm","toluene_fm")
#diff_names<-res$id

diff_fm<-head(res_fm,n=100)
diff_fm

plotMA(res_fm)
plot(res_fm$baseMean,res_fm$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res_fm$padj < .1, "red", "black" ) )


hist(res_fm$pval, breaks=100, col="skyblue", border="slateblue", main="")

res_fmSig <- res_fm[ res_fm$padj < .1, ]

diff_sig_fm<-head( res_fmSig[ order(res_fmSig$pval), ], n=200 )
diff_sig_fm

diff_sig_names_fm<-head( res_fmSig[ order(res_fmSig$pval), ], n=200)
diff_sig_names_fm<-diff_sig_names_fm$id
diff_sig_names_fm
#DOWN
down_fm<-head( res_fmSig[ order( res_fmSig$foldChange, -res_fmSig$baseMean ), ],n=200 )
down_fm
#UP
up_fm<-head( res_fmSig[ order( -res_fmSig$foldChange, -res_fmSig$baseMean ), ],n=200 )
up_fm

write.csv(res_fm, file="toluene_fm.csv")
read.csv("toluene_fm.csv", row.names = 1)

#самцы
res_m <- nbinomTest(cds, "Control_m","toluene_m")

diff_m<-head(res_m,n=100)
diff_m

plotMA(res_m)
plot(res_m$baseMean,res_m$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res_m$padj < .1, "red", "black" ) )


hist(res_m$pval, breaks=100, col="skyblue", border="slateblue", main="")

res_mSig <- res_m[ res_m$padj < .1, ]

diff_sig_m<-head( res_mSig[ order(res_mSig$pval), ],n=200 )
diff_sig_m

diff_sig_names_m<-head( res_mSig[ order(res_mSig$pval), ],n=200)
diff_sig_names_m<-diff_sig_names_m$id
diff_sig_names_m
#DOWN
down_m<-head( res_mSig[ order( res_mSig$foldChange, -res_mSig$baseMean ), ],n=200 )
down_m
#UP
up_m<-head( res_mSig[ order( -res_mSig$foldChange, -res_mSig$baseMean ), ],n=200 )
up_m

write.csv(res_m, file="toluene_m.csv")
read.csv("toluene_m.csv", row.names = 1)

#Множественные сравнения

head(countsTable)
irradiationDesign = data.frame(
  row.names = colnames(countsTable),
  condition = conds,
  libType = c( "males", "males", "females",
               "females", "males", "males", "females", "females",
               "males", "males", "females",
               "females", "males", "males", "females", "females", 
               "males", "males", "females", "females" ) )
irradiationDesign
cdsFull = newCountDataSet( countsTable, irradiationDesign )
cdsFull = estimateSizeFactors( cdsFull )
cdsFull = estimateDispersions( cdsFull )

plotDispEsts( cdsFull )

cdsFullBlind = estimateDispersions( cdsFull, method = "blind" )
 vsdFull = varianceStabilizingTransformation( cdsFullBlind )
 library("RColorBrewer")
 library("gplots")
 select = order(rowMeans(counts(cdsFull)), decreasing=TRUE)[1:30]
 hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
 heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
#heatmap.2(exprs(vsdFull)[select,], labRow=countsTable[select,c('Gene')], col = hmcol, trace="none", margin=c(10, 6))
 dists = dist( t( exprs(vsdFull) ) )
 mat = as.matrix( dists )
 rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(condition, libType, sep=" : "))
 heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
 print(plotPCA(vsdFull, intgroup=c("condition", "libType"), ntop = 500))
plotPCA( vsdFull )
pca<-prcomp(t(exprs(vsdFull)[select, ]))
summary(pca)
plot(pca)

#Венская диаграмма совпадений
all<-read.csv("all_dif.csv", header = T, sep=";")
attach(all)
DFR
names(all)

intersect(IFR, DFR)
intersect(IFR, FFR)
intersect(IFR, TFR)
intersect(DFR, FFR)
intersect(DFR, TFR)
intersect(FFR, TFR)

intersect(IMR, DMR)
intersect(IMR, FMR)
intersect(IMR, TMR)
intersect(FMR, DMR)
intersect(DMR, TMR)
intersect(FMR, TMR)
is.element(DMR, FMR)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    "IFR" = IFR,
    "DFR" = DFR
    #"FFR" = FFR,
    #"FMR" = IMR,
    #"TFR" = IFR
  ),
  filename = "all_dif.tiff",
  na="remove",
  scaled = TRUE,
  ext.text = TRUE,
  ext.line.lwd = 2,
  ext.dist = -0.15,
  ext.length = 0.9,
  ext.pos = -4,
  inverted = TRUE,
  cex = 2.5,
  cat.cex = 2.5,
  rotation.degree = 45,
  main = "Complex Venn Diagram",
  #sub = "Featuring: rotation and external lines",
  main.cex = 2,
  sub.cex = 1
);



________________________________________________________________
d = newCountDataSet( countsTable, irradiationDesign )
d <- newCountDataSet(counttable, meta)

## Estimate library size and dispersion
d <- estimateSizeFactors(d)
d <- estimateDispersions(d)
plotDispEsts(d, main="DESeq: Per-gene dispersion estimates")

## Principal components biplot on variance stabilized data, color-coded by condition-librarytype
print(plotPCA(varianceStabilizingTransformation(d), intgroup=c("condition", "libType")))

## Fit full and reduced models, get p-values
dfit1 <- fitNbinomGLMs(d, count~libType+condition)
dfit0 <- fitNbinomGLMs(d, count~libType)
dpval <- nbinomGLMTest(dfit1, dfit0)
dpadj <- p.adjust(dpval, method="BH")

## Make results table with pvalues and adjusted p-values
dtable <- transform(dfit1, pval=dpval, padj=dpadj)
dtable <- dtable[order(dtable$padj), ]
head(dtable)

