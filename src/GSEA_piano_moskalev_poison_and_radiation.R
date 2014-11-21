# local tags
# TODO -- items for the followup processing


#######
# loadMAdata() --> extractFactors() --> diffExp()

###################################################
###################################################
### Load packages
###################################################
source("http://www.bioconductor.org/biocLite.R")
biocLite("piano", dependencies=TRUE)
biocLite("org.Dm.eg.db")
biocLite("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
biocLite("BSgenome.Dmelanogaster.UCSC.dm3")



library(piano)
if(!try(require(affy))) stop("package affy is missing")
if(!try(require(plier))) stop("package plier is missing")
library("org.Dm.eg.db")
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
library("BSgenome.Dmelanogaster.UCSC.dm3")





###################################################
###################################################
### load and prepare data for processing
###################################################
# The data discussed in this publication have been deposited in NCBIâs Gene Expression Omnibus and are accessible through
# GEO Series accession number GSE50377 http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50377
# all.txt that contains reads is provided by Alexey Moskalev


# setup directory
setwd("C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertibrates/GSEA/GSEA_vc01")
# read data
df_all <- read.csv("./Data/all.txt", header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE)

# df_MAsetup
# rownames(df_MAsetup) <- c("CMR1", "CMR2", "CFR1", "CFR2",
#                           "IMR1", "IMR2", "IFR1", "IFR2",
#                           "DMR1", "DMR2", "DFR1", "DFR2",
#                           "FMR1", "FMR2", "FFR1", "FFR2",
#                           "TMR1", "TMR2", "TFR1", "TFR2")
# 
# colnames(df_MAsetup) <- c("pollutant","sex")

# load data
    setwd("C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertibrates/GSEA/GSEA_vc01")
    dataRaw <- read.csv("./Data/all.txt", header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE)
    setup <- read.csv("./Data/setup.txt", header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE)
    setup <- lapply(setup, as.character)
# Verbose function:, as copied from loadMAdata.r -- a source code of R:Piano
    verbose = TRUE
    .verb <- function(mes, verbose) {
        if(verbose == TRUE) {
            message(mes)
        }
    }

# Normalize the raw data:, as copied and slightly modified from loadMAdata.r -- a source code of R:Piano
    normalization <- "plier"


    if(exists("dataRaw", inherits=FALSE)) {
        # iterplier qubic spline
        if(normalization == "plier") {
            .verb("Preprocessing using PLIER with cubic spline normalization...", verbose)
            dataNorm <- normalize.qspline(data.matrix(dataRaw), verbose=FALSE)
            dataNorm <- as.data.frame(dataNorm)
            colnames(dataNorm) <- colnames(dataRaw)
            rownames(dataNorm) <- rownames(dataRaw)
            .verb("...done", verbose)
        } else {
        .verb("Text file data: No normalization performed.", verbose)
        }
    }

# annotation, as copied and slightly modified from loadMAdata.r -- a source code of R:Piano
# no annotation is used in the current implementation
    
    # flybase transcripts
    txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
    bsdb <- BSgenome.Dmelanogaster.UCSC.dm3

    # according to http://www.bioconductor.org/help/course-materials/2012/Bressanone2012/2012-07-05-Morgan-annotation.pdf
    biocLite("AnnotationDbi", dependencies=TRUE)
    library("AnnotationDbi")
    # gene-centric discovery & selection
        # context: DESeq top table
        library("DESeq")
        countTable <- read.csv("./Data/all.txt", header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE)
        condition <- factor(c("Control_m","Control_m","Control_fm","Control_fm","Irradiation_m","Irradiation_m","Irradiation_fm","Irradiation_fm","Dioxin_m","Dioxin_m","Dioxin_fm","Dioxin_fm","formaldehyde_m","formaldehyde_m","formaldehyde_fm","formaldehyde_fm","toluene_m","toluene_m","toluene_fm","toluene_fm"))
        cds <- newCountDataSet(countTable, condition)
        cds <- estimateSizeFactors( cds )
        res <- nbinomTest(cds, "Control_fm","toluene_fm")
        topTable = res[order(res$pval),]

        # discover & select
        library(org.Dm.eg.db)
        keytypes(org.Dm.eg.db)
        columns(org.Dm.eg.db)
        head(keys(org.Dm.eg.db, keytype="FLYBASE"))
        fbids <- topTable$id[1:5]
        cols <- c("ENTREZID","GENENAME", "CHRLOC")
        anno <- select(x = org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENSEMBLTRANS")
        head(anno)
        # merging top table and annotation
#         fbids <- topTable$id
        fbids <- res$id
        cols <- c("GENENAME", "CHRLOC")
        anno <- select(x = org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENSEMBLTRANS")
        resAnno <- merge(res, anno, by.x="id", by.y = "ENSEMBLTRANS", all.x = TRUE)
        
        # prepare annotation
        # though this line would not work because of non-unique transcript id
        rownames(resAnno) <- resAnno$id

    # flybase
    flyID <- org.Dm.egFLYBASE
    flyID <- toTable(flyID)
    flyIF <- flyID[duplicated(flyID$flybase_id)==FALSE,]
    
    # Gene name
    geneName <- org.Dm.egGENENAME
    geneName <- toTable(geneName)
    # Chromosome location
    chromosome <- org.Dm.egCHRLOC
    chromosome <- toTable(chromosome)
    chromosome <- chromosome[,c(1,3,2)]
    # Annotation data frame
    annot <- merge(flyID,geneName,by.x="gene_id",by.y="gene_id")
    annot <- merge(annot,chromosome,by.x="gene_id",by.y="gene_id")
    rownames(annot) <- annot$flybase_id
    annot <- annot[2:ncol(annot)]
    colnames(annot) <- c("geneName","chromosome","start") # <- remove sys.name?
    .verb("...done", verbose)


# load data with R:Piano:loadMAdata()
# TODO -- deliver parameter independent script. Though getwd() & setup.txt are generic already.


ma_data <- loadMAdata(datadir = getwd(), 
                      setup = "/Data/setup.txt", 
                      dataNorm = dataNorm,
                      platform = "NULL", 
                      normalization = "plier",
                      filter = TRUE, 
                      verbose = TRUE)


ma_datam <- list(dataNorm=dataNorm, setup=setup)
class(ma_datam) = "ArrayData"

extractFactors(ma_datam)


###################################################
### quality control
runQC(ma_datam)


options(device = "windows")
# To only run the PCA:
runQC(ma_data, rnaDeg=FALSE, nuseRle=FALSE, hist=FALSE, boxplot=FALSE, pca=TRUE)
# Additionally, for the PCA you can specify other colors:
runQC(ma_data, rnaDeg=FALSE, nuseRle=FALSE, hist=FALSE, boxplot=FALSE, pca=TRUE, colors=c("cyan","orange"))

runQC(ma_data, rnaDeg=TRUE, nuseRle=FALSE, hist=TRUE, boxplot=TRUE, pca=TRUE)


###################################################
### Differential expression analysis
extractFactors(ma_data)


pfc <- diffExp(ma_data, contrasts=c("control_female", "dioxin_male - dioxin_female","control_male - dioxin_male"))

# Sort genes in "aerobic_Clim - anaerobic_Clim" according to adjusted p-value:
ii <- sort(pfc$resTable[[1]]$adj.P.Val, index.return=TRUE)$ix
pfc$resTable[[1]][ii[1:5],]
# Sort genes in "aerobic_Nlim - anaerobic_Nlim" according to adjusted p-value:
ii <- sort(pfc$resTable[[2]]$adj.P.Val, index.return=TRUE)$ix
pfc$resTable[[2]][ii[1:5],]


###################################################
###################################################
### Gene set analysis
###################################################
### GSA input data

# Get p-values from the aerobic_Clim vs anaerobic_Clim comparison:
myPval <- pfc$pValues["aerobic_Clim - anaerobic_Clim"]
# Display the first values and gene IDs:
head(myPval)

# Custom gene to gene set mapping:
genes2genesets <- cbind(paste("gene",c("A","A","A","B","B","C","C","C","D"),sep=""),
                        paste("set",c(1,2,3,1,3,2,3,4,4),sep=""))
genes2genesets
# Load into correct format:
myGsc <- loadGSC(genes2genesets)
# View summary:
myGsc
# View all gene sets:
myGsc$gsc

myStats <- c(-1.5,-0.5,1,2)
names(myStats) <- paste("gene",c("A","B","C","D"),sep="")
myStats


gsaRes <- runGSA(myStats, gsc=myGsc)

gsaRes

names(gsaRes)




#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
options(device = "windows")
library(DESeq)
file<-"C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertibrates/GSEA/GSEA_vc01/Data/all.txt"
countsTable <- read.delim(file=file, header=TRUE, stringsAsFactors=TRUE )
colnames(countsTable)
rownames(countsTable) <- countsTable$Transcript
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


#################################################
#самки
#################################################


res_fm <- nbinomTest(cds, "Control_fm","toluene_fm")
#diff_names<-res$id

diff_fm<-head(res_fm,n=100)
diff_fm

plotMA(res_fm)

#     baseMean          The base mean (i.e., mean of the counts divided by the size factors) for the counts for both conditions
#     log2FoldChange    The log2 of the fold change
#     padj              The adjusted p values (adjusted with ’p.adjust( pval, method="BH")’)
plot(res_fm$baseMean,res_fm$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res_fm$padj < .1, "red", "black" ) )

#     pval              The p value for rejecting the null hypothesis ’meanA==meanB’
hist(res_fm$pval, breaks=100, col="skyblue", border="slateblue", main="")

#     padj              The adjusted p values (adjusted with ’p.adjust( pval, method="BH")’)
res_fmSig <- res_fm[ res_fm$padj < .1, ]

#     pval              The p value for rejecting the null hypothesis ’meanA==meanB’
diff_sig_fm<-head( res_fmSig[ order(res_fmSig$pval), ], n=200 )
diff_sig_fm

#     pval              The p value for rejecting the null hypothesis ’meanA==meanB’
diff_sig_names_fm<-head( res_fmSig[ order(res_fmSig$pval), ], n=200)

#     id                The ID of the observable, taken from the row names of the counts slots.
diff_sig_names_fm<-diff_sig_names_fm$id
diff_sig_names_fm
#DOWN
#     foldChange        The ratio meanB/meanA
#     baseMean          The base mean (i.e., mean of the counts divided by the size factors) for the counts for both conditions
down_fm<-head( res_fmSig[ order( res_fmSig$foldChange, -res_fmSig$baseMean ), ],n=200 )
down_fm
#UP
#     foldChange        The ratio meanB/meanA
#     baseMean          The base mean (i.e., mean of the counts divided by the size factors) for the counts for both conditions
up_fm<-head( res_fmSig[ order( -res_fmSig$foldChange, -res_fmSig$baseMean ), ],n=200 )
up_fm

write.csv(res_fm, file="toluene_fm.csv")
read.csv("toluene_fm.csv", row.names = 1)


#################################################
#самцы
#################################################
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

#################################################
#Множественные сравнения
#################################################

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

