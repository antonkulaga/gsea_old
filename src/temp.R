library("AnnotationDbi")
library("GenomicFeatures")
library("DESeq")
library("org.Dm.eg.db")
library("GO.db")


setwd("C:/Users/user/Google Drive/Work/Semantic research framework/RFBR 2014, Aging gene ontology on invertibrates/GSEA/GSEA_vc01")

countTable <- read.csv("./Data/all.txt", header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE)
condition <- factor(c("Control_m","Control_m","Control_fm","Control_fm","Irradiation_m","Irradiation_m","Irradiation_fm","Irradiation_fm","Dioxin_m","Dioxin_m","Dioxin_fm","Dioxin_fm","formaldehyde_m","formaldehyde_m","formaldehyde_fm","formaldehyde_fm","toluene_m","toluene_m","toluene_fm","toluene_fm"))
cds <- newCountDataSet(countTable, condition)
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
res <- nbinomTest(cds, "Control_fm","toluene_fm")
topTable = res[order(res$pval),]
# 
# discover & select
keytypes(org.Dm.eg.db)
columns(org.Dm.eg.db)
head(keys(org.Dm.eg.db, keytype="ENSEMBLTRANS"))
fbids <- topTable$id[1:5]
cols <- c("ENTREZID", "GO")
anno <- select(org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENSEMBLTRANS")
head(anno)
merging top table and annotation
fbids <- topTable$id
# fbids <- res$id
cols <- c("ENTREZID")
anno <- select(x = org.Dm.eg.db, keys = fbids, columns = cols, keytype = "ENSEMBLTRANS")
resAnno <- merge(res, anno, by.x="id", by.y = "ENSEMBLTRANS", all.x = TRUE)

# prepare annotation
# though this line would not work because of non-unique transcript id
rownames(resAnno) <- resAnno$id




    if (fitMethod == "nbinom") {
        dataForLimma <- arrayData$dataCount
        res <- nbinomTest(cds, "Control_fm","toluene_fm")
        
        
        countTable <- read.csv("./Data/all.txt", header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE)
        condition <- factor(c("Control_m","Control_m","Control_fm","Control_fm","Irradiation_m","Irradiation_m","Irradiation_fm","Irradiation_fm","Dioxin_m","Dioxin_m","Dioxin_fm","Dioxin_fm","formaldehyde_m","formaldehyde_m","formaldehyde_fm","formaldehyde_fm","toluene_m","toluene_m","toluene_fm","toluene_fm"))
        cds <- newCountDataSet(countTable[1:100,], condition)
        cds <- estimateSizeFactors( cds )
        cds <- estimateDispersions( cds )
        res <- nbinomTest(cds, "Control_fm","toluene_fm")
        topTable = res[order(res$pval),]
        
        
        
    }
	
> head(res)
           id   baseMean  baseMeanA  baseMeanB foldChange log2FoldChange      pval padj
1 FBtr0005009 11447.7315 12171.4427 10724.0202  0.8810804    -0.18265435 0.9000647    1
2 FBtr0005088  8483.8044  8355.3823  8612.2264  1.0307400     0.04368041 0.9774713    1
3 FBtr0005673  1117.5051  1135.5110  1099.4991  0.9682858    -0.04649520 0.5560763    1
4 FBtr0005674  3189.5546  3232.8322  3146.2770  0.9732262    -0.03915294 0.7730619    1
5 FBtr0006151 27784.9850 26455.4474 29114.5225  1.1005114     0.13817414 0.7552802    1
6 FBtr0070000   274.6608   304.3641   244.9576  0.8048174    -0.31326655 0.9529320    1
> str(dexp)
List of 4
 $ pValues    :'data.frame':	25415 obs. of  4 variables:
  ..$ irradiation_female - control_female : num [1:25415] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ dioxin_female - control_female      : num [1:25415] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ formaldehyde_female - control_female: num [1:25415] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ toluene_female - control_female     : num [1:25415] 1 1 1 1 1 1 1 1 1 1 ...
 $ foldChanges:'data.frame':	25415 obs. of  4 variables:
  ..$ irradiation_female - control_female : num [1:25415] -0.158 -0.134 -0.177 -0.158 0.172 ...
  ..$ dioxin_female - control_female      : num [1:25415] 0.1377 -0.0129 -0.2716 -0.1782 0.7712 ...
  ..$ formaldehyde_female - control_female: num [1:25415] -0.0708 -0.0975 -0.3076 -0.6961 0.0497 ...
  ..$ toluene_female - control_female     : num [1:25415] -0.0724 0.2018 0.0381 0.0583 0.306 ...
 $ resTable   :List of 4
  ..$ irradiation_female - control_female :'data.frame':	25415 obs. of  7 variables:
  .. ..$ ProbesetID: Factor w/ 25415 levels "FBtr0005009",..: 1 2 3 4 5 6 7 8 9 10 ...
  .. ..$ logFC     : num [1:25415] -0.158 -0.134 -0.177 -0.158 0.172 ...
  .. ..$ AveExpr   : num [1:25415] 14.1 12.8 9.4 10.4 16 ...
  .. ..$ t         : num [1:25415] -0.88 -1.744 -0.59 -0.805 0.591 ...
  .. ..$ P.Value   : num [1:25415] 0.398 0.109 0.567 0.438 0.566 ...
  .. ..$ adj.P.Val : num [1:25415] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ B         : num [1:25415] -5.77 -4.77 -5.98 -5.83 -5.98 ...
  ..$ dioxin_female - control_female      :'data.frame':	25415 obs. of  7 variables:
  .. ..$ ProbesetID: Factor w/ 25415 levels "FBtr0005009",..: 1 2 3 4 5 6 7 8 9 10 ...
  .. ..$ logFC     : num [1:25415] 0.1377 -0.0129 -0.2716 -0.1782 0.7712 ...
  .. ..$ AveExpr   : num [1:25415] 14.1 12.8 9.4 10.4 16 ...
  .. ..$ t         : num [1:25415] 0.768 -0.167 -0.904 -0.911 2.655 ...
  .. ..$ P.Value   : num [1:25415] 0.459 0.8701 0.3855 0.3822 0.0226 ...
  .. ..$ adj.P.Val : num [1:25415] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ B         : num [1:25415] -6.04 -6.33 -5.93 -5.92 -3.49 ...
  ..$ formaldehyde_female - control_female:'data.frame':	25415 obs. of  7 variables:
  .. ..$ ProbesetID: Factor w/ 25415 levels "FBtr0005009",..: 1 2 3 4 5 6 7 8 9 10 ...
  .. ..$ logFC     : num [1:25415] -0.0708 -0.0975 -0.3076 -0.6961 0.0497 ...
  .. ..$ AveExpr   : num [1:25415] 14.1 12.8 9.4 10.4 16 ...
  .. ..$ t         : num [1:25415] -0.395 -1.269 -0.984 -2.619 0.171 ...
  .. ..$ P.Value   : num [1:25415] 0.7006 0.2311 0.3466 0.0241 0.8672 ...
  .. ..$ adj.P.Val : num [1:25415] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ B         : num [1:25415] -6.48 -5.77 -6.03 -3.5 -6.55 ...
  ..$ toluene_female - control_female     :'data.frame':	25415 obs. of  7 variables:
  .. ..$ ProbesetID: Factor w/ 25415 levels "FBtr0005009",..: 1 2 3 4 5 6 7 8 9 10 ...
  .. ..$ logFC     : num [1:25415] -0.0724 0.2018 0.0381 0.0583 0.306 ...
  .. ..$ AveExpr   : num [1:25415] 14.1 12.8 9.4 10.4 16 ...
  .. ..$ t         : num [1:25415] -0.363 2.397 0.118 0.252 1.054 ...
  .. ..$ P.Value   : num [1:25415] 0.7236 0.0357 0.9084 0.8055 0.315 ...
  .. ..$ adj.P.Val : num [1:25415] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..$ B         : num [1:25415] -6.15 -3.83 -6.24 -6.13 -5.76 ...
 $ vennMembers:List of 15
  ..$ Uniquely in A   : Factor w/ 25415 levels "FBtr0005009",..: 11228 15419 15580 16151 22972
  ..$ Uniquely in B   : Factor w/ 25415 levels "FBtr0005009",..: 1870 4268 4879 5369 6347 6591 12667 16741 16778 17261 ...
  ..$ Uniquely in AB  : Factor w/ 25415 levels "FBtr0005009",..: 
  ..$ Uniquely in C   : Factor w/ 25415 levels "FBtr0005009",..: 5219 5220 8116 10275 14429 15914
  ..$ Uniquely in AC  : Factor w/ 25415 levels "FBtr0005009",..: 
  ..$ Uniquely in BC  : Factor w/ 25415 levels "FBtr0005009",..: 
  ..$ Uniquely in ABC : Factor w/ 25415 levels "FBtr0005009",..: 4878 14412
  ..$ Uniquely in D   : Factor w/ 25415 levels "FBtr0005009",..: 3484 4550 4551 18682 19678 19837 20251 21376
  ..$ Uniquely in AD  : Factor w/ 25415 levels "FBtr0005009",..: 
  ..$ Uniquely in BD  : Factor w/ 25415 levels "FBtr0005009",..: 11239 17154
  ..$ Uniquely in ABD : Factor w/ 25415 levels "FBtr0005009",..: 13274 13275
  ..$ Uniquely in CD  : Factor w/ 25415 levels "FBtr0005009",..: 
  ..$ Uniquely in ACD : Factor w/ 25415 levels "FBtr0005009",..: 
  ..$ Uniquely in BCD : Factor w/ 25415 levels "FBtr0005009",..: 
  ..$ Uniquely in ABCD: Factor w/ 25415 levels "FBtr0005009",..: 
  
  
  
  
  
  
  
  
  
  
fct <- rownames(subset(extractFactors(l_da)$factors, factors == "Control_fm" | factors == "Irradiation_fm"))

