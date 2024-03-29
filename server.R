library(shiny)

packages <- c("DESeq","affy","affyPLM","plier","limma","biomaRt","org.Dm.eg.db","igraph","marray","AnnotationDbi","gplots")
baseFolder <- "./src/"
source("http://bioconductor.org/biocLite.R")     
sources <- c("loadData.R","extractFactors.R","runQC.R","diffExp.R","loadGSC.R","checkLoadArg.R","GSCstatBatch.R","runGSA.R","runGSAhyper.R","networkPlot.R")

loadPackages <- function(pkgs,base,src) {
  library()  
    for(lib in pkgs){     
      if( !require(character.only = TRUE,package = lib, quietly = TRUE) ) {
        print(paste("trying to install",lib,"..."))
        biocLite(lib, dependencies=TRUE)
        #library(lib)
        #if(!try(biocLite(lib))) 
        #stop(paste("package ",lib,"is missing"))                          
      }
    }              
    for(s in src) {
      path <- paste(base,s,sep="")
      source(path)    
    }
}

#setwd("/var/shiny-server/www/GSEA_vc01/")

loadPackages(packages,baseFolder,sources)

lda <- loadData(datadir = getwd(), countfile = "./Data/all.txt", setup="./Data/setup.txt", verbose=TRUE)

shinyServer(
    function(input, output){
        output$osetup <- renderPrint({
            if(input$contr1 != "None" && input$contr2 != "None" && input$contr1 != input$contr2) {
                paste(input$contr1, input$contr2, sep=" compared to ")
            } else {
                "Please select two different contrasts."
            }            
        })
        
        
        output$histPlot <- renderPlot({
            if(input$contr1 != "None" && input$contr2 != "None" && input$contr1 != input$contr2) {
                
                fct <- rownames(subset(extractFactors(lda)$factors, factors == input$contr1 | factors == input$contr2))
                sublda <- lda$dataNorm[,fct]
                
                plot(NULL,
                     main="Histogram of normalized data",
                     xlim = c(0,20), ylim=c(0,0.21),
                     xlab=expression(log[2] ~ intensity),ylab="Density")
                for(i in 1:dim(sublda)[2]){
                    lines(density(sublda[,i]),
                          lwd=2, col=i)
                }
                
                legend('topright',
                       ncol=2,
                       legend=colnames(sublda),
                       fill=1:dim(sublda)[2]
                    )
            } else {
                plot(NULL, main="Histogram of normalized data", xlim = c(0,20), ylim=c(0,0.21))
                legend('topleft',
                       ncol=2,
                       legend='Please select two different contrasts.')
            }
        })
        
        
        output$boxplotPlot <- renderPlot({
            if(input$contr1 != "None" && input$contr2 != "None" && input$contr1 != input$contr2) {
                fct <- rownames(subset(extractFactors(lda)$factors, factors == input$contr1 | factors == input$contr2))
                sublda <- lda$dataNorm[,fct]
                
                boxplot(sublda, main="Boxplot of normalized data",
                        ylab=expression(log[2] ~ intensity), las=2)
            } else {
                plot(NULL,main="Boxplot of normalized data", xlim = c(0,20), ylim=c(0,0.21))
                legend('topleft',
                       ncol=2,
                       legend='Please select two different contrasts.')
            }
        })
        
        
        output$networkPlot <- renderPlot({
            if(input$contr1 != "None" && input$contr2 != "None" && input$contr1 != input$contr2) {
                
                
                dxp <- diffExp(lda, 
                               contrasts=paste(input$contr1, input$contr2, sep=" - "),
                               fitMethod=input$dexpmeth, 
                               adjustMethod="BY", significance=input$dexpsig,
                               plot=FALSE, heatmapCutoff=1e-10, volcanoFC=2
                    )
                
                mGsc <- loadGSC(dxp, group=input$gscgroup)
                
                contrast <- paste(input$contr1, input$contr2, sep=" - ")
                if (file.exists(paste("./Data/GSA_hyper_",contrast,"_",input$gscgroup,"_",input$GSAcutoff,".rds",sep=""))) {
                    gsah<- readRDS(paste("./Data/GSA_hyper_",contrast,"_",input$gscgroup,"_",input$GSAcutoff,".rds",sep=""))
                }
                else {
                    
                    dt <- dxp$resTable[[contrast]][complete.cases(dxp$resTable[[contrast]]),]
                    gsah <- runGSAhyper(genes = dt[,1], 
                                        pvalues = dt[,5], 
                                        pcutoff = input$GSAcutoff, 
                                        gsc = mGsc, 
                                        gsSizeLim=c(1,Inf), 
                                        adjMethod="fdr")
                    saveRDS(gsah,paste("./Data/GSA_hyper_",contrast,"_",input$gscgroup,"_",input$GSAcutoff,".rds",sep=""))
                    
                }
                
                
                #############
                gsaRes = gsah
                class = "non"
                significance=input$GSAsign
                adjusted=FALSE
                geneSets=NULL
                overlap=1
                lay=1 
                label="names"
                cexLabel=0.9
                ncharLabel=25
                cexLegend=1
                nodeSize=c(10,40)
                edgeWidth=c(1,15)
                edgeColor=NULL
                scoreColors=NULL
                
                test <- 1 # Which contrast? Currently only one allowed!
                
                #*********************************************
                # Check arguments:
                #*********************************************
                
                tmp <- try(pValue <- match.arg(class, c("distinct","mixed","non"), several.ok=FALSE), silent=TRUE)
                if(class(tmp) == "try-error") {
                    stop("argument class is not valid")
                }
                if(pValue == "non") {
                    direction <- "none"
                } else {
                    tmp <- try(direction <- match.arg(direction, c("up","down","both"), several.ok=FALSE), silent=TRUE)
                    if(class(tmp) == "try-error") {
                        stop("argument direction is not valid")
                    }
                }
                
                if(pValue == "non") {
                    maintext <- paste(pValue,"-directional",sep="")
                } else {
                    maintext <- paste(pValue,"-directional (",direction,")",sep="")
                }
                
                
                if(pValue == "mixed" & direction == "both") stop("for pValue='mixed', direction can not be 'both'")
                
                if(pValue == "distinct" & direction == "both") pValue <- "dirupdn"
                if(pValue == "distinct" & direction == "up") pValue <- "dirup"
                if(pValue == "distinct" & direction == "down") pValue <- "dirdn"
                if(pValue == "non") pValue <- "mix"
                if(pValue == "mixed" & direction == "up") pValue <- "subup"
                if(pValue == "mixed" & direction == "down") pValue <- "subdn"
                
                #if(test < 1 | test > ncol(gsaRes$nGenes)) stop("argument test is out of bounds")
                if(significance < 0 | significance >1) stop("argument significance is out of bounds")
                if(overlap <= 0) stop("argument overlap has to be at least 1")
                if(length(nodeSize) != 2) stop("argument nodeSize has to have length 2")
                if(length(edgeWidth) != 2) stop("argument edgeWidth has to have length 2")
                if(class(adjusted) != "logical") stop("argument adjusted has to be TRUE or FALSE")
                                
                #########################################################
                # Adds possibility to use output object from runGSAhyper:
                if(length(gsaRes) == 5) {
                    if(all(names(gsaRes) == c("pvalues","p.adj","resTab","contingencyTable","gsc"))) {
                        if(pValue != "mix") stop("When using result from runGSAhyper, only class='non' is allowed")
                        gsaRes$pNonDirectional <- matrix(gsaRes$pvalues,ncol=1)
                        gsaRes$pAdjNonDirectional <- matrix(gsaRes$p.adj,ncol=1) 
                        gsaRes$geneSetStat <- "Fisher's exact test"
                        if(adjusted) {
                            maintext <- paste("p.adj<",significance,sep="")
                        } else {
                            maintext <- paste("p<",significance,sep="")
                        }
                    }
                }
                #########################################################
                
                
                #*********************************************
                # Prepare network:
                #*********************************************
                
                # Extract values:
                gsc <- gsaRes$gsc
                if(adjusted) { # ...if adjusted:
                    if(pValue == "dirup") pValues <- gsaRes$pAdjDistinctDirUp[,test]
                    if(pValue == "dirdn") pValues <- gsaRes$pAdjDistinctDirDn[,test]
                    if(pValue == "dirupdn") {
                        pValues <- apply(abs(cbind(gsaRes$pAdjDistinctDirUp[,test],gsaRes$pAdjDistinctDirDn[,test])),1,min,na.rm=TRUE)
                        tmp     <- apply(abs(cbind(gsaRes$pAdjDistinctDirUp[,test],gsaRes$pAdjDistinctDirDn[,test])),1,which.min)==2
                        pValues[pValues == 0] <- 1e-100
                        pValues[tmp] <- -pValues[tmp]
                    }
                    if(pValue == "mix") pValues <- gsaRes$pAdjNonDirectional[,test]
                    if(pValue == "subup") pValues <- gsaRes$pAdjMixedDirUp[,test]
                    if(pValue == "subdn") pValues <- gsaRes$pAdjMixedDirDn[,test]
                    
                } else { # ...if un-adjusted:
                    if(pValue == "dirup") pValues <- gsaRes$pDistinctDirUp[,test]
                    if(pValue == "dirdn") pValues <- gsaRes$pDistinctDirDn[,test]
                    if(pValue == "dirupdn") {
                        pValues <- apply(abs(cbind(gsaRes$pDistinctDirUp[,test],gsaRes$pDistinctDirDn[,test])),1,min,na.rm=TRUE)
                        tmp     <- apply(abs(cbind(gsaRes$pDistinctDirUp[,test],gsaRes$pDistinctDirDn[,test])),1,which.min)==2
                        pValues[pValues == 0] <- 1e-100
                        pValues[tmp] <- -pValues[tmp]
                    }
                    if(pValue == "mix") pValues <- gsaRes$pNonDirectional[,test]
                    if(pValue == "subup") pValues <- gsaRes$pMixedDirUp[,test]
                    if(pValue == "subdn") pValues <- gsaRes$pMixedDirDn[,test]
                }
                
                # Get gene set names:
                geneSetNames <- names(gsc)
                
                # Alt. 1. Select user defined gene sets:
                if(!is.null(geneSets)) {
                    indSignificant <- which(geneSetNames %in% geneSets)
                                        
                    # Alt. 2. Select significant gene sets:   
                } else {
                    indSignificant <- which(abs(pValues) < significance)
                }
                
                # Check if at least two gene sets are selected:
                if(length(indSignificant) < 3) stop("less than three gene sets were selected, can not plot (tip: adjust the significance cutoff)")
                
                # Get the relevant p-values:
                pSignificant <- pValues[indSignificant]
                
                # Generate overlap matrix for significant gene sets:
                overlapMat <- matrix(nrow=length(indSignificant),ncol=length(indSignificant))
                for(i in 1:nrow(overlapMat)) {
                    for(j in i:ncol(overlapMat)) {
                        tmp <- sum(gsc[[indSignificant[i]]] %in% gsc[[indSignificant[j]]])
                        overlapMat[i,j] <- tmp 
                        overlapMat[j,i] <- tmp
                    }
                }
                
                # Gene set size:
                gsSize <- diag(overlapMat)
                
                # Remove gene sets with only small overlap (defined by argument 'overlap'):
                overlapMat[overlapMat < overlap] <- 0
                
                # Adjecency matrix:
                adjMat <- overlapMat > 0
                
                # Create igraph object:
                tmp <- adjMat
                diag(tmp) <- 0 # For some reason diag=FALSE below stopped working...
                if(all(!tmp)) stop("no overlap between gene sets found, try to decrease argument overlap or increase argument significance")
                g <- graph.adjacency(tmp, mode="undirected", diag=FALSE)
                
                
                #*********************************************
                # Set plotting parameters:
                #*********************************************
                
                # Edge width, according to shared genes:
                edgeOverlap <- rep(NA,ecount(g))
                for(iEdge in 1:ecount(g)) {
                    tmp <- get.edge(g,iEdge)
                    edgeOverlap[iEdge] <- overlapMat[tmp[1],tmp[2]]
                }
                eWidth <- (edgeOverlap-min(edgeOverlap))/(max(edgeOverlap)-min(edgeOverlap))*(edgeWidth[2]-edgeWidth[1])+edgeWidth[1]
                
                # Edge color:
                if(is.null(edgeColor)) edgeColor = c("gray80","gray70","gray60")
                tmp <- seq(min(edgeOverlap), max(edgeOverlap), length.out=length(edgeColor)+1)
                eColor <- rep(edgeColor[1],ecount(g))
                for(i in 2:length(edgeColor)) {
                    eColor[edgeOverlap > tmp[i]] <- edgeColor[i]  
                }
                
                # Node size, according to number of genes:
                vSize <- (gsSize-min(gsSize))/(max(gsSize)-min(gsSize))*(nodeSize[2]-nodeSize[1])+nodeSize[1] 
                
                # Node color:
                if(pValue == "dirupdn") {
                    if(is.null(scoreColors)) {
                        tmp1 <- c('red','tomato','mistyrose')
                        tmp2 <- c('blue','cornflowerblue','azure')
                        gradColorsUp <- colorRampPalette(tmp1,interpolate="linear")(100)
                        gradColorsDn <- colorRampPalette(tmp2,interpolate="linear")(100)
                    } else {
                        if(length(scoreColors)%%2 != 0 | length(scoreColors) < 4) stop("argument scoreColors should contain at least four and an even number of colors")  
                        tmp1 <- scoreColors[1:(length(scoreColors)/2)]
                        tmp2 <- scoreColors[(length(scoreColors)/2+1):length(scoreColors)]
                        gradColorsUp <- colorRampPalette(tmp1,interpolate="linear")(100)
                        gradColorsDn <- colorRampPalette(tmp2,interpolate="linear")(100)
                    }
                    vColor <- rep(NA,length(pSignificant))
                    tmp <- pSignificant[pSignificant > 0 & !is.na(pSignificant)]
                    vColor[pSignificant > 0 & !is.na(pSignificant)] <- gradColorsUp[floor(tmp*99/max(tmp))+1]
                    tmp <- abs(pSignificant[pSignificant < 0 & !is.na(pSignificant)])
                    vColor[pSignificant < 0 & !is.na(pSignificant)] <- gradColorsDn[floor(tmp*99/max(tmp))+1]
                } else {
                    if(is.null(scoreColors)) {
                        tmp <- c('red','orange','yellow')
                    } else {
                        if(length(scoreColors) < 2) stop("argument scoreColors should contain at least two colors")
                        tmp <- scoreColors
                    }
                    gradColors <- colorRampPalette(tmp,interpolate="linear")(100)
                    if(max(pSignificant,na.rm=TRUE)==0) {
                        vColor <- gradColors[rep(50,length(pSignificant))]   
                    } else {
                        vColor <- gradColors[floor(pSignificant*99/max(pSignificant,na.rm=TRUE))+1]  
                    }
                } 
                vColor[is.na(vColor)] <- "#CCCCCC"
                
                # Node labels:
                tmp <- try(label <- match.arg(label, c("names","numbers","numbersAndSizes","namesAndSizes"), several.ok=FALSE), silent=TRUE)
                if(class(tmp) == "try-error") {
                    stop("argument label has to be set to either 'names' or 'numbers'")
                }
                tmp <- names(gsc)[indSignificant]
                for(i in 1:length(tmp)) {
                    if(nchar(tmp[i])>ncharLabel) tmp[i] <- paste(substr(tmp[i],1,ncharLabel),"...",sep="")
                }
                if(label == "names") vLabels <- tmp
                else if(label == "numbers") vLabels <- 1:length(indSignificant)
                else if(label == "numbersAndSizes") vLabels <- paste(1:length(indSignificant)," (",gsSize,")",sep="")
                else if(label == "namesAndSizes") vLabels <- paste(tmp," (",gsSize,")",sep="")
                
                
                #*********************************************
                # Set plotting layout:
                #*********************************************
                
                if(class(lay) == "numeric") {
                    
                    # Edge weight, inverse to number of edges and size (of the larger of the node pairs)
                    # used by Fruchterman Reingold algorithm:
                    eWeight <- rep(NA,ecount(g))
                    for(iEdge in 1:ecount(g)) {
                        tmp <- get.edge(g,iEdge)
                        tmp1 <- gsSize[tmp[1]]
                        tmp2 <- gsSize[tmp[2]]
                        if(tmp1 > tmp2) {
                            eWeight[iEdge] <- 1/(tmp1*sum(adjMat[,tmp[1]]))
                        } else {
                            eWeight[iEdge] <- 1/(tmp2*sum(adjMat[,tmp[2]]))
                        }
                    }
                    if(min(eWeight) != max(eWeight)) {
                        eWeight <- (eWeight-min(eWeight))/(max(eWeight)-min(eWeight))*(100-0)+0
                    } else {
                        eWeight <- rep(50,length(eWeight))
                    }
                    
                    # Predefined layouts:
                    if(lay == 1) lay <- layout.fruchterman.reingold(g,area=vcount(g)^4,repulserad=vcount(g)^5, weights=eWeight)      
                    else if(lay == 2) lay <- layout.fruchterman.reingold(g,weights=eWeight)
                    else if(lay == 3) lay <- layout.spring(g, repulse=TRUE, repeqdis=0.5, mass=2, equil=2)
                    else if(lay == 4) lay <- layout.graphopt(g, charge=vSize)
                    else if(lay == 5) lay <- layout.circle(g)
                    
                    # User defined layout function:   
                } else if(class(lay) == "function"){
                    lay <- lay(g)  
                }
                
                
                #*********************************************
                # Plot:
                #*********************************************
                
                if(pValue == "dirupdn") {
                    layout(matrix(c(1,4,2,4,3,3),ncol=3), widths=c(1,1,5), heights=c(2,3))
                    par(cex.axis=cexLegend)
                    maColorBar(seq(0,max(abs(pSignificant),na.rm=TRUE),by=max(abs(pSignificant),na.rm=TRUE)/100),FALSE,gradColorsUp, 
                               k=7, main="p-values up", cex.main=cexLegend)
                    maColorBar(seq(0,max(abs(pSignificant),na.rm=TRUE),by=max(abs(pSignificant),na.rm=TRUE)/100),FALSE,gradColorsDn, 
                               k=7, main="p-values down", cex.main=cexLegend)
                    
                } else {
                    layout(matrix(c(1,3,2,2),ncol=2), widths=c(1,5), heights=c(1,1))
                    maColorBar(seq(0,max(pSignificant,na.rm=TRUE),by=max(pSignificant,na.rm=TRUE)/100),FALSE,gradColors, 
                               k=7, main="p-values", cex.main=cexLegend)   
                }   
                main <- paste(gsaRes$geneSetStat,", ",maintext,sep="")
                plot(g, vertex.size=vSize, vertex.color=vColor, vertex.label=vLabels, vertex.label.family="sans",
                     vertex.label.color="black", vertex.label.cex=cexLabel, edge.width=eWidth, 
                     edge.color=eColor, layout=lay, main=main)
                plot.new()
                par(usr=c(-1,1,-1,1))
                text(0,0, cex=cexLegend,
                     labels=paste("Max node size:\n",max(gsSize)," genes\n\nMin node size:\n",min(gsSize)," genes\n\n",
                                  "Max edge width:\n",max(edgeOverlap)," genes\n\nMin edge width:\n",min(edgeOverlap),
                                  " genes",sep=""), font=1, adj=0.5)
                
                
            } else {
                plot(NULL,main="Network plot", xlim = c(0,20), ylim=c(0,0.21))
                legend('topleft',
                       ncol=2,
                       legend='Please select two different contrasts.')
            }
        })
    }
)