##############
# unfortunatelly, the original Piano is not fully suitable for my case
# so I re-used Piano scripts making them less generic yet fitting my data better

# here I modify Piano:::runQS(). Mostly because I do not have dataRaw



# removed the following functions
# - saving plots to files
# - function for rna degradation based on AffyRNAdeg() -- no dataRaw 
#       read more in http://watson.nci.nih.gov/bioc_mirror/packages/2.13/bioc/manuals/affy/man/affy.pdf
# - function for affyPLM: fitPLM(), RLE(), NUSE(), -- no dataRaw
#       read more in http://www.bioconductor.org/packages/release/bioc/manuals/affyPLM/man/affyPLM.pdf 
#       http://www.bioconductor.org/packages/release/bioc/vignettes/affyPLM/inst/doc/QualityAssess.pdf

runQC <- function(arrayData, 
                  hist=TRUE, boxplot=TRUE, pca=TRUE, 
                  colorFactor=1,
                  colors=c("red","green","blue","yellow","orange",
                           "purple","tan","cyan","gray60","black", "white"),
                  verbose=TRUE) {
    
    # Verbose function:
    .verb <- function(message, verbose) {
        if(verbose == TRUE) {
            message(message)
        }
    }
    
    # First the functions for the various QCs (the code that calls the functions is below):
    
    # Function for distribution:
    .runHist <- function(arrayData) {
        # Plot file:
        .verb("Generating normalized data distribution plot...", verbose)
        dev.new()
        plot(density(arrayData$dataNorm[,1]), main="Histogram of normalized data",
             xlab=expression(log[2] ~ intensity),ylab="Density")
        for(i in 2:dim(arrayData$dataNorm)[2]){
            lines(density(arrayData$dataNorm[,i]))
        }
        .verb("...done", verbose)
    }
    
    # Function for boxplot:
    .runBoxplot <- function(arrayData) {
        # Plot file:
        .verb("Generating normalized data boxplot...", verbose)
        dev.new()
        par(mar=c(10, 4, 4, 2))
        boxplot(arrayData$dataNorm, main="Boxplot of normalized data",
                ylab=expression(log[2] ~ intensity), las=2)
        .verb("...done", verbose)
    }
    
    # Function for PCA:
    .runPCA <- function(arrayData, colorFactor, colors) {
        .verb("Generating PCA...", verbose)
        dataForPCA <- arrayData$dataNorm
        # Centralize data
        dataForPCA <- dataForPCA - rowMeans(dataForPCA)
        
        # PCA
        dataPrcomp <- prcomp(dataForPCA)
        # Order by pc3 (to plot in size order)
        pc3 <- dataPrcomp$rotation[,3]
        sizeOrder <- order(pc3, decreasing=TRUE)
        pc3 <- pc3[sizeOrder]
        pc1 <- dataPrcomp$rotation[,1]
        pc1 <- pc1[sizeOrder]
        pc2 <- dataPrcomp$rotation[,2]
        pc2 <- pc2[sizeOrder]
        # Scale pc3 size to better interval
        MinSize <- 1.5
        MaxSize <- 4.5
        pc3Size <- MinSize + (MaxSize - MinSize) * (pc3 - min(pc3))/(max(pc3) - min(pc3))
        # Coloring
        Colors <- colors
        mainFactors <- unique(arrayData$setup[,colorFactor])
        colorKey <- arrayData$setup[,colorFactor]
        colorKeyFactors <- mainFactors
        for(i in 1:length(mainFactors)) {
            colorKey[colorKey == mainFactors[i]] <- Colors[i]
            colorKeyFactors[colorKeyFactors == mainFactors[i]] <- Colors[i]
        }
        orderIndex <- NA
        for(i in 1:length(names(pc1))) {
            orderIndex[i] <- which(rownames(arrayData$setup) == names(pc1)[i])
        }
        colorKey <- colorKey[orderIndex]
        .verb("...done", verbose)
        
        
        # Plot variance
        dev.new()
        mp <- barplot(summary(dataPrcomp)$importance[2,],
                      names.arg=colnames(dataPrcomp$rotation),
                      las=3,
                      ylab="Proportion of variance",
                      ylim=c(0,1.25*summary(dataPrcomp)$importance[2,][1]),
                      main="PC importance")
        text(x=cbind(mp,summary(dataPrcomp)$importance[2,]+summary(dataPrcomp)$importance[2,1]*0.15),
             labels=paste(round(summary(dataPrcomp)$importance[3,]*100,1),"%",sep=""),
             srt=90)
        
        # PCA plot
        dev.new()
        layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(7,1), widths=c(1,1))
        par(mar=c(4,4,2,1))
        plot(cbind(pc1,pc2), 
             pch=21, 
             col="black", 
             bg=colorKey, 
             cex=pc3Size, 
             main="PCA plot", 
             xlab="PC1", 
             ylab="PC2")
        mtext("PC3 (dot size)", side=4)
        text(cbind(pc1,pc2), labels=names(pc1), cex=0.8)
        
        par(mar=c(0,3,1,1))
        plot.new()
        title(main="Legend")
        legend("bottomleft", legend=mainFactors,fill=colorKeyFactors,ncol=length(mainFactors)/2+1)

    }
    
    
    # Below is the code that runs the selected functions:
    if(class(arrayData) != "ArrayData") {
        stop("argument arrayData is not of class ArrayData")
    }
    
    # Run the selected QCs:
    if(hist == TRUE) {
        .runHist(arrayData)
    }
    if(boxplot == TRUE) {
        .runBoxplot(arrayData)
    }
    if(pca == TRUE) {
        .runPCA(arrayData, colorFactor=colorFactor, colors=colors)
    }
    
}
