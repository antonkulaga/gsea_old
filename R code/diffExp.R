##############
# unfortunatelly, the original Piano is not fully suitable for my case
# so I re-used Piano scripts making them less generic yet fitting my data better

# here I modify Piano:::diffExp(). Mostly because annotation is not provided.

diffExp <- function(arrayData, contrasts, 
					fitMethod="ls", adjustMethod="fdr", significance=0.001, 
					plot=TRUE, heatmapCutoff=1e-10, volcanoFC=2,
					colors=c("red","green","blue","yellow","orange",
							 "purple","tan","cyan","gray60","black"),
					verbose=TRUE) {
	
	# Argument check:
	if(missing(contrasts)) stop("argument contrasts is not defined")
	if(!fitMethod %in% c("ls","robust","nbinom")) {
		stop("incorrect value of argument fitMethod")
	}
	if(!adjustMethod %in% c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none")) {
		stop("incorrect value of argument adjustMethod")
	}
	if(class(plot) == "logical") {
		if(plot) {
			venn <- heatmap <- volcano <- TRUE  
		} else {
			venn <- heatmap <- volcano <- FALSE
		}
	} else if(class(plot) == "character") {
		venn <- heatmap <- polarPlot <- volcano <- FALSE
		if("venn" %in% plot) venn <- TRUE
		if("heatmap" %in% plot) heatmap <- TRUE
		if("volcano" %in% plot) volcano <- TRUE
	} else {
		stop("argument plot has to be either TRUE, FALSE or a character string")
	}
	
	# Verbose function:
	.verb <- function(message, verbose) {
		if(verbose == TRUE) {
			message(message)
		}
	}
	
	# Get factors from setup and sort according to dataNorm columns
	factors <- extractFactors(arrayData)

	# Calculate p-values for each contrast:
	pValues <- NA
	foldChanges <- NA
	topTabList <- list()
        
    # run difference expression according to the method, nbinomTest is different to the rest
    if (fitMethod != "nbinom") {
        # Run lmFit
        .verb("Fitting linear models...", verbose)
        factors <- factor(factors$factors[,1])
        designMatrix <- model.matrix(~0+factors)
        colnames(designMatrix) <- levels(factors)
        dataForLimma <- arrayData$dataNorm
        fitLm <- lmFit(dataForLimma, design=designMatrix, method=fitMethod, maxit=200)
        
        # Run ebayes
        contrastMatrix <- makeContrasts(contrasts=contrasts, levels=levels(factors))
        fitContrasts <- contrasts.fit(fitLm,contrasts=contrastMatrix)
        fitContrasts <- eBayes(fitContrasts)
        .verb("...done", verbose)
        
        
        for(i in 1:length(contrasts)) {
            .verb(paste("Calculating p-values for contrast ",colnames(fitContrasts)[i],"...",sep=""), verbose)
            topTab <- topTable(fitContrasts,coef=i,adjust.method=adjustMethod,sort.by="none",number=nrow(dataForLimma))
            
            # For limma ersion compatability:
            if(!"ID" %in% colnames(topTab)) {
                topTab <- cbind(rownames(topTab),topTab)
                colnames(topTab)[1] <- "ID"
            }
            
            colnames(topTab)[1] <- c("ProbesetID") 
            
            # Check that topTab rows are consistent with arrayData$dataNorm rows:
            if(nrow(topTab)!=nrow(arrayData$dataNorm)) {
                stop("failed to order the genes correctly")
            }
            if(!all(topTab[,1]==rownames(arrayData$dataNorm))) {
                stop("failed to order the genes correctly")
            }
            
            topTabList[[i]] <- topTab
            pValues <- cbind(pValues,topTab$adj.P.Val)  # P-values for all contrasts, to be further used below
            foldChanges <- cbind(foldChanges,topTab$logFC)  # FC for all contrasts
            .verb("...done", verbose)
        }
        pValues <- as.data.frame(pValues[,2:ncol(pValues)],stringsAsFactors=FALSE)
        rownames(pValues) <- topTab[,1]
        colnames(pValues) <- colnames(fitContrasts)
        foldChanges <- as.data.frame(foldChanges[,2:ncol(foldChanges)],stringsAsFactors=FALSE)
        rownames(foldChanges) <- topTab[,1]
        colnames(foldChanges) <- colnames(fitContrasts)
    } else {

        for (i in 1:length(contrasts)) {
            a <- unlist(strsplit(contrasts[i], split=" - "))
            
            if (file.exists(paste("./Data/de_nbinom_",a[1],"_",a[2],".rds",sep=""))) {
                
                .verb(paste("Loading negative binom for ", contrasts[i], "...", sep=""), verbose)
                res<- readRDS(paste("./Data/de_nbinom_",a[1],"_",a[2],".rds",sep=""))
                
            }
            else {
                .verb(paste("Running negative binom for ", contrasts[i], "...", sep=""), verbose)
                res <- nbinomTest(arrayData$dataCount, a[1], a[2])
                saveRDS(res,paste("./Data/de_nbinom_",a[1],"_",a[2],".rds",sep=""))
            }
            
                pValues <- cbind(pValues,res$padj)  # P-values for all contrasts, to be further used below
                foldChanges <- cbind(foldChanges,res$log2FoldChange)  # FC for all contrasts
                topTabList[[i]] <- res[,c(1,6,2,3,7,8,4)] # full statistics
                colnames(topTabList[[i]]) <- c("ProbesetID", "logFC", "baseMean", "baseMeanA", "P.Value", "adj.P.Val", "baseMeanB")
                .verb(paste("Finished with negative binom for ", contrasts[i], "...", sep=""), verbose)
        }
        pValues <- as.data.frame(pValues[,2:ncol(pValues)],stringsAsFactors=FALSE)
        rownames(pValues) <- res[,1]
        colnames(pValues) <- contrasts
        foldChanges <- as.data.frame(foldChanges[,2:ncol(foldChanges)],stringsAsFactors=FALSE)
        rownames(foldChanges) <- res[,1]
        colnames(foldChanges) <- contrasts
        
    }
	
	# Heatmap
	if(heatmap == TRUE) {
	    .verb("Generating heatmap...", verbose)
	    genesInHeatmap <- vector()
	    for(i in 1:ncol(pValues)) {
	        genesInHeatmap <- c(genesInHeatmap, rownames(pValues)[pValues[,i] < heatmapCutoff])
	    }
	    if(length(unique(genesInHeatmap)) < 2) {
	        .verb("No genes were selected, change the heatmapCutoff. Omitting heatmap.", verbose)
	        # warning("No genes were selected, change the heatmapCutoff. Omitting heatmap.")
	    } else {
	        genesInHeatmap <- unique(genesInHeatmap)
	        heatmapMatrix <- as.matrix(arrayData$dataNorm[rownames(arrayData$dataNorm) %in% genesInHeatmap,])
	        if("annotation" %in% attributes(arrayData)$names) {
	            geneNamesHeatmap <- arrayData$annotation[rownames(arrayData$annotation) %in% genesInHeatmap,]
	            for(j in 1:nrow(heatmapMatrix)) {
	                tmp <- geneNamesHeatmap[rownames(geneNamesHeatmap) == rownames(heatmapMatrix)[j],1]
	                if(length(tmp) == 1) {
	                    rownames(heatmapMatrix)[j] <- tmp
	                }
	            }
	        }
	        dev.new()
	        heatmap.2(heatmapMatrix, trace="none", margins=c(10,12))
	        .verb("...done", verbose)
	    }
	}

	# Output
	names(topTabList) <- contrasts
	return(list(pValues=pValues,foldChanges=foldChanges,resTable=topTabList))
}