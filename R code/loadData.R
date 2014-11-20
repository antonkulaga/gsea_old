##############
# unfortunatelly, the original Piano is not fully suitable for my case
# so I re-used Piano scripts making them less generic yet fitting my data better

# here I modify Piano:::loadMAdata()


loadData <- function(datadir=getwd(), 
                     countfile,
                     setup, 
                     verbose=TRUE) {
    
    # Verbose function:
    .verb <- function(mes, verbose) {
        if(verbose == TRUE) {
            message(mes)
        }
    }
    
	if (file.exists("./Data/datafile.rds")) {
		l_arrayData <- readRDS("./Data/datafile.rds")
	} 
	else {
		    # load data I do no quality checks for now assuming that everything is OK
		countTable <- read.csv(paste(getwd(),countfile,sep=""), header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE)
		setup <- read.csv(paste(getwd(),setup,sep=""), header = TRUE, row.names=1, sep = "",  dec = ",", fill = TRUE, as.is=TRUE)
		#
		factorAnnot <- setup[,1]
		if(dim(setup)[2] > 1) {
			for(i in 1:(dim(setup)[2]-1)) {
				factorAnnot <- paste(factorAnnot,setup[,i+1], sep="_")
			}
		}
		condition <- factor(factorAnnot)

		# normalize data
		cds <- newCountDataSet(countTable, condition)
		cds <- estimateSizeFactors( cds )
		cds <- estimateDispersions( cds )
		vsd <- getVarianceStabilizedData(cds)

		# Construct ArrayData object as return:
		l_arrayData <- list(dataNorm=as.data.frame(vsd), setup=setup, dataCount = cds)
		class(l_arrayData) = "ArrayData"
		
		saveRDS(l_arratData,"./Data/datafile.rds")
	}
	
    return(l_arrayData)
}