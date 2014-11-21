loadGSC <- function(gdata, group, addInfo) {

   # Initial argument checks:
   if(missing(addInfo)) {
      addUserInfo <- "skip"
      addInfo <- "none"
   } else {
      addUserInfo <- "yes"
   }
   
  
   fbids <- as.character(gdata$resTable[[1]][,1])
   
   if (file.exists(paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))) {

       res <- readRDS(paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))
   
    } else {
       
       
       ensembl = useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
       x <- getBM(attributes=c("flybase_transcript_id", "name_1006", "namespace_1003"),
                  # attributes=c("flybase_transcript_id", "name_1006")
                  values=fbids, 
                  mart=ensembl)
       genes2genesets <- unique(x[which(complete.cases(x) & x[,3]==group),1:2])
       
       tmp <- try(gsc <- as.data.frame(genes2genesets, stringsAsFactors=FALSE), silent=TRUE)
       if(class(tmp) == "try-error") {
           stop("argument gdata could not be converted into a data.frame")
       }
       
       
       
       # Get rid of factors:
       for(i in 1:ncol(gsc)) {
           gsc[,i] <- as.character(gsc[,i])
       }
       
       # Check gsc for two columns:
       if(ncol(gsc)!=2) {
           stop("argument gdata has to contain exactly two columns")  
       }
       
       # Remove redundant rows:
       tmp <- nrow(gsc)
       gsc <- unique(gsc)
       #info$redundantGS <- tmp - nrow(gsc)
       
       # Convert to list object:
       geneSets <- unique(gsc[,2])
       gscList <- list()
       for(iGeneSet in 1:length(geneSets)) {
           gscList[[iGeneSet]] <- gsc[gsc[,2] == geneSets[iGeneSet],1]
       }
       names(gscList) <- geneSets
       gsc <- gscList
       
       
       #***************************
       # AddInfo
       #***************************
       
       # Additional info as data.frame:
       if(addUserInfo == "yes") {
           tmp <- try(addInfo <- as.data.frame(addInfo, stringsAsFactors=FALSE), silent=TRUE)
           if(class(tmp) == "try-error") {
               stop("failed to convert additional info in argument 'addInfo' into a data.frame")
           }
       }
       
       if(class(addInfo) == "data.frame") {
           
           # Check for 2 columns:
           if(ncol(addInfo) != 2) stop("additional info in argument 'gdata' or 'addInfo' has to contain 2 columns")
           
           # Check addInfo correlation to gsc:
           tmp <- nrow(addInfo)
           addInfo <- unique(addInfo[addInfo[,1] %in% names(gsc),])
           #info$unmatchedAddInfo <- tmp - nrow(addInfo)
       } else {
           #info$unmatchedAddInfo <- 0     
       }
       
       #********************************
       # Return values:
       #********************************
       
       res <- list(gsc,addInfo)
       names(res) <- c("gsc","addInfo")
       class(res) <- "GSC"
       
       saveRDS(res,paste("./Data/dmelanogaster_getBM","_",group,".rds",sep=""))
       # listAttributes(mart)[c(27:32,37:50,73:74,91:100,123:130,414:419,1004:1009,1017:1024,1079:1085),]
   }


   return(res)
   
}