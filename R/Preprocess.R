##' Check Clinical Info
##' 
##' Check Event and Time-to-event; remove NAs and time-to-event <=0
##' 
##'
##' @param info a data.frame in which the 1st column is event 
##' and the 2nd column is days-to-event.
##' @importFrom survival Surv
##' @return filtered clinical info  
##' @export 
##' @examples
##' \dontrun{
##'   data(SKCM06info)
##'   dim(SKCM06info)
##'   info <- TimiCheckEvent(SKCM06info)
##'   dim(info) 
##' }
##' @author Chenyang Skylar Li

TimiCheckEvent <-  function(info = NULL){
  
  # examine required parameters ------------------------------------------------
  if (is.null(info)){
    stop('The parameter "info" is required. ')
  }

  # Remove NAs -----------------------------------------------------------------
  se1 <- which(is.na(info[1]))
  se2 <- which(is.na(info[2]))
  y=Surv(as.numeric(info[,2]), as.numeric( info[,1]))
  # Remove time <=0-------------------------------------------------------------
  se3 <-  which(y[, "time"] <= 0) 
  se4 <- which(is.na(y))
  se <- unique(c(se1,se2,se3,se4))
  if (length(se) != 0){
    info <- info[-se,]
  }
  
  message(length(se), " individuals were filtered out due to NAs or time <= 0")
  return(info)
}




##' Preprocess Transctiptomic profile
##' 
##' Given the transcriptomic profile, marker gene list, clinical information,
##' this function will return the marker gene expression in selected cohort,
##' after log transformation and  gene wise median normalization 
##' 
##'
##' @param marker a vector of marker gene list.
##' @param rna a data.frame of transcriptomic profile, 
##' in which row names are genes and colnames are individuals.
##' @param cohort a vector of selected individuals.
##' @param log logical value: 
##' if TRUE, log transformation will be performed; The default is TRUE.
##' @param GMNorm logical value: 
##' if TRUE, gene wise median normalization will be performed; The default is TRUE.
##' @return Preprocessed transctiptomic profile 
##' @export 
##' @examples
##' \dontrun{
##'   data("SKCM06info")
##'   data("SKCM06rna")
##'   data("Immune_Marker_n1293")
##'   dim(SKCM06info)
##'   info <- TimiCheckEvent(SKCM06info)
##'   dim(info)
##'   dim(SKCM06rna)
##'   rna <- TimiPrePropress(gene = Immune_Marker_n1326,rna = SKCM06rna,cohort = rownames(info), log = T,GMNorm = T)
##'   dim(rna)
##' }
##' @author Chenyang Skylar Li
TimiPrePropress <- function(marker,
                            rna = NULL,
                            cohort,
                            log=TRUE,
                            GMNorm=TRUE){
  
  # examine required parameters-------------------------------------------------
  
  if (is.null(rna)){
    stop('The parameter "rna" is required. ')
  }
  
  
  # Extract cohort--------------------------------------------------------------
  if (is.null(cohort)){
    warning('No selected cohort. Use the cohort in transcriptiomic profile')
  } else {
    se1 <- which(colnames(rna) %in% cohort)
    if (length(se1) == 0) {
      warning('No individual was found in the transcriptome. Use the enrire cohort.')
    } else if (length(se1) > 0 & length(se1) < 5) {
      warning('No.Individuals < 5 was found in the transcriptome. Use the enrire cohort.')
    } else {
      rna <- rna[,se1]
    }
  }
 
  
  # Extract marker--------------------------------------------------------------
  if (is.null(cohort)){
    warning('No markers. Use the gene in transcriptiomic profile')
  } else {
    se2 <- which(rownames(rna) %in% marker)
    if (length(se2) == 0) {
      stop('No markers was found in the transcriptome. Please check marker ID')
    } else {
      rna <- rna[se2,]
    }
  }
  
  
  # remove genes with low expressions-------------------------------------------
  xx <- apply(rna>0, 1, sum)
  se3 <- which(xx>=ncol(rna)*0.5) 
  if (length(se3) == 0){
    stop('All markers were  in a low expression level')
  }
  rna <- rna[se3,]
  message(nrow(rna)-length(se3)," markers with low expressions were filtered out")
  
  # log transformatio-----------------------------------------------------------
  if(log == TRUE){
    rna <- log10(rna+1)
  }
  
  # gene wise median normalization ---------------------------------------------
  if(GMNorm == TRUE){
    xx <- apply(rna, 1, median)
    rna <- rna-xx
  }
  idx <- order(row.names(rna))
  rna <- rna[idx,]
  return(rna)
}
