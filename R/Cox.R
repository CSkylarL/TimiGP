##' Find marker pairs associated with favorable prognosis
##' 
##' Perform univariate cox regression that fits each marker pair as a variable.
##' If a Pair A_B associated with poor prognosis(HR > 1), even not significant, 
##' it will be changed to B_A and reverse its value in the matrix of pair
##'
##' @param mps a matrix of Marker Pair Score
##' @param info a data.frame in which the 1st column is event 
##' and the 2nd column is days-to-event.
##' @param p.adj p.adjust.methods. 
##' One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
##' The default value is "BH".
##' @param parallel a logical value. If TRUE, enable parallel computing. 
##' The default value if FALSE
##' @param core a numeric value shows the number of cores for parallel execution. 
##' The default value is 1.
##' @return A list of two results: 
##' 1. Result of cox regression(cox_res)
##' (HR: Hazard.Ratio,PV: P-Value,QV: Adjust P-value). 
##' 2. modified marker pair score(mps)
##' @import survival
##' @import dplyr
##' @import foreach
##' @import doParallel
##' @export 
##' @examples
##' \dontrun{
##'   data("SKCM06info")
##'   data("SKCM06rna")
##'   data("Immune_Marker_n1293")
##'   info <- TimiCheckEvent(SKCM06info)
##'   rna <- TimiPrePropress(gene = Immune_Marker_n1293,rna = SKCM06rna,cohort = rownames(info))
##'   mps <- TimiGenePair(rna)
##'   dim(mps)
##'   # TimiCOX 
##'   res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
##'   mps <- res[["mps"]] 
##'   cox_res <- res[["cox_res"]]
##' }
##' @author Chenyang Skylar Li

TimiCOX <-  function(mps = NULL,
                     info = NULL,
                     p.adj = "BH",
                     parallel = FALSE,
                     core = 1){
  # examine required parameters-------------------------------------------------
  if (is.null(mps)){
    stop('The parameter "mps" is required. ')
  }
  
  if (is.null(info)){
    stop('The parameter "info" is required. ')
  }
  
  if (is.null(p.adj)){
    stop('The parameter "p.adj" is required. Please choose one of c("holm", "hochberg", "hommel", "bonferroni","BH", "BY","fdr", "none")')
  }
  
  
  if (all(colnames(mps) == rownames(info))  == FALSE){
    comSam <- intersect(row.names(info), colnames(mps))
    if (length(comSam) > 0) {
      mps <- mps[,comSam]
      info <- info[comSam,]
    } else {
      stop('Samples are different between Gene pair score and clinical info')
    }
    if (all(colnames(mps) == rownames(info))  == FALSE)  {
      stop('Samples are different between Gene pair score and clinical info')
    }
    
  }
  
  # cox function ---------------------------------------------------------------
  MyCox <- function(mytag = NULL,
                    info = NULL,
                    logic = T # T: 0 or 1; F: numeric
                    ){
    if (logic == T){
      mytag <- as.factor(mytag)
     
    } else if (logic == F) {
      mytag <- as.numeric(mytag)
    }
    
    tmp <- as.data.frame(cbind(info, mytag))
    test <- coxph(Surv(as.numeric(info[,2]),
                       as.numeric( info[,1]))~mytag, 
                  tmp) 
    test <- summary(test)
    res <- c( test$conf.int[1], test$coefficients[5])
    names(res) <- c("HR","PV")
    return (res)
  }
  # perform cox ----------------------------------------------------------------
  # examine the mps is generated from logical relation or continuous relation
  
  if (parallel == F) { # No parallel, use apply
    message("------------- Calculating -------------") 
    message("Please be patient. This step usually takes approximately 5 minutes.")  
    message("The computation time may be increased if there are a larger number of cell type markers (e.g., > 1000 markers for > 30 mins).")
    message('To speed up, you can enable parallel computing by setting "parallel = T" and the number of cores (core).')
    if (is.logical(mps) == T) {
      
      # perform cox logical ------------------------------------------------------
      cox_res <- apply(mps * 1, 1, function (x) 
        MyCox(mytag = x,info = info, logic = T))
      
    } else if (is.numeric(mps) == T) {
      # perform cox continuous ---------------------------------------------------
      cox_res <- apply(mps, 1, function (x) 
        MyCox(mytag = x,info = info, logic = F))
      
    } else {
      stop('The matrix are neither logical nor numeric values')
    }
  } else if (parallel == T) { # parallel, select core
    message("Enable parallel computing.")
    if (core == 1) {
      warning("Only 1 core specified. You may want to change core > 1")
    }
    
    message("------------- Calculating -------------") 
    message("Please be patient.")  
    message("The computation time may be increased if there are a larger number of cell type markers.")
    registerDoParallel(cores=core)

    if (is.logical(mps) == T) {
      
      # perform cox logical ------------------------------------------------------
      cox_res <- foreach (i = 1:nrow(mps), .combine=cbind) %dopar% {
        MyCox(mytag = mps[i,]*1,info = info, logic = T)
      }
      
    } else if (is.numeric(mps) == T) {
      # perform cox continuous ---------------------------------------------------
      cox_res <- foreach (i = 1:nrow(mps), .combine=cbind) %dopar% {
        MyCox(mytag = mps[i,]*1,info = info, logic = F)
      }
      
    } else {
      stop('The matrix are neither logical nor numeric values')
    }
    
  } 
  
  cox_res <- cox_res %>%
    t() %>% 
    data.frame() %>%
    mutate(QV = p.adjust(PV, method=p.adj))
  row.names(cox_res) <- row.names(mps)
  # change pair direction of cox result-----------------------------------------
  se <- which(cox_res$HR>1)
  message("The direction of ",length(se)," marker pairs were reversed")
  
  xx <- row.names(cox_res[se,])
  xx <- unlist(strsplit(xx, "_"))
  xx1 <- xx[seq(1,length(xx),2)]
  xx2 <- xx[seq(2,length(xx),2)]
  xx <- paste(xx2, xx1, sep="_")
  
  cox_res$HR[se] <- 1/cox_res$HR[se]
  rownames(cox_res)[se] <- xx
  
  # change pair direction and reverse value of marker pair score----------------
  rownames(mps)[se] <- xx
  
  if (is.logical(mps) == T) {
  mps[se,] <- !mps[se,]
  } else if (is.numeric(mps) == T) {
    mps[se,] <- -mps[se,]
  }  else {
    stop('The matrix are neither logical nor numeric values')
  }
  
  cox_res <-cox_res[order(cox_res[,2]),]
  res <- list("mps" = mps,"cox_res" = cox_res)
  
  return(res)
  
}


