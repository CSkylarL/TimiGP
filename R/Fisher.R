##' Find marker pairs associated with response to therapy
##' 
##' Perform fisher test (two sided) to look for marker pair associated with response to therapy.
##'
##' @param mps a matrix of Marker Pair Score
##' @param info a data.frame in which the 1st column contains 
##' two groups on response to therapy, "Responder == 1" & "Non-Responder == 0" 
##' and row name is patient ID.
##' @param p.adj p.adjust.methods. 
##' One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"	
##' The default value is "BH".
##' @param parallel a logical value. If TRUE, enable parallel computing. 
##' The default value if FALSE
##' @param core a numeric value shows the number of cores for parallel execution. 
##' The default value is 1.
##' @return  Result of fisher test
##' (Rowname: Marker Pair,PV: P-Value,QV: Adjust P-value). 
##' @import dplyr
##' @import foreach
##' @import doParallel
##' @export 
##' @examples 
##' \dontrun{ 
##' # See example07 for more details
##' # Load the example data
##' data(TNBCaPD1info)
##' data(TNBCaPD1rna)
##' # Load the example cell type annotation
##' data("CellType_TNBC_aPDL1")
##' geneset <- CellType_TNBC_aPDL1
##' marker <- unique(geneset$Gene)
##' # Preprocess
##' info <- TNBCaPD1info
##' rna <- TimiPrePropress(marker = marker, cohort = rownames(info),
##'                        log = F, GMNorm = T, rna = TNBCaPD1rna)
##' # Generate marker pair score
##' mps <- TimiGenePair(rna)
##' # Use the function TimiFisher 
##' res <- TimiFisher(mps = mps,info = info,p.adj = "BH")
##' # Assign results
##' mps <- res$mps
##' fisher_res <- res$fisher_res
##' }
##' @author Chenyang Skylar Li

TimiFisher <-  function(mps = NULL,
                        info = NULL,
                        p.adj = "BH",
                        parallel = FALSE,
                        core = 1){
  # examine required parameters
  if (is.null(mps)){
    stop('The parameter "mps" is required. ')
  }
  
  if (is.null(info)){
    stop('The parameter "info" is required. ')
  }
  
  if (is.null(p.adj)){
    stop('The parameter "p.adj" is required. Please choose one of c("holm", "hochberg", "hommel", "bonferroni","BH", "BY","fdr", "none")')
  }
  
  
  if (all((colnames(mps) == rownames(info)) != TRUE) ) {
    comSam <- intersect(row.names(info), colnames(mps))
    if (length(comSam) > 0) {
      mps <- mps[,comSam]
      info <- info[comSam,]
    } else {
      stop('Samples are different between Gene pair score and clinical info')
    }
    if (all((colnames(mps) == rownames(info)) != TRUE) ) {
      stop('Samples are different between Gene pair score and clinical info')
    }
    
  }
  
  # fisher function ------------------------------------------------------------
  Myfun <- function(mytag = NULL,
                    info = NULL,
                    logic = T # T: 0 or 1; F: numeric
  ){
    if (logic == T){
      mytag <- as.factor(mytag)
      
      # Perform fisher test
      tmp <- as.data.frame(cbind(info, mytag))
      stat <- as.matrix(table(tmp))
      stat <-t(stat[c("1","0"),c("1","0")])
      test <- fisher.test(stat,alternative = "two.sided")
      res <- c(  test$estimate, test$p.value)
      names(res) <- c("OR","PV")
      return (res)
      
    } else if (logic == F) {
      
      stop("The function has not been relesed and please stay tuned. 
           For the current version, please set `cont = FALSE` in `TimiGenePair()`.")
      
    }
    
    
  }
  # perform fisher test /logistical regression ---------------------------------
  # examine the mps is generated from logical relation or continuous relation
  
  if (parallel == F) { # No parallel, use apply
    message("------------- Calculating -------------") 
    message("Please be patient. This step usually takes approximately 1 minutes.")  
    message("The computation time may be increased if there are a larger number of cell type markers (e.g., > 1000 markers for > 30 mins).")
    message('To speed up, you can enable parallel computing by setting "parallel = T" and the number of cores (core).')
    if (is.logical(mps) == T) {
      
      # perform fisher test <- logical ------------------------------------------------------
      Myres <- apply(mps * 1, 1, function (x) 
        Myfun(mytag = x,info = info, logic = T))
      
    } else if (is.numeric(mps) == T) {
      # perform logistical regression <- continuous ---------------------------------------------------
      Myres <- apply(mps, 1, function (x) 
        Myfun(mytag = x,info = info, logic = F))
      
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
      
      # perform fisher test <- logical ------------------------------------------------------
      Myres <- foreach (i = 1:nrow(mps), .combine=cbind) %dopar% {
        Myfun(mytag = mps[i,]*1,info = info, logic = T)
      }
      
    } else if (is.numeric(mps) == T) {
      # perform logistical regression <- continuous ---------------------------------------------------
      Myres <- foreach (i = 1:nrow(mps), .combine=cbind) %dopar% {
        Myfun(mytag = mps[i,]*1,info = info, logic = F)
      }
      
    } else {
      stop('The matrix are neither logical nor numeric values')
    }
    
  } 
  
  Myres <- Myres %>%
    t() %>% 
    data.frame() %>%
    mutate(QV = p.adjust(PV, method=p.adj))
  row.names(Myres) <- row.names(mps)
  
  
  # change pair direction of fisher result: ------------------------------------
  # all OR > 1 == MPS score associated with responders
  # reverse OR < 1
  se <- which(Myres$OR < 1)
  
  message("The direction of ",length(se)," marker pairs were reversed")
  
  xx <- row.names(Myres[se,])
  xx <- unlist(strsplit(xx, "_"))
  xx1 <- xx[seq(1,length(xx),2)]
  xx2 <- xx[seq(2,length(xx),2)]
  xx <- paste(xx2, xx1, sep="_")
  
  Myres$OR[se] <- 1/Myres$OR[se]
  rownames(Myres)[se] <- xx
  
  # change pair direction and reverse value of marker pair score ---------------
  rownames(mps)[se] <- xx
  
  if (is.logical(mps) == T) {
    mps[se,] <- !mps[se,]
  } else if (is.numeric(mps) == T) {
    mps[se,] <- -mps[se,]
  }  else {
    stop('The matrix are neither logical nor numeric values')
  }
  
  Myres <-Myres[order(Myres[,2]),]
  res <- list("mps" = mps,"fisher_res" = Myres)
  
  return(res)
  
}


