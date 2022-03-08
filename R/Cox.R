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
##' One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"	
##' @return 1. Result of cox regression
##' (HR: Hazard.Ratio,PV: P-Value,QV: Adjust P-value). 
##' 2. modified marker pair score
##' @import survival
##' @export 
##' @examples
##' \dontrun{
##'   data("SKCM06info")
##'   data("SKCM06rna")
##'   data("Immune_Marker_n1326")
##'   info <- TimiCheckEvent(SKCM06info)
##'   rna <- TimiPrePropress(gene = Immune_Marker_n1326,rna = SKCM06rna,cohort = rownames(info))
##'   mps <- TimiGenePair(rna)
##'   dim(mps)
##'   # TimiCOX 
##'   res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
##'   mps <- res[[1]]
##'   cox_res <- res[[2]]
##' }
##' @author Chenyang Skylar Li

TimiCOX <-  function(mps = NULL,
                     info = NULL,
                     p.adj = c("holm", "hochberg", "hommel", "bonferroni", 
                               "BH", "BY","fdr", "none")){
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
    stop('Samples are different between Gene pair score and clinical info')
  }
  # perform cox
  pval <- hr  <- rep(0, nrow(mps))
    for(k in 1:nrow(mps)){
      cat("\rCOX:", k,"/",nrow(mps))
      mytag <- ifelse(mps[k,], 1, 0)
      xx <- as.data.frame(cbind(info, mytag))
      mycox <- coxph(Surv(as.numeric(info[,2]), as.numeric( info[,1]))~mytag, xx) 
      mycox <- summary(mycox)
      pval[k] <- mycox$coefficients[5]
      tmp <- mycox$conf.int
      hr[k] <- tmp[1]
    }
  cat("\n")
  
  QV <- p.adjust(pval, method=p.adj)
  
  cox_res <- data.frame(HR=hr, PV=pval, QV=QV)
  row.names(cox_res) <- row.names(mps)
  # change pair direction of cox result
  se <- which(cox_res$HR>1)
  message("The direction of ",length(se)," marker pairs were reversed")
  
  xx <- row.names(cox_res[se,])
  xx <- unlist(strsplit(xx, "_"))
  xx1 <- xx[seq(1,length(xx),2)]
  xx2 <- xx[seq(2,length(xx),2)]
  xx <- paste(xx2, xx1, sep="_")
  
  cox_res$HR[se] <- 1/cox_res$HR[se]
  rownames(cox_res)[se] <- xx
  
  # change pair direction and reverse value of marker pair score
  rownames(mps)[se] <- xx
  mps[se,] <- !mps[se,]
  
  cox_res <-cox_res[order(cox_res[,2]),]
  res <- list(mps,cox_res)
  
  return(res)
  
}


