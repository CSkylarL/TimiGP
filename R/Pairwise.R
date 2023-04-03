##' Capture logical/continuous relation between marker pairs
##' 
##' Default: Capture logical relation of any two marker pairs, 
##' and generate a matrix of Marker Pair Score:
##' 1 or TRUE = the expression of gene A > that of gene B, 
##' 0 or FALSE = the expression of gene A < that of gene B. 
##' \cr
##' Optional:  Capture continuous relation of any two marker pairs, 
##' and generate a matrix of Marker Pair Score:
##' the expression of gene A - that of gene B,
##'
##' @param rna a data.frame of preprocessed transcriptomic profile
##' @param cont a logical value. If TRUE, capture the continuous relation. 
##' The default is FALSE
##' @return a matrix of Marker Pair Score
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
##'   rna <- TimiPrePropress(gene = Immune_Marker_n1293,rna = SKCM06rna,cohort = rownames(info))
##'   # TimiGenePair
##'   mps <- TimiGenePair(rna)
##'   dim(mps)
##' }
##' @author Chenyang Skylar Li

TimiGenePair <-  function(rna = NULL,
                          cont = FALSE){
  
  # examine required parameters-------------------------------------------------
  if (is.null(rna)){
    stop('The parameter "rna" is required. ')
  }
  
  if (is.logical(cont) == F){
    stop('cont is a logical value. Please set it to TRUE or FALSE')
  }
  
  if (cont == T & (all(apply(rna,1,median) == 0) == T)) {
    stop('If you want to capture continuous relation between gene pairs,
         please do not use gene-wise median normalization. 
         Please set "GMNorm = F" in "TimiPrePropress" then re-run this function')
  }
  
  message("Generating marker pairs")
  ### S1: pairwise comparision -------------------------------------------------
  irg <- rownames(rna)
  irgp <- outer(irg, irg, paste, sep="_")
  if (cont == FALSE) {
    message("The default is to capture logical relation")
    myfun <- function(myvec){as.vector(outer(myvec, myvec, ">"))}
  } else if (cont == TRUE) {
    message("You choose to capture continuous relation")
    myfun <- function(myvec){as.vector(outer(myvec, myvec, "-"))}
  }
  
  irgp_res <- apply(rna, 2, myfun)
  row.names(irgp_res) <- irgp	 
  
  ### S2: Remove B_A, A_A pairsm, only keep A_B---------------------------------
  xx <- row.names(irgp_res)
  nn <- length(xx)
  tmp <- unlist(strsplit(xx, "_"))
  tmp1 <- tmp[(1:nn)*2-1]
  tmp2 <- tmp[(1:nn)*2-0]
  se <- which(tmp1>tmp2)
  irgp_res <- irgp_res[se,] 
  dim(irgp_res) 
  
  ### S3: filter GPs
  if (cont == FALSE) {
    ### choose GPs with at least 10% samples in both groups-------------------
    n.thr <- round(ncol(rna)*0.1)
    xx1 <- apply(irgp_res==1, 1, sum)
    xx2 <- apply(irgp_res==0, 1, sum)
    se <- which(xx1>=n.thr & xx2>=n.thr)
    
    message(nrow(irgp_res)-length(se)," marker pairs were filtered out")
    irgp_res <- irgp_res[se,]
    dim(irgp_res) 
    
    message(nrow(irgp_res)," marker pairs were produced")
  } else if (cont == TRUE) {
    ### choose GPs with difference > 0.1 across patients -----------------------
    min_range <- 0.1
    se <- apply(irgp_res, 1, function(x) { (max(x) - min(x)) > min_range }) 
    message(nrow(irgp_res)-length(se)," marker pairs were filtered out (CONT)")
    irgp_res <- irgp_res[se,]
    message(nrow(irgp_res)," marker pairs were produced (CONT)")
  }
    
  return(irgp_res)
}



##' Generate cell pair representing potential interaction using 
##' the corresponding marker pair annotation
##' 
##' Given any two different cell types,
##' Cell A has markers a1 and a2 and cell B has markers b1 and b2.
##' Cell pair A_B includes marker pairs: a1_b1, a1_b2, a2_b1, a2_b2;
##' Cell pair B_A includes marker pairs: b1_a1, b1_a2, b2_a1, b2_a2
##' 
##'
##' @param geneset a data.frame of cell markers, 
##' in which the 1st column is cell type, 
##' the 2nd column is the marker gene 
##' and the 3rd column is the name of the dataset (optional)
##' @param core a numeric value shows the number of cores 
##' to use for parallel execution. The default value is 1.
##' @param dataset specified at least 1 dataset occured in the 3rd column of geneset
##' @import foreach
##' @import doParallel
##' @return A dataframe of Cell pair and the corresponding marker pair
##' @export 
##' @examples
##' \dontrun{
##'   data(CellType_Galon2013_cancer)
##'   geneset <- CellType_Galon2013_cancer
##'   cell_pair <- TimiCellPair(geneset = geneset,core = 2)
##' }
##' @author Chenyang Skylar Li
TimiCellPair <- function(geneset = NULL,
                         dataset = NULL,
                         core = 1){
  
  # examine required parameters-------------------------------------------------
  if (is.null(geneset)){
    stop('The parameter "geneset" is required. ')
  }
  # Parallel 
  registerDoParallel(cores=core)
  # optional 3rd column dataset-------------------------------------------------
  if (ncol(geneset) == 3) {
    
    num <- table(geneset[3])
    if(length(num) > 1 & is.null(dataset)){
      stop('There are more than 1 datasets. Please choose one.')
    } else if (length(num) == 1 & is.null(dataset)){
      dataset <- names(num)
    }
    
    if (!is.null(dataset)) {
      se <- which(geneset[,3] %in% dataset)
      if (length(se) == 0){
        stop('Given dataset is different from 3rd column of geneset.')
      } else{
        set <- geneset[se, 1:2]
      }
      
    }
  }
  
  if (ncol(geneset) == 2){
    set <- geneset
  }
  
  
  
  # Cell Pair ------------------------------------------------------------------
  ann <- unique(set[, 1])
  if(length(ann) <=1){
    stop('Less than 2 cell types were found.')
  }
  ann_int <- outer(ann, ann, paste, sep="_")
  nn <- length(ann_int)
  tmp <- unlist(strsplit(ann_int, "_"))
  tmp1 <- tmp[(1:nn)*2-1]
  tmp2 <- tmp[(1:nn)*2-0]
  se <- which(tmp1 != tmp2)
  ann_int <- ann_int[se] 
  # Marker pair annoation 
  geneset_cell_pair <- foreach (i = 1:length(ann_int), .combine=rbind) %dopar% {
    

    ann_i <- unlist(strsplit(ann_int[i], "_"))
    
    se1 <- which(set[, 1]==ann_i[1])
    se2 <- which(set[, 1]==ann_i[2])
    
    xx1 <- set[se1, 2]
    xx2 <- set[se2, 2]
    
    xx <- outer(xx1, xx2, paste, sep="_")
    nn <- length(xx)
    tmp <- unlist(strsplit(xx, "_"))
    tmp1 <- tmp[(1:nn)*2-1]
    tmp2 <- tmp[(1:nn)*2-0]
    se <- which(tmp1 != tmp2)
    xx <- xx[se] 
    nn <- length(xx)
    data.frame("Cell.Pair" = rep(ann_int[i],nn),
               "Marker.Pair" = xx)
  
  }
  return(geneset_cell_pair)
}
