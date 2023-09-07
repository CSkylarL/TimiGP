##' Generate marker pair background
##' 
##' Given a gene pair A_B, it will generate 
##' the background gene pairs A_B and B_A
##'
##' @param marker.pair a vector of marker pairs connected by "_"
##' @return A vector of background gene pairs
##' @export 
##' @examples
##' \dontrun{
##'   data(Bindea2013c_COX_MP_SKCM06)
##'   cox_res <- Bindea2013c_COX_MP_SKCM06
##'   background <- TimiBG(marker.pair = row.names(cox_res))
##' }
##' @author Chenyang Skylar Li
TimiBG <- function(marker.pair = NULL){
  
  # Examine required parameters-------------------------------------------------
  if (is.null(marker.pair)){
    stop('The parameter "marker.pair" is required.')
  }
  
  xx <- marker.pair
  tmp <- unlist(strsplit(xx, "_"))
  nn <- length(xx)
  tmp1 <- tmp[(1:nn)*2-1]
  tmp2 <- tmp[(1:nn)*2-0]
  rxx <- paste(tmp2, tmp1, sep="_")
  background <- unique(c(xx,rxx)) # Total gene pairs
  message(length(background)," Background Pairs were generated from given ",length(marker.pair)," pairs") 
  return(background)
}

##' Enrichment Analysis
##' 
##' It statistically determines the inter-interactions
##' by performing over-representation (enrichment analysis) 
##'
##' @param gene a vector of gene symbol
##' @param geneset a data frame, column 1 shows cell (pair), 
##' column 2 shows gene (pair)
##' @param background background genes or pairs
##' @param p.adj p.adjust.methods. 
##' One of "holm", "hochberg", "hommel", "bonferroni", 
##' "BH", "BY", "fdr", "none". The default is "BH".
##' @param core a numeric value shows the number of cores for parallel execution. 
##' The default value is 1.
##' @param pair a logical value. If TRUE (default), 
##' perform enrichment analysis in cell pair. 
##' If using the function for general cell type enrichment analysis, 
##' please choose FALSE.
##' @import foreach
##' @import doParallel
##' @return A dataframe of enrichment results including cell-cell interactions
##' @export 
##' @examples
##' \dontrun{
##'   # Generate cell interaction Annotation: TimiCellPair
##'   data(CellType_Bindea2013_cancer)
##'   geneset <- CellType_Bindea2013_cancer
##'   cell_pair <- TimiCellPair(geneset = geneset,core = 20)
##'   # Select marker pairs A_B=1 associated with good prognosis
##'   data(Bindea2013c_COX_MP_SKCM06)
##'   cox_res <- Bindea2013c_COX_MP_SKCM06
##'   GP <- rownames(cox_res)[which(cox_res$QV<0.05)]
##'   # generate background: TimiBG
##'   background <- TimiBG(marker.pair = row.names(cox_res))
##'   # Enrichment Analysis: TimiEnrich
##'   res <- TimiEnrich(gene = GP, background = background, 
##'                     geneset = cell_pair, p.adj = "BH",core=2)
##'   
##' }
##' @author Chenyang Skylar Li

TimiEnrich <- function(gene = NULL,
                       background = NULL,
                       geneset = NULL,
                       p.adj =  "BH",
                       core = 1,
                       pair=TRUE){
  # examine required parameters------------------------------------------------
  if (is.null(gene)){
    stop('The parameter "gene" is required.')
  }
  if (is.null(background)){
    stop('The parameter "background" is required.')
  }
  if (is.null(geneset)){
    stop('The parameter "geneset" is required.')
  }
  if (sum(p.adj %in% 
          c("holm", "hochberg", "hommel", 
            "bonferroni","BH", "BY","fdr", "none")) !=1 ){
    stop('The parameter "p.adj" is required.' ,
         'Please choose one of c("holm", "hochberg",' ,
         ' "hommel", "bonferroni", "BH", "BY", "fdr", "none")')
  }
  # Parallel 
  registerDoParallel(cores=core)

  ## prepare parameter for hypergeometric test---------------------------------
  ##               Annotated Cell Interaction X->Y | Not Annotated | Row Sum
  ## Prognostic Pair |     k                       |     n-k          |    n
  ## Not Prognostic  |    M-k                      |    N+k-n-m       |
  ## column Sum      |     M                      |       N-M        |    N
  gene <- as.character(unique(gene))
  bg <- as.character(unique(background))
  ann <- unique(geneset[, 1])
  gene_ann <- as.character(unique(geneset[, 2]))
  
  ## N is the total number of genes in the background distribution
  ## The background distribution by default is all the genes that have annotation.
  bg_ann<- intersect(bg,gene_ann)
  N <- length(bg_ann)
  # output
  
  res <- foreach (i = 1:length(ann), .combine=rbind) %dopar% {
    se <- which(geneset[, 1]==ann[i])
    ## M is the number of genes within that distribution that are annotated 
    ## (either directly or indirectly) to the gene set of interest
    M_gene <- unique(geneset[se, 2])
    se <- which(M_gene %in% bg_ann)
    M_gene <- M_gene[se]
    M <- length(M_gene)
    ## n is the size of the list of genes of interest
    ## those genes that have no annotation should drop.
    n <- sum(gene %in% bg_ann)
    ## k is the number of genes within intereseting gene list which 
    ## are annotated to the gene set.
    k_gene <- gene[which(gene %in% M_gene)]
    k <- length(k_gene)
    data.frame("Cell.Interaction" = ann[i],
               "No.Total.IMGP" = M, 
               "No.Shared.IMGP" = k, 
               "Enrichment.Ratio" = round(k*N/M/n,3), 
               "P.Value" = phyper(k-1, M, N-M, n, lower.tail=FALSE), 
               "Shared.IMGP" =  paste(k_gene,collapse="/"),
               "Total.IMGP" = paste(M_gene,collapse="/"))             
  }
  ## P adjust------------------------------------------------------------------
  res$Adjust.P.Value <- p.adjust(res$P.Value, method=p.adj)
  res$Rank <- rank(res$P.Value,ties.method = "min")
  res <- res[order(res$P.Value,res$Adjust.P.Value), ]
  res$Index <- 1:nrow(res)
  if (pair == TRUE){
    xx <- unlist(strsplit(res$Cell.Interaction,split = "_"))
    res$Favorable.Cell.Type <- xx[seq(1,length(xx),2)]
    res$Unfavorable.Cell.Type <- xx[seq(2,length(xx),2)]
    res <- res[, c(10,9,1,11:12,2:5,8,6:7)]
  } else {
    res <- res[, c(10,9,1:5,8,6:7)]
    
  }
  
  
  
  
  return (res)
}


