##' Permutation of marker genes to control FDR
##' 
##' It conservatively permutates marker gene sets for all cell types 
##' by shuffling their marker genes while preserving the sizes of their sets. 
##' This permutation breaks the underlying association 
##' between the genes and the cell types, 
##' which enables the examination of cell-type marker quality 
##' and the control of false positives in identifying cell interactions. 
##' It exports the FDR (Q value) of each interaction.
##'
##' @param resdata TimiGP enrichment result generated from TimiEnrich
##' @param geneset a data.frame of cell markers used in TimiCellPair,
##' in which the 1st column is cell type, 
##' the 2nd column is the marker gene 
##' and the 3rd column is the name of the dataset (optional)
##' @param gene a vector of gene symbol used in TimiEnrich
##' @param background background genes or pairs used in TimiEnrich
##' @param niter A numeric value used to determine the number of permutation. 
##' The default cutoff is 100.
##' @param core a numeric value shows the number of cores 
##' for parallel execution. The default value is 1.
##' @return A new column "Permutation.FDR" is appended to 
##' dataframe of cell-cell interactions generated from TimiEnrich
##' @export 
##' @examples
##' \dontrun{
##'   data("Bindea2013c_enrich")
##'   res <- Bindea2013c_enrich
##'   res <-  TimiPermFDR(resdata = res, geneset = geneset, gene = GP,
##'                        background = background, niter = 100, core = 20)
##'   head(score)
##' }
##' @author Chenyang Skylar Li

TimiPermFDR <-  function(resdata = NULL,
                         geneset = NULL,
                         gene = NULL,
                         background = NULL,
                         niter = 100,
                         core = 1){
                    
  # Examine required parameters-------------------------------------------------
  if (is.null(resdata)){
    stop('The parameter "resdata" is required.')
  }
  if (is.null(geneset)){
    stop('The parameter "geneset" is required. ')
  }
  if (is.null(gene)){
    stop('The parameter "gene" is required.')
  }
  if (is.null(background)){
    stop('The parameter "background" is required.')
  }
  
  
  # Parallel 
  registerDoParallel(cores=core)
  
  # Permutation ----------------------------------------------------------------
  message("--------------- Run ", niter, 
          " permutations. Please be patient ---------------" )
  permutation_res <- data.frame()
  for (seed in 1:niter) {
    set.seed(seed)
    
    # Shuffle marker gene set 
    perm.set <- geneset
    mygen <-  perm.set$Gene
    perm.set$Gene <- mygen[sample(1:length(mygen))]
    
    message("Permeutation#",seed, 
            ": Different / Total markers: ",
            sum(perm.set$Gene != geneset$Gene),"/", nrow(geneset) )
    
    # Run TimiGP using shuffled cell-type marker sets
    cell_pair <- TimiCellPair(geneset = perm.set,core = core)

    res <- TimiEnrich(gene = gene, background = background,
                      geneset = cell_pair, p.adj = "none",core=core)
    
    res$seed <- seed
    permutation_res <- rbind(permutation_res,res)
  }
  # Permutation FDR ------------------------------------------------------------
  message("--------------- Permutation done. Calculating FDR ---------------")
  resdata$Permutation.FDR <- NULL
  
  for (i in 1:nrow(resdata)) {
    k <- resdata$Rank[i]
    p_k <- resdata$P.Value[i]
    resdata$Permutation.FDR[i] <- min(
      sum(permutation_res$P.Value <= p_k)/k/niter,1) 
  }
  resdata <- resdata[c(1:10,13,11:12)]
  
  message('The permutation has been succesfully perfomed.',
  ' Please check the new column "Permutation.FDR"')
  
  return(resdata)
}


