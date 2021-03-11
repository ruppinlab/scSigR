#' Get the final signature matrix
#'
#' @param agg_ct A matrix of aggregated gene expression profiles
#' @param diff_expr A list of cell types and results for differential gene expr
#' @param cell_types A vector of cell types of interest
#' @param qvalue qvalue threshold for differential expression
#' @param log2fc log2fc threshold for differential expression
#' @param optimalG The optimal number of genes for the signature matrix
#' @param num_gep A number - number of geps
#'
#' @return A dataframe with cell types and average gene expression values
#' @export
#'
#' @examples
#' signature <- scSigR::getSignatureMatrix(
#' agg_ct=geps, 
#' diff_expr = degs, 
#' cell_types = names(table(livnat_meta$Cell_Type)),
#' qvalue=0.01,
#' log2fc=0.5,
#' optimalG=optimalG,
#' num_gep=5)
getSignatureMatrix <- function(
  agg_ct,
  diff_expr, 
  cell_types,
  qvalue=0.01,
  log2fc=0.5,
  optimalG,
  num_gep=5,
  verbose=FALSE){
  
  sub_diff <- list()
  
  for(i in 1:length(cell_types)){
    
    if(verbose) { print(cell_types[i]) }
    
    diff_ct <- diff_expr[[cell_types[i]]]
    
    # Genes with a q value <0.01 (false discovery rate) 
    # and logfc < 0.5 were considered significant
    diff_sig <- diff_ct[diff_ct$BH < qvalue & diff_ct$log2fc > log2fc, ]
    
    # Sort by logfc
    diff_sorted <- diff_sig[order(diff_sig$log2fc, decreasing=TRUE), ]
    
    # Select optimal number of genes
    diff_top <- head(diff_sorted, optimalG)
    sub_diff[[cell_types[i]]] <- diff_top
  }
  
  # Combined dataframes to single matrix
  
  # Get genes
  genes <- c()
  
  for(i in 1:length(cell_types)){
    
    if(verbose) { print(cell_types[i]) }
    
    genes <- c(genes, rownames(sub_diff[[cell_types[i]]]))
  }
  
  if(verbose){print("Number of genes: ", length(genes))}

  # Summarize matrix to single average expression
  av_expr <- data.frame(row.names = row.names(agg_ct))
  
  for(i in 1:length(cell_types)){
    col_ct <- paste(rep(cell_types[i], num_gep), seq(1:num_gep), sep=".")
    av_expr[[cell_types[i]]] <- rowMeans(agg_ct[, col_ct])
  }
  
  av_expr_sel <- av_expr[genes, ]
  return(av_expr_sel)
  
}
  