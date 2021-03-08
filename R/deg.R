#' Perform differential gene expression. One cell type against all others.
#'
#' @param gep A matrix of aggregated gene expression profiles
#' @param cell_types A vector with cell types to consider
#' @param num_gep 
#' @param test_alternative What alternative to use c("two.sided", "less", "greater")
#'
#' @return A list object.
#' @export
#'
#' @examples
#' deg(
#' gep=as.matrix(gep), 
#' cell_types=c("T-cell", "Malignant"),
#' num_gep=5,
#' test_alternative="greater")
deg <- function(
  gep,
  cell_types,
  num_gep=5,
  test_alternative="greater",
  verbose=FALSE){
  
  agg_ct_nz <- gep[rowSums(gep) != 0, ]
  
  diff_expr <- list()
  
  for (i in 1:length(cell_types)){
    
    # Print cell types
    if(verbose){print(cell_types[i])}
    
    col_ct <- paste(rep(cell_types[i], num_gep), seq(1:num_gep), sep=".")
    # Var b has the cell types to compare to (eg. all T cells)
    b <- colnames(gep) %in% col_ct
    
    # Run Wilcoxon tests (eg. T cells against all other cell types)
    suppressWarnings({
      results <- t(apply(agg_ct_nz, 1, function(x) wilcox_class(
        x, b, test_alternative=test_alternative)))
    })
    
    colnames(results) <- c('p_value','log2fc')
    results <- as.data.frame(results)
    
    # Perform FDR correction using Benjamini Hochberg
    results$BH <- p.adjust(results$p_value, "BH")
    diff_expr[[cell_types[i]]] <- results
  }
  
  return(diff_expr)
}


#' Perform wilcox test.
#'
#' @param v A gep matrix 
#' @param b A vector with column names used for test
#' @param test_alternative What alternative to use c("two.sided", "less", "greater")
#'
#' @return
#' @export
#'
#' @examples
wilcox_class <- function(
  v, 
  b, 
  test_alternative="two.sided"){
  
  # Get pvalue
  p <- wilcox.test(v[b],v[!b], alternative = test_alternative)$p.value
  
  # Get fold change
  # add small pseudocount
  vpseudo <- v + 1
  log2fc <- log2(mean(na.omit(vpseudo[b])/mean(na.omit(vpseudo[!b]))))
  
  return(c(p, log2fc))
}