#' This function calculates a condition number for each incremenal increase of
#' genes used for the signature matrix. Thee minimum number of genes to start
#' with is 300 and the max is 500. Genes with a q value < 0.01 and a log2fc < 0.5
#' are considered signfificant
#'
#' @param agg_ct A matrix of aggregated gene expression profiles
#' @param diff_expr List object of differential expression results
#' @param cell_types A vector with cell types of interest
#' @param qvalue qvalue threshold for differential expression
#' @param log2fc log2fc threshold for differential expression
#' @param G_min A number - the miminum number of genes to start with
#' @param G_max A number - Eq to G_min
#' @param num_gep A number - number of geps
#'
#' @return A dataframe with the condition number for each iteration.
#' @export
#'
#' @examples
#' kappa_df <- get_kappa(
#'  agg_ct=gep, 
#'  diff_expr = deg_res, 
#'  cell_types = cell_types,
#'  qvalue=0.01,
#'  log2fc=0.5,
#'  num_gep=5,
#'  G_min=300,
#'  G_max=500)
getKappa <- function(
  agg_ct,
  diff_expr, 
  cell_types,
  qvalue,
  log2fc,
  num_gep=5,
  G_min=300,
  G_max=500,
  verbose=TRUE){
  
  if(verbose){
    message("Get optimal condition number.")
    print(paste("qvalue: ", qvalue, sep=""))
    print(paste("log2fc: ", log2fc, sep=""))
    print(paste("G_min: ", G_min, sep=""))
    print(paste("G_max: ", G_max, sep=""))
  }
  
  # Condition number
  kappa_cond_num <- c()
  
  for (G in G_min:G_max){
    
    sub_diff <- list()
    
    for(i in 1:length(cell_types)){
      
      diff_ct <- diff_expr[[cell_types[i]]]
      
      # Genes with a q value <0.01 (false discovery rate) 
      # and a log2fc > 0.5 were considered significant.
      diff_sig <- diff_ct[diff_ct$BH < qvalue & diff_ct$log2fc > log2fc, ]
      
      # Sort by log2fc
      diff_sorted <- diff_sig[order(diff_sig$log2fc, decreasing=TRUE), ]
      
      # Select top G genes
      diff_top <- head(diff_sorted, G)
      sub_diff[[cell_types[i]]] <- diff_top
    }
    
    if(sum(unlist(lapply(sub_diff, function(x) {dim(x)[1]})) == 0) > 0){
      stop("No significant genes for one or more cell types. Try lowering qvalue.")
    }
    
    # Combined dataframes to single matrix
    # Replace missing genes with 1
    
    # Get genes
    genes <- c()
    
    for(i in 1:length(cell_types)){
      genes <- c(genes, rownames(sub_diff[[cell_types[i]]]))
    }
    
    genes_uniq <- unique(sort(genes))
    
    # Summarize matrix to single average expression
    av_expr <- data.frame(row.names = row.names(agg_ct))
    
    for(i in 1:length(cell_types)){
      col_ct <- paste(rep(cell_types[i], num_gep), seq(1:num_gep), sep=".")
      av_expr[[cell_types[i]]] <- rowMeans(agg_ct[, col_ct])
    }
    
    av_expr_sel <- av_expr[genes, ]
    
    mat_kappa <- as.matrix(av_expr_sel)
    kappa_cond_num <- c(kappa_cond_num, kappa(mat_kappa))
  }
  
  kappa_df <- data.frame(
    gene_iteration = G_min:G_max,
    condition_number = kappa_cond_num
  )
  
  # Reorder dataframe
  kappa_df <- kappa_df[order(kappa_df$condition_number), ]
  return(kappa_df)
}
