#' Aggregate gene expression profile GEP
#' ------------------------------------
# 'For each cell type represented by at least three single cells, we selected 
# '50% of all available single-cell GEPs using random sampling without 
#' replacement (fractional sample sizes were rounded up such that two cells 
#' were sampled if only three were available). We then aggregated the profiles 
#' by summation in non-log linear space and normalized each population-level 
#' GEP into TPM. This process was repeated in order to generate five aggregated 
#' transcriptome replicates per cell type assume availabe cell types are > 3
#'
#' @param matrix 
#' @param meta 
#' @param cell_types 
#' @param scale_factor 
#' @param seed 
#'
#' @return A matrix of aggregated GEP profiles
#' @export
#'
#' @examples
aggregateGEP <- function(
  matrix,
  meta,
  meta_column="Celltype (major-lineage)",
  cell_types,
  scale_factor=10000,
  seed = 1,
  num_gep = 5,
  verbose=FALSE){
  
  # TODO:
  # add quality checks for input
  
  # set seed
  set.seed(seed)
  
  # Aggregated cell types
  agg_ct <- matrix(nrow=dim(matrix)[1], ncol=length(cell_types) * 5)
  matCol <- 1
  
  for (i in 1:length(cell_types)){
    
    if(verbose) {print(cell_types[i])}
    
    # Get cell type
    cell_type_specific <- meta[meta[[meta_column]]==cell_types[i], ]
    
    for(j in 1:num_gep){
      
      # Randomly sample 50% 
      sampled_cells <- base::sample(
        cell_type_specific$Cell, 
        size=round(dim(cell_type_specific)[1] * 0.5, digits = 0),
        replace=FALSE
      )
      
      # Aggregate the profiles by summation in non-log linear space
      cell_gep <- Matrix::rowSums(matrix[, sampled_cells])
      
      # Write result to new matrix
      agg_ct[, matCol] <- cell_gep
      matCol <- matCol + 1
    }
  }
  
  # Add meaningful column and row names
  colnames(agg_ct) <- paste(sort(rep(cell_types, num_gep)), seq(1:num_gep), sep=".")
  row.names(agg_ct) <- row.names(matrix)
  
  # ..and normalized each population-level GEP to scale factor 
  agg_ct_Colsum <- t(scale_factor * (t(agg_ct) / Matrix::colSums(agg_ct)))
  
  return(agg_ct_Colsum)
}