#' Aggregate gene expression profile GEP
#' ------------------------------------
# 'CIBERSORTx methods: "For each cell type represented by at least three single 
#' cells, we selected 50% of all available single-cell GEPs using random 
#' sampling without replacement (fractional sample sizes were rounded up such 
#' that two cells were sampled if only three were available). We then aggregated 
#' the profiles by summation in non-log linear space and normalized each 
#' population-level GEP into TPM. This process was repeated in order to generate 
#' five aggregated transcriptome replicates per cell type assume availabe cell 
#' types are > 3".
#'
#' @param scgep A matrix of normalized gene expression values, rows are genes,
#' columns are cells
#' @param meta A data.frame with meta data information for each cell. A mandatory
#' column name 'Cell' matches with the column names of scgep
#' @param meta_column Name of the column that contains the cell type information
#' @param cell_types A vector indicating the cell types to use for generating
#' the signature
#' @param scale_factor Scaling factor used for normalization
#' @param seed A seed number
#'
#' @return A matrix of aggregated GEP profiles
#' @export
#'
#' @examples
#' aggregateGEP(
#'  livnat_mat, 
#'  livnat_meta, 
#'  cell_types=c("T-cell", "Malignant"),
#'  meta_column="Cell_type",
#'  scale_factor=10000,
#'  seed=1,
#'  num_gep=5)
aggregateGEP <- function(
  scgep,
  meta_data,
  cell_types,
  meta_column="Celltype (major-lineage)",
  scale_factor=10000,
  seed = 1,
  num_gep = 5,
  verbose=FALSE){
  
  # TODO:
  # add quality checks for input
  if(is.data.frame(scgep)) {scgep <- as.matrix(scgep)}
  if(verbose){
    print("Aggregate GEP: ")
    print(str(scgep))
    print(scgep[1:5, 1:5])
    print(paste("Cell type column: ", meta_column, sep=""))
    print(paste("Seed: ", seed, sep=""))
    print(head(meta_data))
  }
  
  # Set seed
  set.seed(seed)
  
  # Aggregated cell types
  agg_ct <- matrix(nrow=dim(scgep)[1], ncol=length(cell_types) * 5)
  matCol <- 1
  
  for (i in 1:length(cell_types)){
    
    if(verbose) {print(cell_types[i])}
    
    # Get cells
    cell_type_specific <- meta_data[meta_data[[meta_column]]==cell_types[i], ]
    
    for(j in 1:num_gep){
      
      # Randomly sample 50% 
      sampled_cells <- base::sample(
        cell_type_specific$Cell, 
        size=round(dim(cell_type_specific)[1] * 0.5, digits = 0),
        replace=FALSE
      )
      
      # Aggregate the profiles by summation in non-log linear space
      cell_gep <- Matrix::rowSums(scgep[, sampled_cells])
      
      # Write result to new matrix
      agg_ct[, matCol] <- cell_gep
      matCol <- matCol + 1
    }
  }
  
  # Add meaningful column and row names
  colnames(agg_ct) <- paste(sort(rep(cell_types, num_gep)), seq(1:num_gep), sep=".")
  row.names(agg_ct) <- row.names(scgep)
  
  # ..and normalized each population-level GEP to scale factor 
  agg_ct_Colsum <- t(scale_factor * (t(agg_ct) / Matrix::colSums(agg_ct)))
  
  return(agg_ct_Colsum)
}