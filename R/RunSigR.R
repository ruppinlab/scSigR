#' This function takes a seurat object and returns a signature matrix
#'
#' @param object A seurat object 'seurat_object'
#' @param cell_types A character vector with cell types present in celltype_column
#' @param celltype_column A character string which indicates the column name in 
#' seurat_object meta.data that contains the cell identities. Default: assign.ident
#' @param cellid_column A character string which indicates the column name in
#' seurat_object meta.data that contains the the cell ids. if NULL row.names of
#' meta.data will be used
#' @param scale_factor A number. Default: Signature matrix will be scaled to 10.000
#' @param seed A number.
#' @param num_gep A number. Number of GEP gene expression profiles to generate.
#' @param test_alternative ['wilcox', 'ttest']. Test to use for differential expression.
#' @param qvalue A number. Default: 0.01 Qvalue cutoff after differential expression.
#' @param log2fc A number. Default: 0.5. Log2fc cutoff afer differential expression.
#' @param G_min A number. Default 300. Minimum number of genes to start with and 
#' include in signature matrix
#' @param G_max A number. Default 500. Maximum number of genes to include in
#' signature matrix
#' @param verbose Boolean.
#'
#' @return
#' @export
#'
#' @examples
RunSigR <- function(
  object,
  cell_types,
  celltype_column="assign.ident",
  cellid_column=NULL,
  scale_factor=10000,
  seed=NULL,
  num_gep=5,
  test_alternative="greater",
  qvalue=0.01,
  log2fc=0.5,
  G_min=300,
  G_max=500,
  verbose=TRUE
){
  
  if(class(seurat_obj)[1]!="Seurat"){
    stop("This is not a Seurat object.")
  }
  
  if(is.null(cellid_column)){
    object@meta.data$Cell <- row.names(object@meta.data)
  }else{
    object@meta.data$Cell <- object@meta.data[[cellid_column]]
  }
  
  if(sum(cell_types %in% unique(object@meta.data[[celltype_column]])) != length(cell_types)){
    stop("Unknown cell types.")
  }
  
  extr_matrix <- as.matrix(object@assays$RNA@data)
  
  geps <- scSigR::aggregateGEP(
    extr_matrix, 
    object@meta.data, 
    cell_types=cell_types,
    celltype_column=celltype_column,
    scale_factor=10000,
    seed=seed,
    num_gep=5,
    verbose = TRUE)
  
  degs <- scSigR::deg(
    gep=as.matrix(geps), 
    cell_types=cell_types,
    num_gep=5,
    test_alternative="greater",
    verbose=TRUE)
  
  kappas <- scSigR::getKappa(
    agg_ct=geps, 
    diff_expr=degs, 
    cell_types=cell_types,
    qvalue=0.01,
    log2fc=0.5, #fc=1.5
    num_gep=5,
    G_min=300,
    G_max=500,
    verbose=TRUE)
  
  # Get lowest condition number
  optimalG <- kappas[kappas$condition_number==min(kappas$condition_number), ]$gene_iteration
  
  signature <- scSigR::getSignatureMatrix(
    agg_ct=geps, 
    diff_expr = degs, 
    cell_types = cell_types,
    qvalue=0.01,
    log2fc=0.5,
    optimalG=optimalG,
    num_gep=5)
  
  return(signature)
  
}
