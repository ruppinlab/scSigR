# This file tests several functions.

# Data
livnat_mat <- readRDS("~/datasets/nbl-integrative-analysis/data/GBM/GBM_forSignature_CV.rds")
livnat_meta <- readRDS("~/datasets/nbl-integrative-analysis/data/GBM/GBM_meta_forSignature_CV.rds")
sig <- readRDS("~/datasets/nbl-integrative-analysis/data/GBM/")


# test aggregate GEP
geps <- scSigR::aggregateGEP(
  livnat_mat, 
  livnat_meta, 
  cell_types=names(table(livnat_meta$Cell_Type)),
  meta_column="Cell_Type",
  scale_factor=10000,
  seed=1,
  num_gep=5,
  verbose = TRUE)

# test differential gene expression
degs <- scSigR::deg(
  gep=as.matrix(gep), 
  cell_types=names(table(livnat_meta$Cell_Type)),
  num_gep=5,
  test_alternative="greater",
  verbose=TRUE)

# test get kappas function
kappas <- scSigR::getKappa(
  agg_ct=geps, 
  diff_expr=degs, 
  cell_types=names(table(livnat_meta$Cell_Type)),
  qvalue=0.01,
  log2fc=1,
  num_gep=5,
  G_min=300,
  G_max=500,
  verbose=TRUE)

# Get lowest condition number
optimalG <- kappas[kappas$condition_number==min(kappas$condition_number), ]$gene_iteration

# test get signature matrix
signature <- scSigR::getSignatureMatrix(
  agg_ct=geps, 
  diff_expr = degs, 
  cell_types = names(table(livnat_meta$Cell_Type)),
  qvalue=0.01,
  log2fc=1,
  optimalG=optimalG,
  num_gep=5)
