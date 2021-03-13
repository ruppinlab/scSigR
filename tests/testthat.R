# This file tests several functions.

#Libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(Metrics)

#Directory
data_dir <- "~/datasets/nbl-integrative-analysis/data/GBM/"

# Data GBM ---------------------------------------------------------------------
gbm_mat <- readRDS(paste(data_dir, "GBM_forSignature_CV.rds", sep=""))
gbm_meta <- readRDS(paste(data_dir, "GBM_meta_forSignature_CV.rds", sep=""))

# test aggregate GEP
geps <- scSigR::aggregateGEP(
  gbm_mat, 
  gbm_meta, 
  cell_types=names(table(gbm_meta$Cell_Type)),
  meta_column="Cell_Type",
  scale_factor=10000,
  seed=1,
  num_gep=5,
  verbose = TRUE)

# test differential gene expression
degs <- scSigR::deg(
  gep=as.matrix(geps), 
  cell_types=names(table(gbm_meta$Cell_Type)),
  num_gep=5,
  test_alternative="greater",
  verbose=TRUE)

# test get kappas function
kappas <- scSigR::getKappa(
  agg_ct=geps, 
  diff_expr=degs, 
  cell_types=names(table(gbm_meta$Cell_Type)),
  qvalue=0.05,
  log2fc=0.5, #fc=1.5
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
  cell_types = names(table(gbm_meta$Cell_Type)),
  qvalue=0.05,
  log2fc=0.5,
  optimalG=optimalG,
  num_gep=5)

# Export matrix ----------------------------------------------------------------
gbm_export <- gbm_mat
sorted_meta <- gbm_meta[order(gbm_meta$Cell_Type), ]
gbm_export <- gbm_export[sorted_meta$Cell]

write.table(
  livnat_export, 
  file=(paste(data_dir, "GBM_expression_matrix.txt", sep="\t")),
  quote=FALSE, col.names=as.character(sorted_meta$Cell_Type), sep="\t")

# Compare signatures -----------------------------------------------------------
sig_ciber <- read.csv(paste(data_dir, "CIBERSORTx_Job13_GBM_expression_matrix_inferred_phenoclasses.CIBERSORTx_Job13_GBM_expression_matrix_inferred_refsample.bm.K999.txt",sep=""),
                      sep="\t")
# Scale signature
signaturex100 <- signature * 100
signaturex100$NAME <- row.names(signaturex100)
signaturex100$NAME <- gsub("-", ".", signaturex100$NAME)

mer <- merge(signaturex100, sig_ciber, by="NAME", all.x=TRUE, all.y=TRUE)

bayesbio::jaccardSets(gsub("-", ".", row.names(signature)), sig_ciber$NAME)

# Compute cell fractions -------------------------------------------------------
source("~/src/scSigR/R/SVM-regression_CF.r")

bulk <- readRDS("~/datasets/nbl-integrative-analysis/data/GBM/GBM_bulk_tpm_CV.rds")
frac <- readRDS("~/datasets/nbl-integrative-analysis/data/GBM/GBM_cellfracs_CV.rds")
truth <- readRDS("~/datasets/nbl-integrative-analysis/data/GBM/GBM_groundtruth_tpm_CV.rds")

matchg <- extract_bulk(bulk, row.names(signature))
matSig <- as.matrix(signature * 100)
estfrac <- estimate_cell_fraction(bulk = matchg, mat_sig = matSig)

estfrac_melted <- reshape2::melt(estfrac)
colnames(estfrac_melted) <- c("patients", "cell_types", "predicted_fraction")

true_fraction_melted <- reshape2::melt(frac)
colnames(true_fraction_melted) <- c("patients", "cell_types", "estimated_fraction")

estfrac_melted$estimated_fraction <- true_fraction_melted$estimated_fraction

g <- ggplot(data=estfrac_melted, 
            aes(x=estimated_fraction, y=predicted_fraction, fill=cell_types))
g <- g + geom_smooth(method = "lm", se = FALSE, color="grey", lw=1)
g <- g + geom_point(pch=21, size=2, alpha=0.6)
g <- g + facet_wrap( ~ cell_types, ncol = 4, nrow = 3)
g <- g + theme_bw()
plot(g)

pdf("~/src/RNA-seq-immune-sig/figures/GBM_est_vs_pred_scSigR.pdf", 
    width=10, height=2.5)
plot(g)
dev.off()

data_correlation <- estfrac_melted %>% 
  group_by(cell_types) %>% 
  summarise(
    pearson = cor(predicted_fraction, estimated_fraction, method = "pearson"),
    spearman = cor(predicted_fraction, estimated_fraction, method = "spearman"),
    rmse = rmse(estimated_fraction, predicted_fraction))

# Now test with sig_ciber
ciber <- sig_ciber[2:5]
row.names(ciber) <- gsub("\\.", "-", sig_ciber$NAME)
ciber <- as.matrix(ciber)
matchg <- extract_bulk(bulk, row.names(ciber))
estfrac <- estimate_cell_fraction(bulk = matchg, mat_sig = ciber)

estfrac_melted <- reshape2::melt(estfrac)
colnames(estfrac_melted) <- c("patients", "cell_types", "predicted_fraction")

true_fraction_melted <- reshape2::melt(frac)
colnames(true_fraction_melted) <- c("patients", "cell_types", "estimated_fraction")

estfrac_melted$estimated_fraction <- true_fraction_melted$estimated_fraction

g <- ggplot(data=estfrac_melted, 
            aes(x=estimated_fraction, y=predicted_fraction, fill=cell_types))
g <- g + geom_smooth(method = "lm", se = FALSE, color="grey", lw=1)
g <- g + geom_point(pch=21, size=2, alpha=0.6)
g <- g + facet_wrap( ~ cell_types, ncol = 4, nrow = 3)
g <- g + theme_bw()
plot(g)

pdf("~/src/RNA-seq-immune-sig/figures/GBM_est_vs_pred_cibersortx.pdf", 
    width=10, height=2.5)
plot(g)
dev.off()

data_correlation <- estfrac_melted %>% 
  group_by(cell_types) %>% 
  summarise(
    pearson = cor(predicted_fraction, estimated_fraction, method = "pearson"),
    spearman = cor(predicted_fraction, estimated_fraction, method = "spearman"),
    rmse = rmse(estimated_fraction, predicted_fraction))

# Get signature from livnat ----------------------------------------------------
livnat_mat <- readRDS(paste(data_dir, "livnat_forSignature_CV.rds", sep=""))
livnat_meta <- readRDS(paste(data_dir, "livnat_meta_forSignature_CV.rds", sep=""))
livnat_bulk <- readRDS(paste(data_dir, "livnat_bulk_tpm_CV.rds", sep=""))
livnat_frac <- readRDS(paste(data_dir, "Livnat_cellfracs_CV.rds", sep=""))

colnames(livnat_meta) <- c("Cell", "celltype", "patient")

# test aggregate GEP
geps <- scSigR::aggregateGEP(
  livnat_mat, 
  livnat_meta, 
  cell_types=names(table(livnat_meta$celltype)),
  meta_column="celltype",
  scale_factor=1000000,
  seed=1,
  num_gep=5,
  verbose = TRUE)

# test differential gene expression
degs <- scSigR::deg(
  gep=as.matrix(geps), 
  cell_types=names(table(livnat_meta$celltype)),
  num_gep=5,
  stat_test="wilcox",
  test_alternative="greater",
  verbose=TRUE)

# test get kappas function
kappas <- scSigR::getKappa(
  agg_ct=geps, 
  diff_expr=degs, 
  cell_types=names(table(livnat_meta$celltype)),
  qvalue=0.05,
  log2fc=0.5, #fc=1.5
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
  cell_types = names(table(livnat_meta$celltype)),
  qvalue=0.05,
  log2fc=0.5,
  optimalG=optimalG,
  num_gep=5)

signature <- as.matrix(signature)
matchg <- livnat_bulk[match(intersect(row.names(signature), rownames(livnat_bulk)), rownames(livnat_bulk)),]
matchS <- signature[match(intersect(row.names(signature), rownames(livnat_bulk)), rownames(signature)),]
estfrac <- scSigR::estimate_cell_fraction(bulk = matchg, mat_sig = matchS)

livnat_frac <- livnat_frac[, colnames(estfrac)]

estfrac_melted <- reshape2::melt(estfrac)
colnames(estfrac_melted) <- c("patients", "cell_types", "predicted_fraction")

true_fraction_melted <- reshape2::melt(livnat_frac)
colnames(true_fraction_melted) <- c("patients", "cell_types", "estimated_fraction")

estfrac_melted$estimated_fraction <- true_fraction_melted$estimated_fraction

g <- ggplot(data=estfrac_melted, 
            aes(x=estimated_fraction, y=predicted_fraction, fill=cell_types))
g <- g + geom_smooth(method = "lm", se = FALSE, color="grey", lw=1)
g <- g + geom_point(pch=21, size=2, alpha=0.6)
g <- g + facet_wrap( ~ cell_types, ncol = 4, nrow = 3)
g <- g + theme_bw()
plot(g)

pdf("~/src/RNA-seq-immune-sig/figures/SKCM_est_vs_pred_scSigR.pdf", 
    width=8, height=5.5)
plot(g)
dev.off()

data_correlation <- estfrac_melted %>% 
  group_by(cell_types) %>% 
  summarise(
    pearson = cor(predicted_fraction, estimated_fraction, method = "pearson"),
    spearman = cor(predicted_fraction, estimated_fraction, method = "spearman"),
    rmse = rmse(estimated_fraction, predicted_fraction))

# Livnat export ----------------------------------------------------------------
livnat_export <- livnat_mat
sorted_meta <- livnat_meta[order(livnat_meta$celltype), ]
sorted_meta$celltype <- gsub("\\.", "", sorted_meta$celltype)
livnat_export <- livnat_export[as.character(sorted_meta$Cell)]

write.table(
  livnat_export, 
  file=(paste(data_dir, "livnat_expression_matrix.txt", sep="")),
  quote=FALSE, col.names=as.character(sorted_meta$celltype), sep="\t")

# Compare signature from CIBERSORTx --------------------------------------------
sig_ciber <- read.csv(paste(data_dir, "CIBERSORTx_Job15_livnat_expression_matrix_inferred_phenoclasses.CIBERSORTx_Job15_livnat_expression_matrix_inferred_refsample.bm.K999.txt",sep=""),
                      sep="\t")

ciber <- sig_ciber[2:11]
row.names(ciber) <- gsub("\\.", "-", sig_ciber$NAME)
ciber <- as.matrix(ciber)
matchg <- extract_bulk(livnat_bulk, row.names(ciber))
c_estfrac <- scSigR::estimate_cell_fraction(bulk = matchg, mat_sig = ciber)
c_estfrac_melted <- reshape2::melt(c_estfrac)
colnames(c_estfrac_melted) <- c("patients", "cell_types", "predicted_fraction")

colnames(livnat_frac) <- gsub("\\.", "", colnames(livnat_frac))
livnat_frac <- livnat_frac[, colnames(c_estfrac)]
true_fraction_melted <- reshape2::melt(livnat_frac)
colnames(true_fraction_melted) <- c("patients", "cell_types", "estimated_fraction")

c_estfrac_melted$estimated_fraction <- true_fraction_melted$estimated_fraction

g <- ggplot(data=c_estfrac_melted, 
            aes(x=estimated_fraction, y=predicted_fraction, fill=cell_types))
g <- g + geom_smooth(method = "lm", se = FALSE, color="grey", lw=1)
g <- g + geom_point(pch=21, size=2, alpha=0.6)
g <- g + facet_wrap( ~ cell_types, ncol = 4, nrow = 3)
g <- g + theme_bw()
plot(g)

pdf("~/src/RNA-seq-immune-sig/figures/SKCM_est_vs_pred_cibersortx.pdf", 
    width=8, height=5.5)
plot(g)
dev.off()

c_data_correlation <- c_estfrac_melted %>% 
  group_by(cell_types) %>% 
  summarise(
    pearson = cor(predicted_fraction, estimated_fraction, method = "pearson"),
    spearman = cor(predicted_fraction, estimated_fraction, method = "spearman"),
    rmse = rmse(estimated_fraction, predicted_fraction))

# Other comparision
mer <- merge(signature, sig_ciber, by="NAME", all.x=TRUE, all.y=TRUE)

bayesbio::jaccardSets(gsub("-", ".", row.names(signature)), sig_ciber$NAME)
         
