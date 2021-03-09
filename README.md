# scSigR
scSigR is an R package for generation of signature matrices from single cell data.

<img src="https://github.com/ruppinlab/scSigR/blob/main/data/images/Figure1.png" alt="grouping" width="1121" height="570">

## Installing scSigR
Install the development version directly from github (requires devtools)

```r
require(devtools)
devtools::install_github("ruppinlab/scSigR")

# If you don't want to install the package, load the function using devtools
devtools::load_all()
```

## Usage

```r
library(scSigR)
?scSigR
```

We use the following input file as an example:
```r
livnat_mat <- readRDS("GBM_forSignature_CV.rds")
livnat_meta <- readRDS("GBM_meta_forSignature_CV.rds")

> livnat_mat[1:3, 1:3]
         MGH102-P1-A01 MGH102-P1-A02 MGH102-P1-A03
A1BG                 0             0             0
A1BG-AS1             0             0             0
A1CF                 0             0             0

> livnat_meta[1:3, 1:4]
           Cell Sample GBM_Type Cell_Type
1 MGH102-P1-A01 MGH102    Adult Malignant
2 MGH102-P1-A02 MGH102    Adult Malignant
3 MGH102-P1-A03 MGH102    Adult Malignant
```

Step 1: Aggregate Gene expression profiles (default 5, num_gep).
```r
geps <- scSigR::aggregateGEP(
  livnat_mat, 
  livnat_meta, 
  cell_types=names(table(livnat_meta$Cell_Type)),
  meta_column="Cell_Type",
  scale_factor=10000,
  seed=1,
  num_gep=5,
  verbose = TRUE)
  
> geps[1:5, 1:3]
         Macrophage.1 Macrophage.2 Macrophage.3
A1BG         0.018044     0.015958     0.028440
A1BG-AS1     0.014134     0.033334     0.045636
A1CF         0.000578     0.001956     0.000988
A2M         17.243004    21.166668    16.667274
A2M-AS1      0.007664     0.006702     0.009444
```

Step 2: Perform differential expression. One cell type again all others.
```r
degs <- scSigR::deg(
  gep=as.matrix(gep), 
  cell_types=names(table(livnat_meta$Cell_Type)),
  num_gep=5,
  test_alternative="greater",
  verbose=TRUE)
  
> str(degs)
List of 4
 $ Macrophage     :'data.frame':	21394 obs. of  3 variables:
  ..$ p_value: num [1:21394] 0.166657 0.269061 0.154511 0.000594 0.126444 ...
  ..$ log2fc : num [1:21394] -0.00547 0.00256 0.000599 4.125409 0.001353 ...
  ..$ BH     : num [1:21394] 0.71124 0.9841 0.66929 0.00798 0.63921 ...
 $ Malignant      :'data.frame':	21394 obs. of  3 variables:
  ..$ p_value: num [1:21394] 0.854381 0.759161 0.006627 0.146536 0.000772 ...
  ..$ log2fc : num [1:21394] -0.0288 -0.01232 0.00155 -2.3772 0.02977 ...
  ..$ BH     : num [1:21394] 0.9029 0.8149 0.0137 0.1879 0.002 ...
```

Step 3: Find matrix with lowest condition number.
```r
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
  
> head(kappas)
   gene_iteration condition_number
28            327         3.522208
30            329         3.522257
29            328         3.522266
31            330         3.522270
1             300         3.522290
2             301         3.522322
```

Get lowest condition number
```r
optimalG <- kappas[kappas$condition_number==min(kappas$condition_number), ]$gene_iteration
```
Step 4: Get signature matrix.
```r
signature <- scSigR::getSignatureMatrix(
  agg_ct=geps, 
  diff_expr = degs, 
  cell_types = names(table(livnat_meta$Cell_Type)),
  qvalue=0.01,
  log2fc=1,
  optimalG=optimalG,
  num_gep=5)
  
> signature[1:5, 1:4]
        Macrophage    Malignant Oligodendrocyte  T-cell
C1QB      69.74890 3.772358e-05       0.0046928 0.00000
CCL3      61.06000 3.427602e-03       0.0000000 0.00000
FCGR3A    42.27122 4.432520e-04       0.0000000 0.00000
HLA-DRA  100.50797 4.126496e-02       0.0000000 4.17466
FCER1G    40.07535 1.700346e-02       0.0000000 0.00000
```
