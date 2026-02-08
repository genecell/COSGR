## COSG in R

Accurate and fast cell marker gene identification with COSG


COSG is a cosine similarity-based method for more accurate and scalable marker gene identification.

* COSG is a general method for cell marker gene identification across different data modalities, e.g., scRNA-seq, scATAC-seq and spatially resolved transcriptome data.
* Marker genes or genomic regions identified by COSG are more indicative and with greater cell-type specificity.
* COSG is ultrafast for large-scale datasets, and is capable of identifying marker genes for one million cells in less than two minutes.

The method and benchmarking results are described in [Dai et al., (2022)](https://academic.oup.com/bib/advance-article-abstract/doi/10.1093/bib/bbab579/6511197?redirectedFrom=fulltext). The preprint is available in [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.15.448484v1).

Here is the R version for COSG, and the python version is hosted in https://github.com/genecell/COSG.

### Installation

```
# install.packages('remotes')
remotes::install_github(repo = 'genecell/COSGR')
```

### Seurat v5 Compatibility

COSGR v1.0.0+ is fully compatible with Seurat v3, v4, and v5.

For **Seurat v5** users, use the `layer` parameter:

```r
marker_cosg <- cosg(
  pbmc_small,
  groups='all',
  assay='RNA',
  layer='data',    # Use 'layer' for Seurat v5
  mu=1,
  n_genes_user=100)
```

For **Seurat v3/v4** users, the `slot` parameter continues to work:

```r
marker_cosg <- cosg(
  pbmc_small,
  groups='all',
  assay='RNA',
  slot='data',     # Use 'slot' for Seurat v3/v4
  mu=1,
  n_genes_user=100)
```

### Usage

Please check out the [vignette](https://github.com/genecell/COSGR/blob/master/vignettes/quick_start.Rmd) and the [PBMC10K tutorial](https://github.com/genecell/COSGR/blob/master/vignettes/pbmc10k_tutorial_cosg.Rmd) to get started.


Note I: we released our Python toolkit, [PIASO](https://piaso.org), in which some methods were built upon COSG.

Note II: we have also recently released [PIASOmarkerDB](https://piaso.org/piasomarkerdb), a cell type marker gene database for the single-cell and spatial transcriptomics community!


```
suppressMessages(library(Seurat))
data('pbmc_small',package='Seurat')
# Check cell groups:
table(Idents(pbmc_small))
#>
#>  0  1  2
#> 36 25 19
#######
# Run COSG (Seurat v5 - recommended):
marker_cosg <- cosg(
 pbmc_small,
 groups='all',
 assay='RNA',
 layer='data',
 mu=1,
 n_genes_user=100)
#######
# Check the marker genes:
 head(marker_cosg$names)
#>       0      1     2
#> 1   CD7 S100A8 MS4A1
#> 2  CCL5   TYMP CD79A
#> 3  GNLY S100A9 TCL1A
#> 4 LAMP1  FCGRT  NT5C
#> 5  GZMA IFITM3 CD79B
#> 6   LCK   LST1 FCER2
 head(marker_cosg$scores)
#>           0         1         2
#> 1 0.6391917 0.8954042 0.6922908
#> 2 0.6391267 0.8312083 0.5832425
#> 3 0.6328148 0.8120045 0.5757478
#> 4 0.6164937 0.7755955 0.5533107
#> 5 0.5846589 0.7413060 0.5163446
#> 6 0.5795238 0.7380483 0.5115180
####### Run COSG for selected groups, i.e., '0' and 2':
#######
marker_cosg <- cosg(
 pbmc_small,
 groups=c('0', '2'),
 assay='RNA',
 layer='data',
 mu=1,
 n_genes_user=100)
```

### Tip
1. If you would like to identify more specific marker genes, you could assign `mu` to larger values, such as `mu=10` or `mu=100`.
2. You could set the parameter `remove_lowly_expressed` to `TRUE` to not consider genes expressed very lowly in the target cell group, and you can use the parameter `expressed_pct` to adjust the threshold for the percentage. For example:
```
marker_region<-cosg(
    seo,
  groups='all',
  assay='peaks',
  layer='data',
  mu=100,
  n_genes_user=100,
  remove_lowly_expressed=TRUE,
  expressed_pct=0.1
)
```

### Citation

If COSG is useful for your research, please consider citing [Dai, M., Pei, X., Wang, X.-J., 2022. Accurate and fast cell marker gene identification with COSG. Brief. Bioinform. bbab579](https://academic.oup.com/bib/advance-article-abstract/doi/10.1093/bib/bbab579/6511197?redirectedFrom=fulltext).
