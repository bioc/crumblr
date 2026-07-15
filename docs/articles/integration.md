# Integration with dreamlet / SingleCellExperiment

## Load and process single cell data

Here we perform analysis of PBMCs from 8 individuals stimulated with
interferon-β [Kang, et al, 2018, Nature
Biotech](https://www.nature.com/articles/nbt.4042). We perform standard
processing with
[dreamlet](https://gabrielhoffman.github.io/dreamlet/index.html) to
compute pseudobulk before applying `crumblr`.

Here, single cell RNA-seq data is downloaded from
[ExperimentHub](https://bioconductor.org/packages/ExperimentHub/).

``` r

library(dreamlet)
library(muscat)
library(ExperimentHub)
library(scater)

# Download data, specifying EH2259 for the Kang, et al study
eh <- ExperimentHub()
sce <- eh[["EH2259"]]

sce$ind <- as.character(sce$ind)

# only keep singlet cells with sufficient reads
sce <- sce[rowSums(counts(sce) > 0) > 0, ]
sce <- sce[, colData(sce)$multiplets == "singlet"]

# compute QC metrics
qc <- perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]

# set variable indicating stimulated (stim) or control (ctrl)
sce$StimStatus <- sce$stim
```

### Aggregate to pseudobulk

Dreamlet creates the pseudobulk dataset:

``` r

# Since 'ind' is the individual and 'StimStatus' is the stimulus status,
# create unique identifier for each sample
sce$id <- paste0(sce$StimStatus, sce$ind)

# Create pseudobulk data by specifying cluster_id and sample_id for aggregating cells
pb <- aggregateToPseudoBulk(sce,
  assay = "counts",
  cluster_id = "cell",
  sample_id = "id",
  verbose = FALSE
)
```

### Process data

Here we evaluate whether the observed cell proportions change in
response to interferon-β.

``` r

library(crumblr)

# use dreamlet::cellCounts() to extract data
cellCounts(pb)[1:3, 1:3]
```

    ##          B cells CD14+ Monocytes CD4 T cells
    ## ctrl101      101             136         288
    ## ctrl1015     424             644         819
    ## ctrl1016     119             315         413

``` r

# Apply crumblr transformation
# cobj is an EList object compatable with limma workflow
# cobj$E stores transformed values
# cobj$weights stores precision weights
cobj <- crumblr(cellCounts(pb))
```

### Analysis

Now continue on with the downstream analysis

``` r

library(variancePartition)

fit <- dream(cobj, ~ StimStatus + ind, colData(pb))
fit <- eBayes(fit)

topTable(fit, coef = "StimStatusstim", number = Inf)
```

    ##                         logFC    AveExpr          t     P.Value  adj.P.Val         B
    ## CD8 T cells       -0.25085170  0.0857175 -4.0787416 0.002436375 0.01949100 -1.279815
    ## Dendritic cells    0.37386979 -2.1849234  3.1619195 0.010692544 0.02738587 -2.638507
    ## CD14+ Monocytes   -0.10525402  1.2698117 -3.1226341 0.011413912 0.02738587 -2.709377
    ## B cells           -0.10478652  0.5516882 -3.0134349 0.013692935 0.02738587 -2.940542
    ## CD4 T cells       -0.07840101  2.0201947 -2.2318104 0.050869691 0.08139151 -4.128069
    ## FCGR3A+ Monocytes  0.07425165 -0.2567492  1.6647681 0.128337022 0.17111603 -4.935304
    ## NK cells           0.10270672  0.3797777  1.5181860 0.161321761 0.18436773 -5.247806
    ## Megakaryocytes     0.01377768 -1.8655172  0.1555131 0.879651456 0.87965146 -6.198336

Given the results here, we see that CD8 T cells at others change
relative abundance following treatment with interferon-β.

## Session Info

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin23.6.0
    ## Running under: macOS Sonoma 14.7.1
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.33/lib/libopenblasp-r0.3.33.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.9.5             forcats_1.0.1               stringr_1.6.0              
    ##  [4] dplyr_1.2.1                 purrr_1.2.2                 readr_2.2.0                
    ##  [7] tidyr_1.3.2                 tibble_3.3.1                tidyverse_2.0.0            
    ## [10] glue_1.8.1                  crumblr_1.4.4               muscData_1.24.0            
    ## [13] scater_1.38.1               scuttle_1.20.0              ExperimentHub_3.0.0        
    ## [16] AnnotationHub_4.0.0         BiocFileCache_3.0.0         dbplyr_2.6.0               
    ## [19] muscat_1.24.0               dreamlet_1.9.3              SingleCellExperiment_1.32.0
    ## [22] SummarizedExperiment_1.40.0 Biobase_2.70.0              GenomicRanges_1.62.1       
    ## [25] Seqinfo_1.0.0               IRanges_2.44.0              S4Vectors_0.48.1           
    ## [28] BiocGenerics_0.56.0         generics_0.1.4              MatrixGenerics_1.22.0      
    ## [31] matrixStats_1.5.0           variancePartition_1.40.2    BiocParallel_1.44.0        
    ## [34] limma_3.66.0                ggplot2_4.0.3               BiocStyle_2.38.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] dichromat_2.0-0.1         GSEABase_1.72.0           progress_1.2.3           
    ##   [4] Biostrings_2.78.0         TH.data_1.1-5             vctrs_0.7.3              
    ##   [7] digest_0.6.39             png_0.1-9                 corpcor_1.6.10           
    ##  [10] shape_1.4.6.1             ggrepel_0.9.8             mixsqp_0.3-54            
    ##  [13] parallelly_1.48.0         MASS_7.3-66               fontLiberation_0.1.0     
    ##  [16] pkgdown_2.2.1             reshape2_1.4.5            SQUAREM_2026.1           
    ##  [19] foreach_1.5.2             withr_3.0.3               xfun_0.60                
    ##  [22] ggfun_0.2.1               survival_3.8-9            memoise_2.0.1            
    ##  [25] ggbeeswarm_0.7.3          emmeans_2.0.4             systemfonts_1.3.2        
    ##  [28] ragg_1.5.2                tidytree_0.4.8            zoo_1.8-15               
    ##  [31] GlobalOptions_0.1.4       gtools_3.9.5              KEGGgraph_1.70.0         
    ##  [34] prettyunits_1.2.0         KEGGREST_1.50.0           otel_0.2.0               
    ##  [37] httr_1.4.8                globals_0.19.1            ashr_2.2-63              
    ##  [40] babelgene_22.9            curl_7.1.0                ScaledMatrix_1.18.0      
    ##  [43] SparseArray_1.10.10       xtable_1.8-8              desc_1.4.3               
    ##  [46] doParallel_1.0.17         evaluate_1.0.5            S4Arrays_1.10.1          
    ##  [49] Rfast_2.1.5.2             hms_1.1.4                 bookdown_0.47            
    ##  [52] irlba_2.3.7               colorspace_2.1-3          filelock_1.0.3           
    ##  [55] magrittr_2.0.5            Rgraphviz_2.54.0          viridis_0.6.5            
    ##  [58] ggtree_4.0.5              lattice_0.22-9            future.apply_1.20.2      
    ##  [61] scattermore_1.2           XML_3.99-0.23             pillar_1.11.1            
    ##  [64] nlme_3.1-170              iterators_1.0.14          caTools_1.18.3           
    ##  [67] compiler_4.5.1            beachmat_2.26.0           stringi_1.8.7            
    ##  [70] rmeta_3.0                 minqa_1.2.8               plyr_1.8.9               
    ##  [73] msigdbr_26.1.0            crayon_1.5.3              abind_1.4-8              
    ##  [76] truncnorm_1.0-9           blme_1.0-7                metadat_1.6-0            
    ##  [79] gridGraphics_0.5-1        locfit_1.5-9.12           bit_4.6.0                
    ##  [82] mathjaxr_2.0-0            sandwich_3.1-2            codetools_0.2-20         
    ##  [85] multcomp_1.4-31           textshaping_1.0.5         BiocSingular_1.26.1      
    ##  [88] bslib_0.11.0              slam_0.1-56               GetoptLong_1.1.1         
    ##  [91] remaCor_0.0.20            splines_4.5.1             circlize_0.4.18          
    ##  [94] Rcpp_1.1.2                sparseMatrixStats_1.22.0  EnrichmentBrowser_2.40.0 
    ##  [97] knitr_1.51                blob_1.3.0                clue_0.3-68              
    ## [100] BiocVersion_3.22.0        lme4_2.0-1                fs_2.1.0                 
    ## [103] listenv_1.0.0             DelayedMatrixStats_1.32.0 Rdpack_2.6.6             
    ## [106] IHW_1.38.0                ggplotify_0.1.3           estimability_2.0.0       
    ## [109] Matrix_1.7-5              statmod_1.5.2             tzdb_0.5.0               
    ## [112] fANCOVA_0.6-1             pkgconfig_2.0.3           tools_4.5.1              
    ## [115] cachem_1.1.0              RhpcBLASctl_0.23-42       rbibutils_2.4.1          
    ## [118] RSQLite_3.53.3            viridisLite_0.4.3         DBI_1.3.0                
    ## [121] numDeriv_2016.8-1.1       zigg_0.0.2                fastmap_1.2.0            
    ## [124] rmarkdown_2.31            scales_1.4.0              grid_4.5.1               
    ## [127] broom_1.0.13              sass_0.4.10               patchwork_1.3.2          
    ## [130] coda_0.19-4.1             BiocManager_1.30.27       graph_1.88.1             
    ## [133] zenith_1.12.0             farver_2.1.2              reformulas_0.4.4         
    ## [136] aod_1.3.3                 mgcv_1.9-4                yaml_2.3.12              
    ## [139] cli_3.6.6                 lifecycle_1.0.5           mashr_0.2.79             
    ## [142] glmmTMB_1.1.14            mvtnorm_1.4-2             backports_1.5.1          
    ## [145] annotate_1.88.0           timechange_0.4.0          gtable_0.3.6             
    ## [148] rjson_0.2.23              metafor_5.0-1             parallel_4.5.1           
    ## [151] ape_5.8-1                 jsonlite_2.0.0            dirmult_0.1.3-5          
    ## [154] edgeR_4.8.2               bitops_1.0-9              bit64_4.8.2              
    ## [157] assertthat_0.2.1          yulab.utils_0.2.4         BiocNeighbors_2.4.0      
    ## [160] RcppParallel_5.1.11-2     jquerylib_0.1.4           pbkrtest_0.5.5           
    ## [163] lazyeval_0.2.3            htmltools_0.5.9           sctransform_0.4.3        
    ## [166] rappdirs_0.3.4            httr2_1.3.0               XVector_0.50.0           
    ## [169] gdtools_0.5.1             RCurl_1.98-1.19           treeio_1.34.0            
    ## [172] gridExtra_2.3.1           EnvStats_3.1.0            boot_1.3-32              
    ## [175] TMB_1.9.21                invgamma_1.2              R6_2.6.1                 
    ## [178] DESeq2_1.50.2             ggiraph_0.9.6             gplots_3.3.0             
    ## [181] fdrtool_1.2.18            cluster_2.1.8.2           aplot_0.3.1              
    ## [184] nloptr_2.2.1              DelayedArray_0.36.1       tidyselect_1.2.1         
    ## [187] vipor_0.4.7               fontBitstreamVera_0.1.1   AnnotationDbi_1.72.0     
    ## [190] future_1.70.0             rsvd_1.0.5                KernSmooth_2.23-26       
    ## [193] S7_0.2.2                  fontquiver_0.2.1          data.table_1.18.4        
    ## [196] htmlwidgets_1.6.4         ComplexHeatmap_2.26.1     RColorBrewer_1.1-3       
    ## [199] rlang_1.3.0               lmerTest_3.2-1            lpsymphony_1.38.0        
    ## [202] beeswarm_0.4.0
