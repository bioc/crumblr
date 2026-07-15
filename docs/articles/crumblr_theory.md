# Normal approximation vs. empirical simulation

### Asymptotic normal approximation

Let the vector \\{\bf p}\\ be the true fractions across \\D\\
categories. Consider \\C\\ total counts sampled from a
Dirichlet-multinomial (DMN) distribution with overdispersion \\\tau\\,
where \\\tau=1\\ reduces to the multinomial distribution. The [centered
log ratio](https://rdrr.io/cran/compositions/man/clr.html) (CLR) of the
\\i^{th}\\ estimated fraction in \\{\bf \hat p}\\ is

\\\begin{equation} \tag{1} \text{clr}\_i({\bf \hat p}) = \log(\hat
p_i) - \frac{1}{D}\sum\_{j=1}^D \log(\hat p_j) \end{equation}\\

and we show that the sampling variance is

\\\begin{equation} \tag{2} \text{var}\[\text{clr}\_i({\bf \hat p})\] =
\frac{\tau}{C} \left\[ \frac{1}{\hat p_i} - \frac{2}{ D \hat p_i} +
\frac{1}{D^2}\sum\_{j=1}^D \frac{1}{\hat p_j} \right\]. \label{eqn2}
\end{equation}\\

### Simulations

The sampling variance is derived from asymototic theory, so we examine
its behavior for finite total counts. Here we evaluate the empirical
variance from \\1,000\\ draws from a Dirichlet-multinomial distribution
while varying \\D\\, \\\tau\\, \\C\\. A pseudocount of 0.5 is added to
the observed counts since the asymptotic theory is not defined for
counts of zero.

Here we plot the standard deviation after CLR transform from the
empirical DMN and the asymptotic normal approximation under a range of
conditions. Results are shown for instances with at least 2 counts.

#### D=2 categories

![](crumblr_theory_files/figure-html/plot1-1.png)

#### D=15 categories

![](crumblr_theory_files/figure-html/plot2-1.png)

#### Interpretation

The asymptotic standard deviation shows good agreement with the
empirical results even for small values of \\C\\, *when at least 2
counts are observed*. In practice, it is often reasonable to assume a
sufficient number of counts before a variable is included in an
analysis. Importantly, with less than 2 counts the asymptotic theory
gives a *larger* standard deviation than the emprical results (results
not shown). Therefore, this approach is conservative and should not
underestimate the true amount of variation. The asymptotic normal
approximation is most accurate for large total counts \\C\\, large
proportions \\p\\, and small overdispersion \\\tau\\.

##### Consideration of overdispersion

Based on Equation (2), the variance of the CLR-transformed proportions
is a *linear* function of \\\tau\\. Importantly, downstream analysis of
the CLR-transformed proportions with a precision-weighted linear (mixed)
model or a variance stabilizing transform depends only on the *relative*
variances. Since relative variances are invariant to the scale of
\\\tau\\, for these applications the value of \\\tau\\ can be set to 1
instead of being estimated from the data.

For other applications, `crumblr` can estimate \\\tau\\ from the data by
using `crumblr(counts, tau=NULL)`. This calls `dmn.mle()` to estimate
the parameters of the DMN distribution and is substantially faster than
alternatives.

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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] dirmult_0.1.3-5  lubridate_1.9.5  forcats_1.0.1    stringr_1.6.0    dplyr_1.2.1     
    ##  [6] purrr_1.2.2      readr_2.2.0      tidyr_1.3.2      tibble_3.3.1     tidyverse_2.0.0 
    ## [11] glue_1.8.1       crumblr_1.4.4    ggplot2_4.0.3    BiocStyle_2.38.0
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3          jsonlite_2.0.0              magrittr_2.0.5             
    ##   [4] farver_2.1.2                nloptr_2.2.1                rmarkdown_2.31             
    ##   [7] fs_2.1.0                    ragg_1.5.2                  vctrs_0.7.3                
    ##  [10] minqa_1.2.8                 ggtree_4.0.5                htmltools_0.5.9            
    ##  [13] S4Arrays_1.10.1             broom_1.0.13                SparseArray_1.10.10        
    ##  [16] gridGraphics_0.5-1          variancePartition_1.40.2    sass_0.4.10                
    ##  [19] KernSmooth_2.23-26          bslib_0.11.0                htmlwidgets_1.6.4          
    ##  [22] desc_1.4.3                  pbkrtest_0.5.5              plyr_1.8.9                 
    ##  [25] cachem_1.1.0                lifecycle_1.0.5             iterators_1.0.14           
    ##  [28] pkgconfig_2.0.3             Matrix_1.7-5                R6_2.6.1                   
    ##  [31] fastmap_1.2.0               rbibutils_2.4.1             MatrixGenerics_1.22.0      
    ##  [34] digest_0.6.39               numDeriv_2016.8-1.1         aplot_0.3.1                
    ##  [37] patchwork_1.3.2             S4Vectors_0.48.1            textshaping_1.0.5          
    ##  [40] GenomicRanges_1.62.1        labeling_0.4.3              timechange_0.4.0           
    ##  [43] abind_1.4-8                 compiler_4.5.1              fontquiver_0.2.1           
    ##  [46] aod_1.3.3                   withr_3.0.3                 S7_0.2.2                   
    ##  [49] backports_1.5.1             BiocParallel_1.44.0         viridis_0.6.5              
    ##  [52] gplots_3.3.0                MASS_7.3-66                 rappdirs_0.3.4             
    ##  [55] DelayedArray_0.36.1         corpcor_1.6.10              gtools_3.9.5               
    ##  [58] caTools_1.18.3              tools_4.5.1                 otel_0.2.0                 
    ##  [61] ape_5.8-1                   remaCor_0.0.20              nlme_3.1-170               
    ##  [64] grid_4.5.1                  reshape2_1.4.5              generics_0.1.4             
    ##  [67] gtable_0.3.6                tzdb_0.5.0                  hms_1.1.4                  
    ##  [70] XVector_0.50.0              BiocGenerics_0.56.0         pillar_1.11.1              
    ##  [73] yulab.utils_0.2.4           limma_3.66.0                splines_4.5.1              
    ##  [76] treeio_1.34.0               lattice_0.22-9              tidyselect_1.2.1           
    ##  [79] fontLiberation_0.1.0        SingleCellExperiment_1.32.0 knitr_1.51                 
    ##  [82] fontBitstreamVera_0.1.1     reformulas_0.4.4            gridExtra_2.3.1            
    ##  [85] bookdown_0.47               IRanges_2.44.0              Seqinfo_1.0.0              
    ##  [88] SummarizedExperiment_1.40.0 RhpcBLASctl_0.23-42         stats4_4.5.1               
    ##  [91] xfun_0.60                   Biobase_2.70.0              statmod_1.5.2              
    ##  [94] matrixStats_1.5.0           stringi_1.8.7               lazyeval_0.2.3             
    ##  [97] ggfun_0.2.1                 yaml_2.3.12                 boot_1.3-32                
    ## [100] evaluate_1.0.5              codetools_0.2-20            gdtools_0.5.1              
    ## [103] BiocManager_1.30.27         ggplotify_0.1.3             cli_3.6.6                  
    ## [106] RcppParallel_5.1.11-2       systemfonts_1.3.2           Rdpack_2.6.6               
    ## [109] jquerylib_0.1.4             dichromat_2.0-0.1           Rcpp_1.1.2                 
    ## [112] zigg_0.0.2                  EnvStats_3.1.0              Rfast_2.1.5.2              
    ## [115] pkgdown_2.2.1               bitops_1.0-9                lme4_2.0-1                 
    ## [118] viridisLite_0.4.3           mvtnorm_1.4-2               tidytree_0.4.8             
    ## [121] ggiraph_0.9.6               lmerTest_3.2-1              scales_1.4.0               
    ## [124] fANCOVA_0.6-1               rlang_1.3.0
