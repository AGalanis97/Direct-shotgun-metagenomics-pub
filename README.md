# Direct shotgun metagenomics captures species abundance of honey samples
This repository contains all files and code to reproduce the results in the publication [add link]

### Overview of folders and description
#### 1. Mocks
This folder contains: 1) The FASTQ files of all mocks used, including the kraken2 evaluation mocks; 2) The true abundance of each organism; 3) The detected abundance using each tool and kraken2 specifically for its evaluation. Only the mocks whose name begins with A,B, or C (example A64) were used for kraken2 evaluation.

#### 2. Figures
It contains a separate folder for each figure. Within the folders, subfolders contain the necessary data as well as the script to produce the figures. The figures will be saved within their folder when the script is run.

#### 3. Output data 
Contains the raw figures that were used in the publication. They are deposited here for direct comparison to any changes that may occur when adjusting the scripts. 

#### 3. Scripts with functions useful in this study
For archive purposes and distribution.

### This repository makes use of the here R package. When running the code from scripts, there is no need to change working directories - all paths are relative to the R project file located in this repository.

> sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252
[4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ALDEx2_1.18.0

loaded via a namespace (and not attached):
  [1] colorspace_1.4-1            ggtree_2.0.4                ggsignif_0.6.0              class_7.3-15                ggridges_0.5.2             
  [6] snakecase_0.11.0            htmlTable_1.13.3            XVector_0.26.0              GenomicRanges_1.38.0        base64enc_0.1-3            
 [11] rstudioapi_0.11             ggpubr_0.2.5                bit64_0.9-7                 lubridate_1.7.4             AnnotationDbi_1.48.0       
 [16] codetools_0.2-16            splines_3.6.3               doParallel_1.0.15           geneplotter_1.64.0          knitr_1.28                 
 [21] Formula_1.2-3               jsonlite_1.6.1              annotate_1.64.0             cluster_2.1.0               png_0.1-7                  
 [26] rgeos_0.5-2                 RPostgreSQL_0.6-2           BiocManager_1.30.10         compiler_3.6.3              httr_1.4.1                 
 [31] rvcheck_0.1.8               backports_1.1.5             ggcorrplot_0.1.3            assertthat_0.2.1            Matrix_1.2-18              
 [36] lazyeval_0.2.2              fasterize_1.0.2             acepack_1.4.1               htmltools_0.4.0             tools_3.6.3                
 [41] gtable_0.3.0                glue_1.3.2                  GenomeInfoDbData_1.2.2      dplyr_0.8.5                 Rcpp_1.0.4                 
 [46] Biobase_2.46.0              raster_3.0-12               vctrs_0.2.4                 ape_5.3                     nlme_3.1-144               
 [51] iterators_1.0.12            xfun_0.12                   stringr_1.4.0               lifecycle_0.2.0             XML_3.99-0.3               
 [56] zlibbioc_1.32.0             scales_1.1.0                parallel_3.6.3              SummarizedExperiment_1.16.1 RColorBrewer_1.1-2         
 [61] memoise_1.1.0               gridExtra_2.3               ggplot2_3.3.0               rpart_4.1-15                reshape_0.8.8              
 [66] latticeExtra_0.6-29         stringi_1.4.6               RSQLite_2.2.0               genefilter_1.68.0           S4Vectors_0.24.3           
 [71] foreach_1.4.8               tidytree_0.3.2              e1071_1.7-3                 checkmate_2.0.0             BiocGenerics_0.32.0        
 [76] BiocParallel_1.20.1         GenomeInfoDb_1.22.1         rlang_0.4.5                 pkgconfig_2.0.3             matrixStats_0.56.0         
 [81] bitops_1.0-6                lattice_0.20-38             purrr_0.3.3                 sf_0.8-1                    treeio_1.10.0              
 [86] htmlwidgets_1.5.1           cowplot_1.0.0               bit_1.1-15.2                tidyselect_1.0.0            GGally_1.5.0               
 [91] BIEN_1.2.4                  plyr_1.8.6                  magrittr_1.5                DESeq2_1.26.0               R6_2.4.1                   
 [96] IRanges_2.20.2              Hmisc_4.4-0                 DelayedArray_0.12.2         DBI_1.1.0                   pillar_1.4.3               
[101] foreign_0.8-75              units_0.6-6                 survival_3.1-8              RCurl_1.98-1.1              sp_1.4-1                   
[106] nnet_7.3-12                 tibble_2.1.3                janitor_2.0.1               crayon_1.3.4                KernSmooth_2.23-16         
[111] plotly_4.9.2.1              jpeg_0.1-8.1                locfit_1.5-9.1              grid_3.6.3                  data.table_1.12.8          
[116] blob_1.2.1                  digest_0.6.25               classInt_0.4-3              xtable_1.8-4                tidyr_1.0.2                
[121] stats4_3.6.3                munsell_0.5.0               viridisLite_0.3.0          
