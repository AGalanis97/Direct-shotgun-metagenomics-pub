# Direct shotgun metagenomics captures species abundance of honey samples
This repository contains all files and code to reproduce the results in the publication [add link]

### Overview of folders and description
#### 1. Mocks
This folder contains: 1) The FASTQ files of all mocks used, including the kraken2 evaluation mocks; 2) The true abundance of each organism; 3) The detected abundance using each tool and kraken2 specifically for its evaluation. Only the mocks whose name begins with A,B, or C (example A64) were used for kraken2 evaluation.

#### 2. Figures
It contains a separate folder for each figure. Within the folders, subfolders contain the necessary data as well as the script to produce the figures. The figures will be saved within their folder when the script is run.

#### 3. Output data 
Contains the raw figures that were used in the publication. They are deposited here for direct comparison to any changes that may occur when adjusting the scripts. 

#### 3. Useful scripts
It contains bash scripts used to generate mocks and evaluate kraken2. Mostly for archive purposes and distribution.

This repository makes use of the here R package. When running the code from scripts, there is no need to change working directories - all paths are relative to the R project file located in this repository.

<summary> <span title='Click to Expand'> current session info </span> </summary>

```r

- Session info -------------------------------------------------------------------------------------------------------------------------------------------
 setting  value                       
 version  R version 3.6.3 (2020-02-29)
 os       Windows 10 x64              
 system   x86_64, mingw32             
 ui       RStudio                     
 language (EN)                        
 collate  English_United Kingdom.1252 
 ctype    English_United Kingdom.1252 
 tz       Europe/Istanbul             
 date     2020-09-30                  

- Packages -----------------------------------------------------------------------------------------------------------------------------------------------
 package              * version  date       lib source        
 acepack                1.4.1    2016-10-29 [1] CRAN (R 3.6.3)
 ALDEx2               * 1.18.0   2019-10-29 [1] Bioconductor  
 annotate               1.64.0   2019-10-30 [1] Bioconductor  
 AnnotationDbi          1.48.0   2019-10-29 [1] Bioconductor  
 ape                    5.3      2019-03-17 [1] CRAN (R 3.6.3)
 assertthat             0.2.1    2019-03-21 [1] CRAN (R 3.6.3)
 backports              1.1.5    2019-10-02 [1] CRAN (R 3.6.1)
 base64enc              0.1-3    2015-07-28 [1] CRAN (R 3.6.0)
 BIEN                   1.2.4    2020-02-27 [1] CRAN (R 3.6.3)
 Biobase                2.46.0   2019-10-29 [1] Bioconductor  
 BiocGenerics           0.32.0   2019-10-29 [1] Bioconductor  
 BiocManager            1.30.10  2019-11-16 [1] CRAN (R 3.6.3)
 BiocParallel           1.20.1   2019-12-21 [1] Bioconductor  
 bit                    1.1-15.2 2020-02-10 [1] CRAN (R 3.6.2)
 bit64                  0.9-7    2017-05-08 [1] CRAN (R 3.6.2)
 bitops                 1.0-6    2013-08-17 [1] CRAN (R 3.6.0)
 blob                   1.2.1    2020-01-20 [1] CRAN (R 3.6.3)
 checkmate              2.0.0    2020-02-06 [1] CRAN (R 3.6.3)
 class                  7.3-15   2019-01-01 [2] CRAN (R 3.6.3)
 classInt               0.4-3    2020-04-07 [1] CRAN (R 3.6.3)
 cli                    2.0.2    2020-02-28 [1] CRAN (R 3.6.3)
 clipr                  0.7.0    2019-07-23 [1] CRAN (R 3.6.3)
 cluster                2.1.0    2019-06-19 [2] CRAN (R 3.6.3)
 codetools              0.2-16   2018-12-24 [2] CRAN (R 3.6.3)
 colorspace             1.4-1    2019-03-18 [1] CRAN (R 3.6.1)
 cowplot                1.0.0    2019-07-11 [1] CRAN (R 3.6.3)
 crayon                 1.3.4    2017-09-16 [1] CRAN (R 3.6.3)
 data.table             1.12.8   2019-12-09 [1] CRAN (R 3.6.3)
 DBI                    1.1.0    2019-12-15 [1] CRAN (R 3.6.3)
 DelayedArray           0.12.2   2020-01-06 [1] Bioconductor  
 desc                   1.2.0    2018-05-01 [1] CRAN (R 3.6.3)
 DESeq2                 1.26.0   2019-10-29 [1] Bioconductor  
 details              * 0.2.1    2020-01-12 [1] CRAN (R 3.6.3)
 digest                 0.6.25   2020-02-23 [1] CRAN (R 3.6.3)
 doParallel             1.0.15   2019-08-02 [1] CRAN (R 3.6.3)
 dplyr                  0.8.5    2020-03-07 [1] CRAN (R 3.6.3)
 e1071                  1.7-3    2019-11-26 [1] CRAN (R 3.6.3)
 fansi                  0.4.1    2020-01-08 [1] CRAN (R 3.6.3)
 fasterize              1.0.2    2020-03-25 [1] CRAN (R 3.6.3)
 foreach                1.4.8    2020-02-09 [1] CRAN (R 3.6.3)
 foreign                0.8-75   2020-01-20 [2] CRAN (R 3.6.3)
 Formula                1.2-3    2018-05-03 [1] CRAN (R 3.6.0)
 genefilter             1.68.0   2019-10-29 [1] Bioconductor  
 geneplotter            1.64.0   2019-10-29 [1] Bioconductor  
 GenomeInfoDb           1.22.1   2020-03-27 [1] Bioconductor  
 GenomeInfoDbData       1.2.2    2020-03-21 [1] Bioconductor  
 GenomicRanges          1.38.0   2019-10-29 [1] Bioconductor  
 GGally                 1.5.0    2020-03-25 [1] CRAN (R 3.6.3)
 ggcorrplot             0.1.3    2019-05-19 [1] CRAN (R 3.6.3)
 ggplot2                3.3.0    2020-03-05 [1] CRAN (R 3.6.3)
 ggpubr                 0.2.5    2020-02-13 [1] CRAN (R 3.6.3)
 ggridges               0.5.2    2020-01-12 [1] CRAN (R 3.6.3)
 ggsignif               0.6.0    2019-08-08 [1] CRAN (R 3.6.3)
 ggtree                 2.0.4    2020-04-13 [1] Bioconductor  
 glue                   1.3.2    2020-03-12 [1] CRAN (R 3.6.3)
 gridExtra              2.3      2017-09-09 [1] CRAN (R 3.6.3)
 gtable                 0.3.0    2019-03-25 [1] CRAN (R 3.6.3)
 Hmisc                  4.4-0    2020-03-23 [1] CRAN (R 3.6.3)
 htmlTable              1.13.3   2019-12-04 [1] CRAN (R 3.6.3)
 htmltools              0.4.0    2019-10-04 [1] CRAN (R 3.6.3)
 htmlwidgets            1.5.1    2019-10-08 [1] CRAN (R 3.6.3)
 httr                   1.4.1    2019-08-05 [1] CRAN (R 3.6.3)
 IRanges                2.20.2   2020-01-13 [1] Bioconductor  
 iterators              1.0.12   2019-07-26 [1] CRAN (R 3.6.3)
 janitor                2.0.1    2020-04-12 [1] CRAN (R 3.6.3)
 jpeg                   0.1-8.1  2019-10-24 [1] CRAN (R 3.6.1)
 jsonlite               1.6.1    2020-02-02 [1] CRAN (R 3.6.3)
 KernSmooth             2.23-16  2019-10-15 [2] CRAN (R 3.6.3)
 knitr                  1.28     2020-02-06 [1] CRAN (R 3.6.3)
 lattice                0.20-38  2018-11-04 [2] CRAN (R 3.6.3)
 latticeExtra           0.6-29   2019-12-19 [1] CRAN (R 3.6.3)
 lazyeval               0.2.2    2019-03-15 [1] CRAN (R 3.6.3)
 lifecycle              0.2.0    2020-03-06 [1] CRAN (R 3.6.3)
 locfit                 1.5-9.1  2013-04-20 [1] CRAN (R 3.6.3)
 lubridate              1.7.4    2018-04-11 [1] CRAN (R 3.6.3)
 magrittr               1.5      2014-11-22 [1] CRAN (R 3.6.3)
 Matrix                 1.2-18   2019-11-27 [2] CRAN (R 3.6.3)
 matrixStats            0.56.0   2020-03-13 [1] CRAN (R 3.6.3)
 memoise                1.1.0    2017-04-21 [1] CRAN (R 3.6.3)
 munsell                0.5.0    2018-06-12 [1] CRAN (R 3.6.3)
 nlme                   3.1-144  2020-02-06 [2] CRAN (R 3.6.3)
 nnet                   7.3-12   2016-02-02 [2] CRAN (R 3.6.3)
 pillar                 1.4.3    2019-12-20 [1] CRAN (R 3.6.3)
 pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 3.6.3)
 plotly                 4.9.2.1  2020-04-04 [1] CRAN (R 3.6.3)
 plyr                   1.8.6    2020-03-03 [1] CRAN (R 3.6.3)
 png                    0.1-7    2013-12-03 [1] CRAN (R 3.6.0)
 purrr                  0.3.3    2019-10-18 [1] CRAN (R 3.6.3)
 R6                     2.4.1    2019-11-12 [1] CRAN (R 3.6.3)
 raster                 3.0-12   2020-01-30 [1] CRAN (R 3.6.3)
 RColorBrewer           1.1-2    2014-12-07 [1] CRAN (R 3.6.0)
 Rcpp                   1.0.4    2020-03-17 [1] CRAN (R 3.6.3)
 RCurl                  1.98-1.1 2020-01-19 [1] CRAN (R 3.6.2)
 reshape                0.8.8    2018-10-23 [1] CRAN (R 3.6.3)
 rgeos                  0.5-2    2019-10-03 [1] CRAN (R 3.6.3)
 rlang                  0.4.5    2020-03-01 [1] CRAN (R 3.6.3)
 rpart                  4.1-15   2019-04-12 [2] CRAN (R 3.6.3)
 RPostgreSQL            0.6-2    2017-06-24 [1] CRAN (R 3.6.3)
 rprojroot              1.3-2    2018-01-03 [1] CRAN (R 3.6.3)
 RSQLite                2.2.0    2020-01-07 [1] CRAN (R 3.6.3)
 rstudioapi             0.11     2020-02-07 [1] CRAN (R 3.6.3)
 rvcheck                0.1.8    2020-03-01 [1] CRAN (R 3.6.3)
 S4Vectors              0.24.3   2020-01-18 [1] Bioconductor  
 scales                 1.1.0    2019-11-18 [1] CRAN (R 3.6.3)
 sessioninfo            1.1.1    2018-11-05 [1] CRAN (R 3.6.3)
 sf                     0.8-1    2020-01-28 [1] CRAN (R 3.6.3)
 snakecase              0.11.0   2019-05-25 [1] CRAN (R 3.6.3)
 sp                     1.4-1    2020-02-28 [1] CRAN (R 3.6.3)
 stringi                1.4.6    2020-02-17 [1] CRAN (R 3.6.2)
 stringr                1.4.0    2019-02-10 [1] CRAN (R 3.6.3)
 SummarizedExperiment   1.16.1   2019-12-20 [1] Bioconductor  
 survival               3.1-8    2019-12-03 [2] CRAN (R 3.6.3)
 tibble                 2.1.3    2019-06-06 [1] CRAN (R 3.6.2)
 tidyr                  1.0.2    2020-01-24 [1] CRAN (R 3.6.3)
 tidyselect             1.0.0    2020-01-27 [1] CRAN (R 3.6.3)
 tidytree               0.3.2    2020-03-12 [1] CRAN (R 3.6.3)
 treeio                 1.10.0   2019-10-29 [1] Bioconductor  
 units                  0.6-6    2020-03-16 [1] CRAN (R 3.6.3)
 vctrs                  0.2.4    2020-03-10 [1] CRAN (R 3.6.3)
 viridisLite            0.3.0    2018-02-01 [1] CRAN (R 3.6.3)
 withr                  2.1.2    2018-03-15 [1] CRAN (R 3.6.3)
 xfun                   0.12     2020-01-13 [1] CRAN (R 3.6.3)
 XML                    3.99-0.3 2020-01-20 [1] CRAN (R 3.6.2)
 xml2                   1.2.5    2020-03-11 [1] CRAN (R 3.6.3)
 xtable                 1.8-4    2019-04-21 [1] CRAN (R 3.6.3)
 XVector                0.26.0   2019-10-29 [1] Bioconductor  
 zlibbioc               1.32.0   2019-10-29 [1] Bioconductor  

[1] C:/Users/user/Documents/R/win-library/3.6
[2] C:/Program Files/R/R-3.6.3/library

```


<br>
