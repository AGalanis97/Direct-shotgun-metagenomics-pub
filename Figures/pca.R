if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('taxonomizr')) install.packages('taxonomizr'); library('taxonomizr')
if (!require('here')) install.packages('here'): library('here')
# Taxonomizr will return the taxonomy for each species. However, this requires that a database is built locally (requires 60 GB of space).
# prepareDatabase('nameNode.sqlite')
# This process will take over 3 hours on a regular laptop/PC.  
# from here: and simply unzip it in the cloned repository. Place it at the top level, honeyDSM-seq and not in the subfolders.
if (!require('DESeq2')) install.packages('DESeq2'); library('DESeq2')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
library(factoextra)
library(prcomp)

normalised_counts_methodseason <- read.csv("./Figures/Figure_3/normalised_methodseason_family.csv")[,-1]
normalised_counts_pca <- normalised_counts_methodseason %>% column_to_rownames(var = "Taxonomic_ID")
normalised_counts_pca <- t(normalised_counts_pca[,c(1:8)])

res.pca <- prcomp(normalised_counts_pca, scale. = F)

# Visualise dimensions
fviz_eig(res.pca)

biplot <- factoextra::fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", select.var = list(contrib = 20)
)

res.pca$x

fviz_contrib(res.pca, choice = "var", axes = 2, top = 20)
