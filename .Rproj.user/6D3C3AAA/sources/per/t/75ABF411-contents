# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figure X for function of metagenome across samples and methods
# This script is written in R version 3.6.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('here')) install.packages('here'): library('here')
if (!require('pheatmap')) install.packages('pheatmap'): library('pheatmap')
if (!require('EnhancedVolcano')) BiocManager::install('EnhancedVolcano'): library('EnhancedVolcano')
if (!require("DEGreport")) BiocManager::install("DEGreport"): library('DEGreport')
if (!require("reshape2")) install.packages("reshape2"): library('reshape2')
if (!require('DESeq2')) install.packages('DESeq2'): library('DESeq2')


df_humann <- read.table('./Figures/Figure_7/go_annotation_notnorm.tsv',
                         row.names=1, header=T, sep='\t', 
                         comment.char='', quote='')

# Remove the suffix from the sample labels
names(df_humann) <- gsub('_Abundance.RPKs', '', names(df_humann))

hives_metadata <- read.csv("./Figures/Figure_3/Data_fig_3/metadata_hives.csv")
hives_metadata_new <- hives_metadata %>% column_to_rownames(var="X")


df_humann2 <- round(df_humann, digits = 0)

Hives_dds_functional <- DESeqDataSetFromMatrix(countData = df_humann2, colData = hives_metadata_new, design = ~Method + Season, tidy = F)
Hives_dds_RLE_functional <- estimateSizeFactors(Hives_dds_functional,type = "ratio")
Hives_normalised_counts_functional <- counts(Hives_dds_RLE_functional, normalized = TRUE)
Hives_counts_vst_functional <- varianceStabilizingTransformation(Hives_dds_RLE_functional, blind = FALSE)

ordered_counts_functional <- order(rowMeans(counts(Hives_dds_RLE_functional,normalized=TRUE)), decreasing=TRUE)

# Draw heatmap of all go functions
final_heatmap_functional <- pheatmap(assay(Hives_counts_vst_functional)[ordered_counts_functional,],cluster_rows = FALSE, show_rownames = FALSE, clustering_distance_cols = "correlation", cluster_cols = TRUE)

# Volcano plot to evaluate whether some differ by method 
get_volcano_plot <- function(dds_rle_object, .contrast) {
  object <- DESeq(dds_rle_object,test = "LRT", reduced = ~Method, full = ~Method+Season)
  res_dds <- results(object, .contrast)
  res <- lfcShrink(object,contrast = .contrast, res=res_dds,type = 'normal')
  EnhancedVolcano(res, lab = NA, x = 'log2FoldChange',y= 'pvalue', FCcutoff = 1,
                  colAlpha = 2,labSize = 2, xlim = c(-2,2))
}

get_volcano_plot(Hives_dds_RLE_functional, .contrast = c('Method','SM','Direct_SM'))

# By season
get_volcano_plot(Hives_dds_RLE_functional, .contrast = c('Season','May','November'))


# Calculate difference between SM and DSM for each hive individually
df_humann_matrix <- as.data.frame(matrix_humann2)
Hive_4_fun <- df_humann_matrix$DSM_H4 - df_humann_matrix$SM_H4
Hive_5_fun <- df_humann_matrix$DSM_H5 - df_humann_matrix$SM_H5
Hive_6_fun <- df_humann_matrix$DSM_H6 - df_humann_matrix$SM_H6
Hive_7_fun <- df_humann_matrix$DSM_H7 - df_humann_matrix$SM_H7


differences <- as.data.frame(cbind(Hive_4_fun,Hive_5_fun,Hive_6_fun,Hive_7_fun))
differences$pathways <- rownames(df_humann_matrix)

differences_for_boxplot <- melt(differences)

differences_boxplot <-  ggplot(differences_for_boxplot, aes(x=pathways, y=value)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11),
    axis.text.y = element_blank()
    ) + xlab("") + ylab ("DSM vs SM difference (%)") + coord_flip()


# Per species FUNCTIONAL diversity 
go_names <- read.table('./Figures/Figure_7/combined_tables_with_go_names.tsv',
                        row.names=1, header=T, sep='\t', 
                        comment.char='', quote='')


# Remove the suffix from the sample labels
names(go_names) <- gsub('_Abundance.RPKs', '', names(go_names))

go_names_deseq <- round(go_names, digits = 0)

Hives_dds_go <- DESeqDataSetFromMatrix(countData = go_names_deseq, colData = hives_metadata_new, design = ~Method + Season, tidy = F)
Hives_dds_RLE_go <- estimateSizeFactors(Hives_dds_go,type = "ratio")
Hives_normalised_counts_go <- counts(Hives_dds_RLE_go, normalized = TRUE)

# Custom function to create dataframe with per species functional diversity so that we can
# later plot it .
species_functional_diversity <- function(normalised_counts) {
  go_names_normalised <- as.data.frame(normalised_counts) %>% rownames_to_column(var="species")
  go_names_normalised$go_function <- go_names_normalised$species 
  go_names_normalised$species <- gsub(".*[__]([^.]+)[.].*", "\\1", go_names_normalised$species)
  go_names_normalised$go_function <- gsub("\\|.*", "", go_names_normalised$go_function)
  go_names_normalised <- go_names_normalised[!grepl("GO:", go_names_normalised$species),]
  go_names_normalised <- go_names_normalised[!grepl("UNMAPPED", go_names_normalised$species),]
  go_names_normalised <- go_names_normalised[!grepl("UNGROUPED", go_names_normalised$species),]
  go_names_normalised <- go_names_normalised[!grepl("UNGROUPED", go_names_normalised$go_function),]
}


functional_diversity_samples <- species_functional_diversity(Hives_normalised_counts_go)

# Custom function to return top X go_functions
get_most_abundant_go <- function(normalised_counts, top_x_gos) {
  go_names_normalised_x <- as.data.frame(normalised_counts) %>% rownames_to_column(var="species")
  go_names_normalised_x <- go_names_normalised_x[!grepl("g_", go_names_normalised_x$species),]
  go_names_normalised_x$rowsum <- rowSums(go_names_normalised_x[,-1])
  go_names_normalised_x <- go_names_normalised_x[-c(1,2),]
  go_names_normalised_x <- go_names_normalised_x[order(go_names_normalised_x$rowsum, decreasing = T),]
  go_top_x <- go_names_normalised_x$species[top_x_gos]
}

top_20_gos <- get_most_abundant_go(Hives_normalised_counts_go, top_x_gos = c(1:20))

go_for_plot <- reshape2::melt(functional_diversity_samples)

go_for_plot3 <- go_for_plot2[go_for_plot2$go_function %in% top_20_gos,]

functional_diversity_per_species <- ggplot(go_for_plot3, aes(x=species, y=value, fill = go_function)) + geom_bar(stat = "identity", position = "fill") + theme_minimal() + theme(legend.position="none", axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, size = 10), axis.title.x = element_blank()) + scale_fill_viridis_d(option = "magma") + ylab("GO function (%)") + labs(fill = "GO function") + scale_y_continuous(labels = c(0,25,50,75,100))
# Process previous honey samples
go_names_previous <- read.table('./Figures/Figure_7/combined_tables_with_go_names_previous.tsv',
                       row.names=1, header=T, sep='\t', 
                       comment.char='', quote='')


# Remove the suffix from the sample labels
names(go_names_previous) <- gsub('_Abundance.RPKs', '', names(go_names_previous))
previous_metadata <- read.csv("./Figures/Figure_7/metadata_previous.csv")
previous_metadata <- previous_metadata %>% column_to_rownames(var="X")

go_names_deseq_previous <- round(go_names_previous, digits = 0)


Hives_dds_go_previous <- DESeqDataSetFromMatrix(countData = go_names_deseq_previous, colData = previous_metadata, design = ~Source, tidy = F)
Hives_dds_RLE_go_previous <- estimateSizeFactors(Hives_dds_go_previous,type = "ratio")
Hives_normalised_counts_go_previous <- counts(Hives_dds_RLE_go_previous, normalized = TRUE)


previous_functional_diversity <- species_functional_diversity(go_names_deseq_previous)
top_30_previous <- get_most_abundant_go(Hives_normalised_counts_go_previous, top_x_gos = c(1:20))

# Get the 30 most abundant functions 


go_for_plot_previous <- reshape2::melt(previous_functional_diversity)

go_for_plot_previous2 <- go_for_plot_previous[go_for_plot_previous$go_function %in% top_30_previous,]

functional_diversity_per_species_previous <- ggplot(go_for_plot_previous2, aes(x=species, y=value, fill = go_function)) + geom_bar(stat = "identity", position = "fill") + theme_minimal() + theme(legend.position="none", axis.text.x = element_text(face = "italic", angle = 45, hjust = 1, size = 10), axis.title.x = element_blank()) + scale_fill_viridis_d(option = "magma") + ylab("GO function (%)") + labs(fill = "GO function") + scale_y_continuous(labels = c(0,25,50,75,100))


# Create venn diagram to compare the gos from our samples to the previously published ones
