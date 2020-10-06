# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figures 3A, 3B (30 most abundant species) and 3C but will also create similar Figures to 3B for Family and Genus levels
# This script is written in R version 3.6.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('taxonomizr')) install.packages('taxonomizr'); library('taxonomizr')
if (!require('here')) install.packages('here'): library('here')
# Taxonomizr will return the taxonomy for each species. However, this requires that a database is built locally (requires 60 GB of space).
# prepareDatabase('nameNode.sqlite')
# This process will take over 3 hours on a regular laptop/PC.  
# from here: and simply unzip it in the cloned repository. Place it at the top level, honeyDSM-seq and not in the subfolders.
if (!require('DESeq2')) install.packages('DESeq2'); library('DESeq2')
if (!require('pheatmap')) install.packages('pheatmap'): library('pheatmap')
if (!require('EnhancedVolcano')) BiocManager::install('EnhancedVolcano'): library('EnhancedVolcano')


data_path <- "./Figures/Figure_3/Data_fig_3"

# Vector with the column names of the kraken2 output
kraken2_output_names <- c("Reads_assigned_rooted_at_taxon", "Reads_assigned_directly_to_taxon","Rank_code","Taxonomic_ID","Name")

# Import the kraken2 output
import_kraken2_files <- list.files(path = data_path,pattern = "\\.kraken2",full.names = T)

kraken2_files <- lapply(import_kraken2_files, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,-c(1)]
})

naming_list <- list.files(path = data_path,pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

kraken2_files <- setNames(kraken2_files, substring(naming_list, first  = 1, last = nchar(naming_list) -8))

kraken2_files = lapply(kraken2_files,setNames,kraken2_output_names)

kraken2_files = lapply(kraken2_files,arrange, Taxonomic_ID)


# We will add the classification information. Plase note that this will take a little while to generate so be patient!
classification_ranks <- function(df) {
  taxids <- df[,4]
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum", 'nameNode.sqlite'))
  Kingdom <- as.data.frame(getTaxonomy(taxids, desiredTaxa = "kingdom", 'nameNode.sqlite'))
  Superkingdom <- as.data.frame(getTaxonomy(taxids, desiredTaxa = "superkingdom", 'nameNode.sqlite'))
  cbind(df, Phylum, Kingdom, Superkingdom)
}

kraken2_files = lapply(kraken2_files, classification_ranks)

# Filter out human and Drosophila contamination. For human contamination, the Chordata clade is removed to ensure that no false negatives appear.
kraken2_files_family <- kraken2_files %>% lapply(filter, Rank_code == "F") %>% lapply(filter, Name != "Drosophilidae") %>% lapply(filter, is.na(Phylum)|Phylum != "Chordata")
kraken2_files_genus <- kraken2_files %>% lapply(filter, Rank_code == "G") %>% lapply(filter, Name != "Drosophila") %>% lapply(filter, is.na(Phylum)|Phylum != "Chordata")
kraken2_files_species <- kraken2_files %>% lapply(filter, Rank_code == "S") %>% lapply(filter, Name != "Drosophila melanogaster") %>% lapply(filter, is.na(Phylum)|Phylum != "Chordata")

# This will generate a table with total reads per sample and level

get_totals_per_sample <- function(df) {
  classified <- df[2,1]
}

totals_classified <- lapply(kraken2_files,get_totals_per_sample) %>% do.call(rbind,.) %>% as.data.frame(.)

count_after_filter <- Hives_comparison_species[,-1]
totals_filter <- apply(count_after_filter,2,sum) %>% as.data.frame()

# merge the dataframes and calculate percent
class_filter <- cbind(totals_classified,totals_filter) 
colnames(class_filter) <- c("Total","Filtered")
class_filter$Percent <- round((class_filter$Filtered/class_filter$Total)*100, digits = 2) 


# This function will grab the number of reads assigned directly to the different levels
get_reads_per_level <- function(df) {
  df_1 <-  sum(df$Reads_assigned_directly_to_taxon);
  return(df_1)
}


family_per_sample <- lapply(kraken2_files_family,get_reads_per_level) %>% do.call(cbind,.) %>% as.data.frame()
genus_per_sample <-lapply(kraken2_files_genus,get_reads_per_level) %>% do.call(cbind,.) %>% as.data.frame()
species_per_sample <- lapply(kraken2_files_species,get_reads_per_level) %>% do.call(cbind,.) %>% as.data.frame()

all_reads_per_sample <- rbind(family_per_sample,genus_per_sample,species_per_sample)
all_reads_per_sample$Level <- c("Family","Genus","Species")


# This will only keep the taxonomic ID and number of reads from the table. For Genus and Family the rooted reads are used.
kraken2_files_filter_genus <- lapply(kraken2_files_genus, "[", c(1,4))
kraken2_files_filter_species <- lapply(kraken2_files_species, "[", c(2,4))
kraken2_files_filter_family <- lapply(kraken2_files_family, "[", c(1,4))

# The functions below rename the columns.
Hives_genus <- lapply(names(kraken2_files_filter_genus), function(x){
  colnames(kraken2_files_filter_genus[[x]]) <- c(x,"Taxonomic_ID")
  kraken2_files_filter_genus[[x]]
})
names(Hives_genus) <- names(kraken2_files_filter_genus) 

Hives_species <- lapply(names(kraken2_files_filter_species), function(x){
  colnames(kraken2_files_filter_species[[x]]) <- c(x,"Taxonomic_ID")
  kraken2_files_filter_species[[x]]
})
names(Hives_species) <- names(kraken2_files_filter_species) 


Hives_family <- lapply(names(kraken2_files_filter_family), function(x){
  colnames(kraken2_files_filter_family[[x]]) <- c(x,"Taxonomic_ID")
  kraken2_files_filter_family[[x]]
})
names(Hives_family) <- names(kraken2_files_filter_family) 

# The functions below collapse the list into a single data frame and all NA values are adjusted to 0
Hives_comparison_genus <- Hives_genus %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.),0))
Hives_comparison_species <- Hives_species %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))
Hives_comparison_family <- Hives_family %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

# Normalisation with DESeq2 requires a metadata file
hives_metadata <- read.csv("./Figures/Figure_3/Data_fig_3/metadata_hives.csv")

# We apply the RLE normalisation but using the Hive + Method to observe the method effect so we can use it for the tree figure
Hives_dds_family <- DESeqDataSetFromMatrix(countData = Hives_comparison_family, colData = hives_metadata, design = ~Hive + Method, tidy = TRUE)
Hives_dds_RLE_family <- estimateSizeFactors(Hives_dds_family,type = "ratio")
Hives_normalised_counts_family <- counts(Hives_dds_RLE_family, normalized = TRUE)
Hives_counts_vst_family <- varianceStabilizingTransformation(Hives_dds_RLE_family, blind = FALSE)


# Export
Hives_normalised_counts_family_export <- as.data.frame(Hives_normalised_counts_family)
Hives_normalised_counts_family_export <- Hives_normalised_counts_family_export %>% rownames_to_column(., var = "Taxonomic_ID")
write.csv(Hives_normalised_counts_family_export, file = "./Figures/Figure_3/normalised_counts_family.csv")


# Create the heatmap
Hives_annotation <- as.data.frame(colData(Hives_dds_family))

get_super_kingdom_or_plant <- function(x) {
  ifelse(getTaxonomy(x, desiredTaxa = "superkingdom", sqlFile = 'nameNode.sqlite') != "Eukaryota", getTaxonomy(x, desiredTaxa = "superkingdom", sqlFile = 'nameNode.sqlite'), ifelse(getTaxonomy(x, desiredTaxa = "kingdom", sqlFile = 'nameNode.sqlite') == "Viridiplantae", getTaxonomy(x, desiredTaxa = "kingdom",sqlFile =  'nameNode.sqlite'), getTaxonomy(x, desiredTaxa = "superkingdom", sqlFile = 'nameNode.sqlite')))
  }

TaxonomicIDs_family <- as.numeric(rownames((Hives_counts_vst_family)))
Taxid_taxonomy_family <- as.data.frame(get_super_kingdom_or_plant(TaxonomicIDs_family)) 
Taxonomic_IDs_dataframe_family <- as.data.frame(TaxonomicIDs_family)
annotationrows_family <- bind_cols(Taxid_taxonomy_family,Taxonomic_IDs_dataframe_family)
colnames(annotationrows_family) <- c("Domain", "Taxonomic_ID")
annotationrows_family <- tibble::column_to_rownames(annotationrows_family, var = "Taxonomic_ID")

ann_colors = list(
  Domain = c(Bacteria = "#0D0887FF", Viruses = "#73D055FF", Eukaryota = "#CC4678FF", Viridiplantae = "#F0F921FF", Archaea = "#ad05f5"),
  Season = c(May = "#5cb300", July = "#ecf547", November = "#8f7c56"))


ordered_counts_family <- order(rowMeans(counts(Hives_dds_RLE_family,normalized=TRUE)), decreasing=TRUE)

final_heatmap_family <- pheatmap(assay(Hives_counts_vst_family)[ordered_counts_family,],cluster_rows = FALSE, show_rownames = FALSE, clustering_distance_cols = "correlation", cluster_cols = TRUE, annotation_col = select(Hives_annotation,Season,Varroa,Population), annotation_row = annotationrows_family, annotation_colors = ann_colors[1:2], main = "Heatmap families", height = 20, scale = 'none')
ggsave(filename="heatmap_family.png", plot=final_heatmap_family, path = "./Figures/Figure_3")

# The same functions as above are applied to Genus and Species levels

Hives_dds_genus <- DESeqDataSetFromMatrix(countData = Hives_comparison_genus, colData = hives_metadata, design = ~Method + Hive, tidy = TRUE)
Hives_dds_RLE_genus <- estimateSizeFactors(Hives_dds_genus,type = "ratio")
Hives_normalised_counts_genus <- counts(Hives_dds_RLE_genus, normalized = TRUE)
Hives_counts_vst_genus <- varianceStabilizingTransformation(Hives_dds_RLE_genus, blind = FALSE)

Hives_annotation <- as.data.frame(colData(Hives_dds_genus))

TaxonomicIDs_genus <- as.numeric(rownames((Hives_counts_vst_genus)))
Taxid_taxonomy_genus <- as.data.frame(get_super_kingdom_or_plant(TaxonomicIDs_genus)) 
Taxonomic_IDs_dataframe_genus <- as.data.frame(TaxonomicIDs_genus)
annotationrows_genus <- bind_cols(Taxid_taxonomy_genus,Taxonomic_IDs_dataframe_genus)
colnames(annotationrows_genus) <- c("Domain", "Taxonomic_ID")
annotationrows_genus <- tibble::column_to_rownames(annotationrows_genus, var = "Taxonomic_ID")

ann_colors = list(
  Domain = c(Bacteria = "#0D0887FF", Viruses = "#73D055FF", Eukaryota = "#CC4678FF", Viridiplantae = "#F0F921FF", Archaea = "#ad05f5"),
  Season = c(May = "#5cb300", July = "#ecf547", November = "#8f7c56"))


ordered_counts_genus <- order(rowMeans(counts(Hives_dds_RLE_genus,normalized=TRUE)), decreasing=TRUE)
final_heatmap_genus <- pheatmap(assay(Hives_counts_vst_genus)[ordered_counts_genus,],cluster_rows = FALSE, show_rownames = FALSE, clustering_distance_cols = "correlation", cluster_cols = TRUE, annotation_col = select(Hives_annotation,"Season"), annotation_row = annotationrows_genus, annotation_colors = ann_colors[1:2], main = "Heatmap genera", height = 15)
ggsave(filename="heatmap_genus.png", plot=final_heatmap_genus, path ="./Figures/Figure_3")

# Species
Hives_dds_species <- DESeqDataSetFromMatrix(countData = Hives_comparison_species, colData = hives_metadata, design = ~Method + Hive, tidy = TRUE)
Hives_dds_RLE_species <- estimateSizeFactors(Hives_dds_species,type = "ratio")
Hives_normalised_counts_species <- counts(Hives_dds_RLE_species, normalized = TRUE)
Hives_counts_vst_species <- varianceStabilizingTransformation(Hives_dds_RLE_species, blind = FALSE)

# Export the counts because it will be needed for the next figure
Hives_normalised_counts_species_export <- as.data.frame(Hives_normalised_counts_species)
Hives_normalised_counts_species_export <- Hives_normalised_counts_species_export %>% rownames_to_column(., var = "Taxonomic_ID")
write.csv(Hives_normalised_counts_species_export, file = "./Figures/Figure_3/normalised_counts_species.csv")


Hives_annotation <- as.data.frame(colData(Hives_dds_species))

TaxonomicIDs_species <- as.numeric(rownames((Hives_counts_vst_species)))
Taxid_taxonomy_species <- as.data.frame(get_super_kingdom_or_plant(TaxonomicIDs_species)) 
Taxonomic_IDs_dataframe_species <- as.data.frame(TaxonomicIDs_species)

generate_species <- getTaxonomy(TaxonomicIDs_species,desiredTaxa = "species", 'nameNode.sqlite')
generate_species = as.data.frame(generate_species)
generate_species$species = as.character(generate_species$species)
generate_species_omit_na <- na.omit(generate_species)
generate_species_omit_na <- rownames_to_column(generate_species_omit_na, var = "Taxonomic_ID")

# 4 species were not identified and therefore we have to replace manually by looking up to the NCBI taxonomy database.
# I create a vector so that it is clear which Taxnomic IDs I am going to replace.
unidentified_species <- subset(generate_species,is.na(generate_species$species))
unidentified_species_vector <- rownames(unidentified_species)
unidentified_species_names_ncbi <- c("Streptomyces sp. WAC08241","Actinomadura sp. WMMB 499","Exaiptasia diaphana","Pseudomonas mediterranea","Xanthomonas hortorum pv. gardneri")
unidentified_species$species <- unidentified_species_names_ncbi
unidentified_species <- rownames_to_column(unidentified_species, var = "Taxonomic_ID")

# And now we will add them back to the original dataframe
generate_species_v2 <- rbind(generate_species_omit_na,unidentified_species)


annotationrows_species <- bind_cols(Taxid_taxonomy_species,Taxonomic_IDs_dataframe_species,generate_species_v2)
annotationrows_species <- annotationrows_species[,-c(2,3)]
colnames(annotationrows_species) <- c("Domain","Species")


# Fix NA values. For some taxonomizr isn't able to return the superkingdom so we simply remove them from inclusion in the heatmap (approx 10 species)
na_values <- subset(annotationrows_species, is.na(annotationrows_species$Domain))
na_values_tax <- getTaxonomy(na_values$Taxonomic_ID, desiredTaxa = "superkingdom",'nameNode.sqlite')
na_values_tax_removed_na <- na.omit(na_values_tax)
na_values_tax <- as.data.frame(na_values_tax) %>% rownames_to_column(var = "Taxonomic_ID")
na_values_tax_removed_na$Taxonomic_ID <- as.numeric(na_values_tax_removed_na$Taxonomic_ID)
colnames(na_values_tax) <- c("Taxonomic_ID","Domain")

annotationrows_species <- tibble::column_to_rownames(annotationrows_species, var = "Species")


ann_colors = list(
  Domain = c(Bacteria = "#0D0887FF", Viruses = "#F0F921FF", Eukaryota = "#CC4678FF", Viridiplantae = "#73D055FF", Archaea = "#ad05f5"),
  Season = c(May = "#5cb300", July = "#ecf547", November = "#8f7c56"))

rownames(Hives_counts_vst_species) <- generate_species_v2$species
rownames(annotationrows_species) <- generate_species_v2$species

ordered_counts_species <- order(rowMeans(counts(Hives_dds_RLE_species,normalized=TRUE)), decreasing=TRUE)[1:30]
final_heatmap_species <- pheatmap(assay(Hives_counts_vst_species)[ordered_counts_species,],clustering_distance_cols = "correlation", cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = TRUE, annotation_col = select(Hives_annotation,Season), labels_col = c("DirectSM H4","DirectSM H5","DirectSM H6" ,"DirectSM H7","SM H4", "SM H5","SM H6","SM H7"), annotation_row = annotationrows_species, annotation_colors = ann_colors[1:2], height = 100, fontsize = 20)

ggsave(filename="top30_heatmap_species.pdf",device="pdf", plot=final_heatmap_species, height = 13, width = 15, path = "./Figures/Figure_3")



# Create MA plots for method
Hives_dds_species_v2 <- DESeqDataSetFromMatrix(countData = Hives_comparison_species, colData = hives_metadata, design = ~Method + Hive, tidy = TRUE)
Hives_dds_RLE_species_v2 <- estimateSizeFactors(Hives_dds_species_v2,type = "ratio")

Hives_dds_genus_v2 <- DESeqDataSetFromMatrix(countData = Hives_comparison_genus, colData = hives_metadata, design = ~Method + Hive, tidy = TRUE)
Hives_dds_RLE_genus_v2 <- estimateSizeFactors(Hives_dds_genus_v2,type = "ratio")

Hives_dds_family_v2 <- DESeqDataSetFromMatrix(countData = Hives_comparison_family, colData = hives_metadata, design = ~Method + Hive, tidy = TRUE)
Hives_dds_RLE_family_v2 <- estimateSizeFactors(Hives_dds_family_v2,type = "ratio")

Hives_normalised_counts_species_v2 <- counts(Hives_dds_RLE_species_v2, normalized = TRUE)
Hives_normalised_counts_family_v2 <- counts(Hives_dds_RLE_family_v2, normalized = TRUE)
Hives_normalised_counts_genus_v2 <- counts(Hives_dds_RLE_genus_v2, normalized = TRUE)

# Export the counts because it will be needed for the next figure
Hives_normalised_counts_family_v2_exprort <- as.data.frame(Hives_normalised_counts_family_v2)
Hives_normalised_counts_family_v2_exprort <- Hives_normalised_counts_family_v2_exprort %>% rownames_to_column(., var = "Taxonomic_ID")
write.csv(Hives_normalised_counts_family_v2_exprort, file = "./Figures/Figure_3/normalised_counts_methodhive_family.csv")

# Also export normalised based on method and season to check seasonal fluctuations 
Hives_normalised_counts_species_v2_export <- as.data.frame(Hives_normalised_counts_species_v2)
Hives_normalised_counts_species_v2_export <- Hives_normalised_counts_species_v2_export %>% rownames_to_column(., var = "Taxonomic_ID")
write.csv(Hives_normalised_counts_species_v2_export, file = "./Figures/Figure_6/Data_fig_6/normalised_counts_methodhive_species.csv")

# Export genus
# Export the counts because it will be needed for the next figure
Hives_normalised_counts_genus_v2_exprort <- as.data.frame(Hives_normalised_counts_genus_v2)
Hives_normalised_counts_genus_v2_exprort <- Hives_normalised_counts_genus_v2_exprort %>% rownames_to_column(., var = "Taxonomic_ID")
write.csv(Hives_normalised_counts_genus_v2_exprort, file = "./Figures/Figure_3/normalised_counts_methodhive_genus.csv")

# Make an annotation table so we can adjust the circle size
annotation_volcano_family <- Taxid_taxonomy_family %>% rownames_to_column(var = "Taxonomic_ID")
annotation_volcano_family$Taxonomic_ID <- as.numeric(annotation_volcano_family$Taxonomic_ID)


# draw volcano plots
get_volcano_plot <- function(dds_rle_object, contrast, annotation) {
  object <- DESeq(dds_rle_object)
  res_dds <- results(object, contrast)
  res <- lfcShrink(object,contrast = contrast, res=res_dds,type = 'normal')
  annotations <- annotation %>% rownames_to_column(var = "Taxonomic_ID")
  annotations$Taxonomic_ID <- as.numeric(annotations$Taxonomic_ID)
  EnhancedVolcano(res, lab = NA, x = 'log2FoldChange',y= 'pvalue', FCcutoff = 1,
                  colAlpha = 2, pointSize = c(ifelse(annotations$superkingdom[match(rownames(res),annotations$Taxonomic_ID)] == "Viridiplantae", 3,1))
                  ,labSize = 2, xlim = c(-2,2))
}


family_volcano_plot <- get_volcano_plot(Hives_dds_RLE_family_v2, contrast = c('Method','Direct_SM','SM'), annotation = Taxid_taxonomy_family)
genus_volcano_plot <- get_volcano_plot(Hives_dds_RLE_genus_v2, contrast = c('Method','Direct_SM','SM'), annotation = Taxid_taxonomy_genus)
species_volcano_plot <- get_volcano_plot(Hives_dds_RLE_species_v2, contrast = c('Method','Direct_SM','SM'), annotation = Taxid_taxonomy_species)


# Save the MA plots
ggsave(plot = family_volcano_plot, filename = "./Figures/Figure_3/volcano_plot_family.pdf")
ggsave(plot = genus_volcano_plot, filename = "./Figures/Figure_3/volcano_plot_genus.pdf")
ggsave(plot = species_volcano_plot, filename = "./Figures/Figure_3/volcano_plot_species.pdf")


# Classification function for the DESeq2 object
classification_deseq <- function(df) {
  taxids <- df$Taxonomic_ID
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum",'nameNode.sqlite'))
  Superkingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "superkingdom",'nameNode.sqlite'))
  Kingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "kingdom", 'nameNode.sqlite'))
  Class <- as.character(getTaxonomy(taxids, desiredTaxa = "class", 'nameNode.sqlite'))
  Order <- as.character(getTaxonomy(taxids, desiredTaxa = "order", 'nameNode.sqlite'))
  Family <- as.character(getTaxonomy(taxids,desiredTaxa = "family", 'nameNode.sqlite'))
  Genus <- as.character(getTaxonomy(taxids,desiredTaxa = "genus", 'nameNode.sqlite'))
  Species <- as.character(getTaxonomy(taxids,desiredTaxa = "species", sqlFile = 'nameNode.sqlite'))
  cbind(df, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, Species)
}


# Plot differentially abundant species per season
Hives_dds_species_v3 <- DESeqDataSetFromMatrix(countData = Hives_comparison_species, colData = hives_metadata, design = ~Season, tidy = TRUE)
Hives_dds_RLE_species_v3 <- estimateSizeFactors(Hives_dds_species_v2,type = "ratio")

Hives_dds_genus_v3 <- DESeqDataSetFromMatrix(countData = Hives_comparison_genus, colData = hives_metadata, design = ~Method + Season, tidy = TRUE)
Hives_dds_RLE_genus_v3 <- estimateSizeFactors(Hives_dds_genus_v2,type = "ratio")

Hives_dds_family_v3 <- DESeqDataSetFromMatrix(countData = Hives_comparison_family, colData = hives_metadata, design = ~Method + Season, tidy = TRUE)
Hives_dds_RLE_family_v3 <- estimateSizeFactors(Hives_dds_family_v2,type = "ratio")


# Function to draw the significantly differentially abundant species/genera/families/whatever you throw at it
plotDiffAbund <- function(colNums, DESeq_RLE_object, title, level = c("species","genus","family","order","phylum")) {
  
  dds_object <- DESeq(Hives_dds_genus_v3, test = "LRT", reduced = ~1)
  rld <- rlog(dds_object, blind=F)
  results <- subset(results(dds_object), padj < 0.05)
  
  
  # make the lists
  upgenes <- rownames(head(results[ order( results$log2FoldChange ), ], n=40))
  downgenes <- rownames(head(results[ order( -results$log2FoldChange ), ], n=40))
  
  # this gives us the rows we want
  rows <- match(upgenes, row.names(rld))
  mat <- assay(rld)[rows,c(1:8)]
  mat <- mat - rowMeans(mat)
  
  TaxonomicIDs_species <- as.numeric(rownames((mat)))
  Taxid_taxonomy_species <- as.data.frame(get_super_kingdom_or_plant(TaxonomicIDs_species)) 
  Taxonomic_IDs_dataframe_species <- as.data.frame(TaxonomicIDs_species)
  
  generate_species <- getTaxonomy(TaxonomicIDs_species,desiredTaxa = level, sqlFile = 'nameNode.sqlite')
  generate_species = as.data.frame(generate_species)
  generate_species[,1] = as.character(generate_species[,1])
  annotationrows_species1 <- bind_cols(Taxid_taxonomy_species,Taxonomic_IDs_dataframe_species,generate_species)
  colnames(annotationrows_species1) <- c("Domain", "Taxonomic_ID", level)
  
  annotationrows_species1 = annotationrows_species1[,-2]
  annotationrows_species1 <- tibble::column_to_rownames(annotationrows_species1, var = level)
  rownames(mat) <- generate_species[,1]
  
  
  
  # the labels are hard coded at the moment :(
  df <- as.data.frame(colData(rld)[c("Season")])
  pheatmap(mat, fontsize=15, annotation_colors = ann_colors, annotation_row = annotationrows_species1, cluster_rows = FALSE, show_rownames = TRUE, height = 100, angle_col = "270", cluster_cols = TRUE, show_colnames = TRUE, border_color = NA, annotation_col = df, labels_col = c("DirectSM H4","DirectSM H5","DirectSM H6" ,"DirectSM H7","SM H4", "SM H5","SM H6","SM H7"), scale = 'row')
}


# Normalise using method and season to get the species that are variable by season
Hives_dds_species_v3 <- DESeqDataSetFromMatrix(countData = Hives_comparison_species, colData = hives_metadata, design = ~Method + Season, tidy = TRUE)
Hives_dds_RLE_species_v3 <- estimateSizeFactors(Hives_dds_species_v3,type = "ratio")
Hives_dds_genus_v3 <- DESeqDataSetFromMatrix(countData = Hives_comparison_genus, colData = hives_metadata, design = ~Season, tidy = TRUE)
Hives_dds_RLE_genus_v3 <- estimateSizeFactors(Hives_dds_genus_v3,type = "ratio")
Hives_dds_family_v3 <- DESeqDataSetFromMatrix(countData = Hives_comparison_family, colData = hives_metadata, design = ~Method + Season, tidy = TRUE)
Hives_dds_RLE_family_v3 <- estimateSizeFactors(Hives_dds_family_v3,type = "ratio")


differentially_abundant_species <- plotDiffAbund(
  DESeq_RLE_object =Hives_dds_species_v3,
  colNums = c(1:8),
  level = "species")

differentially_abundant_genera <- plotDiffAbund(
  DESeq_RLE_object =Hives_dds_genus_v3,
  colNums = c(1:8),
  level = "genus")

differentially_abundant_families <- plotDiffAbund(
  DESeq_RLE_object =Hives_dds_family_v3,
  colNums = c(1:8),
  level = "family")


ggsave(plot= differentially_abundant_species,filename = "./Figures/Figure_3/significantly_abundant_species.pdf",device="pdf", height = 10, width = 14)
ggsave(plot= differentially_abundant_genera,filename = "./Figures/Figure_3/significantly_abundant_genera.pdf",device="pdf", height = 10, width = 14)
ggsave(plot= differentially_abundant_families,filename = "./Figures/Figure_3/significantly_abundant_families.pdf",device="pdf", height = 10, width = 14)


