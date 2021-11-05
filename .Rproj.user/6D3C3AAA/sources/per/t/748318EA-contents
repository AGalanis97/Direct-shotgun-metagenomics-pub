# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.
# This script is written in R version 4.0.3

# This script will reproduce Figures 3 A, B, D 
# It also produces data and figures for family, genus, and species level

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('taxonomizr')) install.packages('taxonomizr'); library('taxonomizr')
if (!require('here')) install.packages('here'): library('here')
# Taxonomizr will return the taxonomy for each species. However, this requires that a database is built locally (requires 60 GB of space).
# prepareDatabase('nameNode.sqlite')
# This process will take over 3 hours on a regular laptop/PC.  
# from here: and simply unzip it in the cloned repository. Place it at the top level, honeyDSM-seq and not in the subfolders.
if (!require('DESeq2')) BiocManager::install('DESeq2'): library('DESeq2')
if (!require('waffle')) devtools::install_github("hrbrmstr/waffle"); library('waffle')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('pheatmap')) install.packages('pheatmap'): library('pheatmap')
if (!require('EnhancedVolcano')) BiocManager::install('EnhancedVolcano'): library('EnhancedVolcano')
if (!require("DEGreport")) BiocManager::install("DEGreport"): library('DEGreport')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('factoextra')) install.packages('factoextra'); library('factoextra')

# if (!require("lasso2")) install.packages("lasso2"): library('lasso2')
# Requires R 4.0 now. Install from archive.

data_path <- "./Figures/Figure_3/Data_fig_3"
output_path <- "./Output_data/Figure_3_output/"

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
Hives_comparison_genus <- Hives_genus %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% dplyr::select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.),0))
Hives_comparison_species <- Hives_species %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% dplyr::select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))
Hives_comparison_family <- Hives_family %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% dplyr::select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

#Hives_comparison_genus <- Hives_comparison_species
#Hives_comparison_genus$Taxonomic_ID <- as.character(Hives_comparison_species$Taxonomic_ID)
#Hives_comparison_genus <- Hives_comparison_genus %>% mutate_if(is.numeric, funs(ifelse(. < 50,0,.)))
#Hives_comparison_genus <- Hives_comparison_genus %>% filter(rowSums(.[,-1]) > 0)

# Before we proceed with normalisation we create a function for taxonomic classification before we export the normalised reads
# Classification function for the DESeq2 object to export the normalised reads
classification_deseq_export <- function(normalised_reads, path_filename) {
  df <- normalised_reads %>% round(.,digits=0) %>% as.data.frame() %>% rownames_to_column(var = "Taxonomic_ID")
  taxids <- df$Taxonomic_ID
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum",'nameNode.sqlite'))
  Superkingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "superkingdom",'nameNode.sqlite'))
  Kingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "kingdom", 'nameNode.sqlite'))
  Class <- as.character(getTaxonomy(taxids, desiredTaxa = "class", 'nameNode.sqlite'))
  Order <- as.character(getTaxonomy(taxids, desiredTaxa = "order", 'nameNode.sqlite'))
  Family <- as.character(getTaxonomy(taxids,desiredTaxa = "family", 'nameNode.sqlite'))
  Genus <- as.character(getTaxonomy(taxids,desiredTaxa = "genus", 'nameNode.sqlite'))
  Species <- as.character(getTaxonomy(taxids,desiredTaxa = "species", sqlFile = 'nameNode.sqlite'))
  df_export <- cbind(df, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, Species)
  write.csv(df_export, file = path_filename)
}


# Normalisation with DESeq2 requires a metadata file
hives_metadata <- read.csv("./Figures/Figure_3/Data_fig_3/metadata_hives.csv", stringsAsFactors = T)

# We apply the RLE normalisation but using the Method + Season so we can observe seasonal effects and do methodological validation
Hives_dds_family <- DESeqDataSetFromMatrix(countData = Hives_comparison_family, colData = hives_metadata, design = ~Method + Season, tidy = TRUE)
Hives_dds_RLE_family <- estimateSizeFactors(Hives_dds_family,type = "ratio")

filtered_family <- rowSums( counts(Hives_dds_RLE_family, normalized = TRUE) >= 50) >=2
Hives_dds_RLE_family_filtered <- Hives_dds_RLE_family[filtered_family,]

Hives_normalised_counts_family <- counts(Hives_dds_RLE_family_filtered, normalized = TRUE)
Hives_counts_vst_family <- varianceStabilizingTransformation(Hives_dds_RLE_family_filtered, blind = FALSE)

# Calculate the p-value with the method as contrast for figure 3 tree
family_deseq <- DESeq(Hives_dds_RLE_family_filtered)
results_family_deseq <- DESeq2::results(family_deseq, contrast = c('Method','Direct_SM','SM'))
pvalues_family <- as.data.frame(results_family_deseq$padj) %>% rename("results_family_deseq$padj" = "padj")
write.csv(pvalues_family,file = "./Figures/Figure_3/padj_family_method.csv")

# Export 
classification_deseq_export(normalised_reads = Hives_normalised_counts_family, path_filename = "./Figures/Figure_3/normalised_methodseason_family_filtered.csv")

# This function will report whether an organism belongs to plants or not
get_super_kingdom_or_plant <- function(x) {
  ifelse(getTaxonomy(x, desiredTaxa = "superkingdom", sqlFile = 'nameNode.sqlite') != "Eukaryota", getTaxonomy(x, desiredTaxa = "superkingdom", sqlFile = 'nameNode.sqlite'), ifelse(getTaxonomy(x, desiredTaxa = "kingdom", sqlFile = 'nameNode.sqlite') == "Viridiplantae", getTaxonomy(x, desiredTaxa = "kingdom",sqlFile =  'nameNode.sqlite'), getTaxonomy(x, desiredTaxa = "superkingdom", sqlFile = 'nameNode.sqlite')))
}

# Create the heatmap
Hives_annotation <- as.data.frame(colData(Hives_dds_family))

TaxonomicIDs_family <- as.numeric(rownames((Hives_counts_vst_family)))
Taxid_taxonomy_family <- as.data.frame(get_super_kingdom_or_plant(TaxonomicIDs_family)) 
Taxonomic_IDs_dataframe_family <- as.data.frame(TaxonomicIDs_family)
annotationrows_family <- bind_cols(Taxid_taxonomy_family,Taxonomic_IDs_dataframe_family)
colnames(annotationrows_family) <- c("Domain", "Taxonomic_ID")
annotationrows_family <- tibble::column_to_rownames(annotationrows_family, var = "Taxonomic_ID")

ann_colors = list(
  Domain = c(Bacteria = "#E6AB02", Viruses = "#eb53a1", Eukaryota = "#7570B3", Viridiplantae = "#1B9E77", Archaea = "#D95F02"),
  Season = c(May = "#5cb300", July = "#e8df2a", November = "#8f7c56"))

ordered_counts_family <- order(rowMeans(counts(Hives_dds_RLE_family_filtered,normalized=TRUE)), decreasing=TRUE)

final_heatmap_family <- pheatmap(assay(Hives_counts_vst_family)[ordered_counts_family,],cluster_rows = FALSE, show_rownames = FALSE, clustering_distance_cols = "correlation", cluster_cols = TRUE, annotation_col = select(Hives_annotation,Season), annotation_row = annotationrows_family, annotation_colors = ann_colors[1:2], main = "Heatmap families", height = 20, scale = 'none')

ggsave(filename="heatmap_family.png", plot=final_heatmap_family, path = output_path)

# The same functions as above are applied to Genus and Species levels
Hives_dds_genus <- DESeqDataSetFromMatrix(countData = Hives_comparison_genus, colData = hives_metadata, design = ~Method + Season, tidy = TRUE)
Hives_dds_RLE_genus <- estimateSizeFactors(Hives_dds_genus,type = "ratio")

filtered_genus <- rowSums( counts(Hives_dds_RLE_genus, normalized = TRUE) >= 50) >=2
Hives_dds_RLE_genus_filtered <- Hives_dds_RLE_genus[filtered_genus,]


Hives_normalised_counts_genus <- counts(Hives_dds_RLE_genus, normalized = TRUE)
Hives_counts_vst_genus <- varianceStabilizingTransformation(Hives_dds_RLE_genus_filtered, blind = FALSE)

# Export reads
classification_deseq_export(normalised_reads = Hives_normalised_counts_genus, path_filename = "./Figures/Figure_3/normalised_methodseason_genus.csv")


Hives_annotation <- as.data.frame(colData(Hives_dds_genus))

TaxonomicIDs_genus <- as.numeric(rownames((Hives_counts_vst_genus)))
Taxid_taxonomy_genus <- as.data.frame(get_super_kingdom_or_plant(TaxonomicIDs_genus)) 
Taxonomic_IDs_dataframe_genus <- as.data.frame(TaxonomicIDs_genus)
annotationrows_genus <- bind_cols(Taxid_taxonomy_genus,Taxonomic_IDs_dataframe_genus)
colnames(annotationrows_genus) <- c("Domain", "Taxonomic_ID")
annotationrows_genus <- tibble::column_to_rownames(annotationrows_genus, var = "Taxonomic_ID")

ann_colors = list(
  Domain = c(Bacteria = "#E6AB02", Viruses = "#eb53a1", Eukaryota = "#7570B3", Viridiplantae = "#1B9E77", Archaea = "#D95F02"),
  Season = c(May = "#5cb300", July = "#e8df2a", November = "#8f7c56"))


ordered_counts_genus <- order(rowMeans(counts(Hives_dds_RLE_genus_filtered,normalized=TRUE)), decreasing=TRUE)
final_heatmap_genus <- pheatmap(assay(Hives_counts_vst_genus)[ordered_counts_genus,],cluster_rows = FALSE, show_rownames = FALSE, clustering_distance_cols = "correlation", cluster_cols = TRUE, annotation_col = select(Hives_annotation,"Season"), annotation_row = annotationrows_genus, annotation_colors = ann_colors[1:2], main = "Heatmap genera", height = 15)
ggsave(filename="heatmap_genus.png", plot=final_heatmap_genus, path =output_path)


# Species
Hives_dds_species <- DESeqDataSetFromMatrix(countData = Hives_comparison_species, colData = hives_metadata, design = ~Method + Season, tidy = TRUE)

Hives_dds_RLE_species <- estimateSizeFactors(Hives_dds_species,type = "ratio")

filtered_species <- rowSums( counts(Hives_dds_RLE_species, normalized = TRUE) >= 50) >=2
Hives_dds_RLE_species_filtered <- Hives_dds_RLE_species[filtered_species,]


Hives_normalised_counts_species <- counts(Hives_dds_RLE_species, normalized = TRUE)
Hives_counts_vst_species <- varianceStabilizingTransformation(Hives_dds_RLE_species, blind = FALSE)

# Export the counts because it will be needed for the next figure
classification_deseq_export(normalised_reads = Hives_normalised_counts_species, path_filename = "./Figures/Figure_3/normalised_methodseason_species.csv")

Hives_annotation <- as.data.frame(colData(Hives_dds_species))

TaxonomicIDs_species <- as.numeric(rownames((Hives_counts_vst_species_filtered)))
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
#unidentified_species_vector <- rownames(unidentified_species)
#unidentified_species_names_ncbi <- c("Streptomyces sp. WAC08241","Actinomadura sp. WMMB 499","Exaiptasia diaphana","Pseudomonas mediterranea","Xanthomonas hortorum pv. gardneri")
#unidentified_species$species <- unidentified_species_names_ncbi
#unidentified_species <- rownames_to_column(unidentified_species, var = "Taxonomic_ID")

# And now we will add them back to the original dataframe
generate_species_v2 <- rbind(generate_species_omit_na,unidentified_species)


annotationrows_species <- bind_cols(Taxid_taxonomy_species,Taxonomic_IDs_dataframe_species,generate_species_v2)
annotationrows_species <- annotationrows_species[,-c(2,3)]
colnames(annotationrows_species) <- c("Domain","Species")
annotationrows_species <- tibble::column_to_rownames(annotationrows_species, var = "Species")


ann_colors = list(
  Domain = c(Bacteria = "#E6AB02", Viruses = "#eb53a1", Eukaryota = "#7570B3", Viridiplantae = "#1B9E77", Archaea = "#D95F02"),
  Season = c(May = "#5cb300", July = "#e8df2a", November = "#8f7c56"))

rownames(Hives_counts_vst_species_filtered) <- generate_species_v2$species
rownames(annotationrows_species) <- generate_species_v2$species

ordered_counts_species <- order(rowMeans(counts(Hives_dds_RLE_species_filtered,normalized=TRUE)), decreasing=TRUE)[1:30]
final_heatmap_species <- pheatmap(assay(Hives_counts_vst_species_filtered)[ordered_counts_species,],clustering_distance_cols = "correlation", cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = TRUE, annotation_col = select(Hives_annotation,Season), labels_col = c("DirectSM H4","DirectSM H5","DirectSM H6" ,"DirectSM H7","SM H4", "SM H5","SM H6","SM H7"), annotation_row = annotationrows_species, annotation_colors = ann_colors[1:2], height = 100, fontsize = 20)

ggsave(filename="top30_heatmap_species.pdf",device="pdf", plot=final_heatmap_species, height = 13, width = 15, path = output_path)


# draw volcano plots
get_volcano_plot <- function(dds_rle_object, contrast, annotation) {
  object <- DESeq(dds_rle_object)
  res_dds <- results(object, contrast=contrast)
  res <- lfcShrink(object,contrast = contrast, res=res_dds,type = 'normal')
  annotations <- annotation %>% rownames_to_column(var = "Taxonomic_ID")
  annotations$Taxonomic_ID <- as.numeric(annotations$Taxonomic_ID)
  EnhancedVolcano(res, lab = NA, x = 'log2FoldChange',y= 'pvalue', FCcutoff = 1,
                  colAlpha = 2,labSize = 2, xlim = c(-2,2))
}

# These MA plots are made using the methods as contrasts to observe effects due to method
family_volcano_plot <- get_volcano_plot(Hives_dds_RLE_family_filtered, contrast = c('Method','Direct_SM','SM'), annotation = Taxid_taxonomy_family)
genus_volcano_plot <- get_volcano_plot(Hives_dds_RLE_genus, contrast = c('Method','Direct_SM','SM'), annotation = Taxid_taxonomy_genus)
species_volcano_plot <- get_volcano_plot(Hives_dds_RLE_species, contrast = c('Method','Direct_SM','SM'), annotation = Taxid_taxonomy_species)


# This function will export the Taxnomoic IDs that have a log2foldchange > 1 or < -1. For this we use abs() - absolute value.
# We transform the Taxonomic IDs to faminly names etc.
# This will create a csv file
extract_lfc_names <- function(volcano.plot,outputname) {
  data_plot <- family_volcano_plot$data
  taxids_data <- rownames(data_plot)
  df_lfc_subset <- cbind(taxids_data,data_plot)
  colnames(df_lfc_subset)[1] <- "Taxonomic_ID"
  taxids <- as.character(df_lfc_subset$Taxonomic_ID)
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum",'nameNode.sqlite'))
  Superkingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "superkingdom",'nameNode.sqlite'))
  Kingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "kingdom", 'nameNode.sqlite'))
  Class <- as.character(getTaxonomy(taxids, desiredTaxa = "class", 'nameNode.sqlite'))
  Order <- as.character(getTaxonomy(taxids, desiredTaxa = "order", 'nameNode.sqlite'))
  Family <- as.character(getTaxonomy(taxids,desiredTaxa = "family", 'nameNode.sqlite'))
  Genus <- as.character(getTaxonomy(taxids,desiredTaxa = "genus", 'nameNode.sqlite'))
  Species <- as.character(getTaxonomy(taxids,desiredTaxa = "species", sqlFile = 'nameNode.sqlite'))
  df_lfc_export <- cbind(df_lfc_subset, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, Species)
  write.csv(df_lfc_export, file = outputname)
}

extract_lfc_names(family_volcano_plot,outputname = "./Output_data/Figure_3_output/family_lfc.csv")
extract_lfc_names(genus_volcano_plot,outputname = "./Output_data/Figure_3_output/genus_lfc.csv")
extract_lfc_names(species_volcano_plot,outputname = "./Output_data/Figure_3_output/species_lfc.csv")

export_significant_lfc <- function(volcano.plot,outputname) {
  data_plot <- volcano.plot$data
  taxids_data <- rownames(data_plot)
  df_lfc_subset <- cbind(taxids_data,data_plot)
  df_lfc_subset <- subset(df_lfc_subset, abs(df_lfc_subset$log2FoldChange) > 1)
  colnames(df_lfc_subset)[1] <- "Taxonomic_ID"
  taxids <- as.character(df_lfc_subset$Taxonomic_ID)
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum",'nameNode.sqlite'))
  Superkingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "superkingdom",'nameNode.sqlite'))
  Kingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "kingdom", 'nameNode.sqlite'))
  Class <- as.character(getTaxonomy(taxids, desiredTaxa = "class", 'nameNode.sqlite'))
  Order <- as.character(getTaxonomy(taxids, desiredTaxa = "order", 'nameNode.sqlite'))
  Family <- as.character(getTaxonomy(taxids,desiredTaxa = "family", 'nameNode.sqlite'))
  Genus <- as.character(getTaxonomy(taxids,desiredTaxa = "genus", 'nameNode.sqlite'))
  Species <- as.character(getTaxonomy(taxids,desiredTaxa = "species", sqlFile = 'nameNode.sqlite'))
  df_lfc_export <- cbind(df_lfc_subset, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, Species)
  write.csv(df_lfc_export, file = outputname)

}

export_significant_lfc(family_volcano_plot,outputname = "./Figures/Figure_3/Significant/Significant_families_volcano_plot.csv")
export_significant_lfc(genus_volcano_plot,outputname = "./Figures/Figure_3/Significant/Significant_genus_volcano_plot.csv")
export_significant_lfc(species_volcano_plot,outputname = "./Figures/Figure_3/Significant/Significant_species_volcano_plot.csv")


# Save the MA plots
ggsave(plot = family_volcano_plot, filename = "volcano_plot_family.pdf", path = output_path)
ggsave(plot = genus_volcano_plot, filename = "volcano_plot_genus.pdf", path = output_path)
ggsave(plot = species_volcano_plot, filename = "volcano_plot_species.pdf", path = output_path)



# Function to draw the significantly differentially abundant species/genera/families/whatever you throw at it
# The full model is the Method+Season and we compare it against a reduced model that includes only the Method
# In this way we can find out which organisms vary according to season.

plotDiffAbund <- function(colNums, DESeq_RLE_object, title, level = c("species","genus","family","order","phylum"), pathcsv) {
  
  dds_object <- DESeq(DESeq_RLE_object, test = "LRT", reduced = ~Method, full = ~Method+Season)
  rld <- rlog(dds_object, blind=F)
  results <- subset(results(dds_object), padj < 0.05)
  
  #results_export <- as.data.frame(rownames(results))
  #colnames(results_export) <- c("Taxonomic_ID")
  
  # make the lists
  upgenes <- rownames(head(results[ order( results$log2FoldChange ), ], n=1000))
  downgenes <- rownames(head(results[ order( -results$log2FoldChange ), ], n=1000))
  
  # this gives us the rows we want
  rows <- match(downgenes, row.names(rld))
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
  write.csv(annotationrows_species1,pathcsv)
  
  annotationrows_species1 = annotationrows_species1[,-2]
  annotationrows_species1 <- tibble::column_to_rownames(annotationrows_species1, var = level)
  rownames(mat) <- generate_species[,1]
  
  
  
  # the labels are hard coded at the moment :(
  df <- as.data.frame(colData(rld)[c("Season")])
  pheatmap(mat, fontsize=15, annotation_colors = ann_colors, annotation_row = annotationrows_species1, cluster_rows = FALSE, show_rownames = TRUE, height = 100, angle_col = "270", cluster_cols = TRUE, show_colnames = TRUE, border_color = NA, annotation_col = df, labels_col = c("DirectSM H4","DirectSM H5","DirectSM H6" ,"DirectSM H7","SM H4", "SM H5","SM H6","SM H7"), scale = 'row')
}

# Draw the heatmaps
differentially_abundant_species <- plotDiffAbund(
  DESeq_RLE_object =Hives_dds_RLE_species,
  colNums = c(1:8),
  level = "species", pathcsv = "./Figures/Figure_3/significant_species.csv")

differentially_abundant_genera <- plotDiffAbund(
  DESeq_RLE_object =Hives_dds_RLE_genus,
  colNums = c(1:8),
  level = "genus", pathcsv = "./Figures/Figure_3/significant_genera.csv")

differentially_abundant_families <- plotDiffAbund(
  DESeq_RLE_object =Hives_dds_RLE_family,
  colNums = c(1:8),
  level = "family",
  pathcsv = "./Figures/Figure_3/significant_families.csv")



ggsave(plot= differentially_abundant_species,filename = "./Figures/Figure_3/significantly_abundant_species.pdf",device="pdf", height = 10, width = 14)
ggsave(plot= differentially_abundant_genera,filename = "./Figures/Figure_3/significantly_abundant_genera.pdf",device="pdf", height = 10, width = 14)
ggsave(plot= differentially_abundant_families,filename = "./Figures/Figure_3/significantly_abundant_families.pdf",device="pdf", height = 10, width = 14)

# Now we can cluster the significantly abundant families/genera/species and observe patterns 
sig_res_LRT <- function(dds_object1, meta, replacegenes) {
  dds_object2 <- DESeq(dds_object1, test = "LRT", reduced = ~Method, full = ~Method+Season)
  res_LRT <- results(dds_object2)
  sig <- res_LRT %>% data.frame() %>% rownames_to_column(var="Taxonomic_ID") %>% as_tibble() %>%  filter(padj < 0.05)  
  clustering <- sig %>% arrange(padj) %>% head(n=1000)
  rld <- rlog(dds_object2, blind = F)
  rld_mat <- assay(rld)
  cluster_rlog <- rld_mat[clustering$Taxonomic_ID, ]
  meta$Method <- c("DirectSM","DirectSM","DirectSM","DirectSM","SM","SM","SM","SM")
  results_deg <- degPatterns(cluster_rlog, metadata = meta, time = "Season", col="Method", minc = 2)
  primary_plot <- degPlotCluster(results_deg$normalized, time = "Season", color = "Method", facet = T) + theme_bw() + scale_x_discrete(limits = c("May","July","November")) + scale_colour_brewer(type = "qual", palette = "Set1") + labs(title = "", y = "Z-score of abundance")
  primary_plot$data$title <- str_replace(primary_plot$data$title, "genes", print(replacegenes));
  return(primary_plot)
  }

hives_metadata_new <- hives_metadata %>% column_to_rownames(var = "X")

family_sig_res <- sig_res_LRT(Hives_dds_RLE_family, meta = hives_metadata_new, replacegenes = "Families")
genus_sig_res <- sig_res_LRT(Hives_dds_RLE_genus, meta = hives_metadata_new, replacegenes = "Genera")
species_sig_res <- sig_res_LRT(Hives_dds_RLE_species, meta = hives_metadata_new, replacegenes = "Species")

ggsave(plot = family_sig_res, filename = "./Figures/Figure_3/DEGpattern_family.pdf", height = 4)
ggsave(plot = genus_sig_res, filename = "./Figures/Figure_3/DEGpattern_genus.pdf", height = 4)
ggsave(plot = species_sig_res, filename = "./Figures/Figure_3/DEGpattern_species.pdf", height = 4)

# Export the Families/Genera/Species per cluster, and return taxonomy
# we use the split function to split the clusters in a list and then export them

# First we define a new classification function
cluster_taxonomy <- function(df2) {
  taxids <- df2$Taxonomic_ID
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum",'nameNode.sqlite'))
  Superkingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "superkingdom",'nameNode.sqlite'))
  Kingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "kingdom", 'nameNode.sqlite'))
  Class <- as.character(getTaxonomy(taxids, desiredTaxa = "class", 'nameNode.sqlite'))
  Order <- as.character(getTaxonomy(taxids, desiredTaxa = "order", 'nameNode.sqlite'))
  Family <- as.character(getTaxonomy(taxids,desiredTaxa = "family", 'nameNode.sqlite'))
  Genus <- as.character(getTaxonomy(taxids,desiredTaxa = "genus", 'nameNode.sqlite'))
  Species <- as.character(getTaxonomy(taxids,desiredTaxa = "species", sqlFile = 'nameNode.sqlite'))
  df_export <- cbind(df2, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, Species)
}


# New export function
export_list <- function(lst, FamGenSp, output_path) {
  names(lst) %>% purrr::walk(~write.csv(lst[[.]], paste0(output_path, FamGenSp,"_filtered_cluster_", ., ".csv")))
}

# Taxonomy list
taxonomy_per_cluster <- function(sig_res_object) {
 df <- sig_res_object$data[,c(1,7)]
 df[,1] <- str_replace(df[,1], "X", "")
 colnames(df) <- c("Taxonomic_ID", "cluster")
 df_new <- distinct(df, df$Taxonomic_ID, .keep_all = TRUE)
 df_list <- split(df_new, as.factor(df_new$cluster))
 new_list <- lapply(df_list,cluster_taxonomy)
}

listfamily <- taxonomy_per_cluster(family_sig_res)
listgenus <- taxonomy_per_cluster(genus_sig_res)
listspecies <- taxonomy_per_cluster(species_sig_res)


export_list(listfamily, FamGenSp = "Family",output_path = "./Figures/Figure_3/")
export_list(listgenus, FamGenSp = "Genus",output_path = "./Figures/Figure_3/")
export_list(listspecies, FamGenSp = "Species",output_path = "./Figures/Figure_3/")


# PCA analysis

# Filter for low counts and contaminants (species level only)
species_ncbi_list <- read_excel("./Non Plants/NCBI.xlsx")
contaminants <- species_ncbi_list[species_ncbi_list$`Relation category` == "Human cross contamination",]
contaminants_vector <- as.vector(contaminants$Taxonomic_ID)


filtered_family <- rowSums( counts(Hives_dds_RLE_family) >= 5) >=2
Hives_dds_RLE_family_filtered <- Hives_dds_RLE_family[filtered_family,]

filtered_genus <- rowSums( counts(Hives_dds_RLE_genus) >= 5) >=2
Hives_dds_RLE_genus_filtered <- Hives_dds_RLE_genus[filtered_genus,]

filtering_contaminants <- !(rownames(Hives_dds_RLE_species) %in% contaminants_vector)
Hives_dds_RLE_species <- Hives_dds_RLE_species[filtering_contaminants,]
# Filter for at least 100 counts in 2 samples
filtered_species <- rowSums( counts(Hives_dds_RLE_species) >= 100) >=2
Hives_dds_RLE_species_filtered <- Hives_dds_RLE_species[filtered_species,]

# set rownames
rownames(mat) <- species_annot2$species

# Draw the plot
pca_biplots <- function(dds_rle_object, annotation, filename) {
  # Prepare data 
  dds_object <- DESeq(dds_rle_object, test = "LRT", reduced = ~Method, full = ~Method+Season)
  rld_object <- rlog(dds_object, blind=T)
  mat <- assay(rld_object)
  colnames(mat) <- str_replace(colnames(mat), pattern = "_", replacement = " ")
  # ... and annotations
  tmp_annot <- annotation
  tmp_annot$Taxonomic_ID <- str_replace_all(tmp_annot$Taxonomic_ID, pattern = " ", replacement = "")
  tmp_annot$Taxonomic_ID <- as.numeric(tmp_annot$Taxonomic_ID)

  tmp_annot2 <- tmp_annot[match(rownames(mat),tmp_annot$Taxonomic_ID),]
  annotation_biplot <- get_super_kingdom_or_plant(tmp_annot2$Taxonomic_ID)
  annotation_biplot <- as.data.frame(annotation_biplot)
  rownames(annotation_biplot) <- tmp_annot2$Taxonomic_ID
  
  #Set rownames for mat
  rownames(mat) <- tmp_annot2[,2]
  
  # PCA time!
  pca<-prcomp(t(mat))
  
  # Export PC comeponents 
  eig <- fviz_eig(pca)
  ggsave(plot = eig, filename = paste0(output_path,filename,"PC_components.pdf"))
  
  # Biplot
  biplot <- factoextra::fviz_pca_biplot(pca, repel = TRUE, select.var = list(contrib = 20),
      col.ind = dds_object$Season, pointsize =2, geom.var = c("text","point"), geom.ind = c("arrow","text"),
      fill.var = annotation_biplot$superkingdom, fill.ind = dds_object$Season, pointshape = 21) + scale_fill_manual(name = c("Domain"), breaks = c("Bacteria","Viridiplantae","Eukaryota","Viruses") ,labels = c("Bacteria","Viridiplantae","Eukaryota","Viruses"), values = c("#E6AB02","#1B9E77","#7570B3","#eb53a1")) + theme_minimal() + scale_color_manual(name=c("Season"), labels = c("May","July","November"), values = c("#5cb300","#e8df2a","#8f7c56"))

  # this will increase the line width
  biplot[["layers"]][[2]][["aes_params"]][["size"]] <- 1
  to_export <- biplot + labs(x = paste("PC1:", paste0(round(eig[["data"]][["eig"]][1], digits = 1),"%")), y = paste("PC2:", paste0(round(eig[["data"]][["eig"]][2], digits = 1),"%")))
  ggsave(to_export,filename = paste0(output_path,filename,"PCA_biplot.pdf"))
  }



# Prepare the annotation for family and genus
annotationrows_family_plot <- getTaxonomy(rownames(annotationrows_family),desiredTaxa = "family", 'nameNode.sqlite')
annotationrows_family_plot <- as.data.frame(annotationrows_family_plot) %>% rownames_to_column(var="Taxonomic_ID")

annotationrows_genus_plot <- getTaxonomy(rownames(annotationrows_genus),desiredTaxa = "genus", 'nameNode.sqlite')
annotationrows_genus_plot <- as.data.frame(annotationrows_genus_plot) %>% rownames_to_column(var="Taxonomic_ID")
annotationrows_genus_plot$genus <- as.character(annotationrows_genus_plot$genus)
# Raplce row 740 with Bacillus walking sticks because it is recognised as a duplicate 
# Remember Bacillus is also a bacterium!
annotationrows_genus_plot[740,2] <- "Bacillus (walking stick)"

# Create the plots
pca_biplots(Hives_dds_RLE_family_filtered, annotation = annotationrows_family_plot, filename = "family_filtered")
pca_biplots(Hives_dds_RLE_genus_filtered, annotation = annotationrows_genus_plot , filename = "genus_filtered")
pca_biplots(Hives_dds_RLE_species_filtered, annotation = generate_species_v2, filename = "species_filtered")
