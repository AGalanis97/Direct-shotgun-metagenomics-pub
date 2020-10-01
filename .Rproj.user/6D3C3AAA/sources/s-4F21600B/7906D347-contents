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

# Taxonomizr will return the taxonomy for each species. For this purpose we have already included the library in the repository.
# If needed, the library can be rebuilt by the user with the following command:
# prepareDatabase('accessionTaxa.sql')
# Note that this will require a lot of time to run and complete. However its use is encouraged if complete reproducibility is of importance.
if (!require('DESeq2')) install.packages('DESeq2'); library('DESeq2')
if (!require('pheatmap')) install.packages('pheatmap'): library('pheatmap')
if (!require('taxize')) install.packages('taxize'): library('taxize')


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

# Optionally add the relative abundance as well
kraken2_files = lapply(kraken2_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[1,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,1] * 100;
  return(df)
})

# We will add the classification information
classification_ranks <- function(df) {
  setwd("C:/Users/user/Desktop/taxdmp")
  taxids <- df[,4]
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum"))
  Kingdom <- as.data.frame(getTaxonomy(taxids, desiredTaxa = "kingdom"))
  Superkingdom <- as.data.frame(getTaxonomy(taxids, desiredTaxa = "superkingdom"))
  cbind(df, Phylum, Kingdom, Superkingdom)
}

kraken2_files = lapply(kraken2_files, classification_ranks)

# Filter out human and Drosophila contamination. For human contamination, the Chordata clade is removed to ensure that no false negatives appear.
kraken2_files_family <- kraken2_files %>% lapply(filter, Rank_code == "F") %>% lapply(filter, Name != "Drosophilidae") %>% lapply(filter, is.na(Phylum)|Phylum != "Chordata")
kraken2_files_genus <- kraken2_files %>% lapply(filter, Rank_code == "G") %>% lapply(filter, Name != "Drosophila") %>% lapply(filter, is.na(Phylum)|Phylum != "Chordata")
kraken2_files_species <- kraken2_files %>% lapply(filter, Rank_code == "S") %>% lapply(filter, Name != "Drosophila melanogaster") %>% lapply(filter, is.na(Phylum)|Phylum != "Chordata")

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

# We apply the RLE normalisation
Hives_dds_family <- DESeqDataSetFromMatrix(countData = Hives_comparison_family, colData = hives_metadata, design = ~Method + Hive, tidy = TRUE)
Hives_dds_RLE_family <- estimateSizeFactors(Hives_dds_family,type = "ratio")
Hives_normalised_counts_family <- counts(Hives_dds_RLE_family, normalized = TRUE)
Hives_counts_vst_family <- varianceStabilizingTransformation(Hives_dds_RLE_family, blind = FALSE)

# Create the heatmap
Hives_annotation <- as.data.frame(colData(Hives_dds_family))

get_super_kingdom_or_plant <- function(x) {
  setwd("C:/Users/user/Desktop/taxdmp")
  ifelse(getTaxonomy(x, desiredTaxa = "superkingdom") != "Eukaryota", getTaxonomy(x, desiredTaxa = "superkingdom"), ifelse(getTaxonomy(x, desiredTaxa = "kingdom") == "Viridiplantae", getTaxonomy(x, desiredTaxa = "kingdom"), getTaxonomy(x, desiredTaxa = "superkingdom")))
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

final_heatmap_family <- pheatmap(assay(Hives_counts_vst_family)[ordered_counts_family,],cluster_rows = FALSE, show_rownames = FALSE, clustering_distance_cols = "correlation", cluster_cols = TRUE, annotation_col = select(Hives_annotation,Season,Varroa,Population), annotation_row = annotationrows_family, annotation_colors = ann_colors[1:2], main = "Heatmap families", height = 20)
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
Hives_dds_RLE_species <- estimateSizeFactors(Hives_dds_species,type = "ratio")
Hives_normalised_counts_species <- counts(Hives_dds_RLE_species, normalized = TRUE)
Hives_counts_vst_species <- varianceStabilizingTransformation(Hives_dds_RLE_species, blind = FALSE)


Hives_annotation <- as.data.frame(colData(Hives_dds_species))

TaxonomicIDs_species <- as.numeric(rownames((Hives_counts_vst_species)))
Taxid_taxonomy_species <- as.data.frame(get_super_kingdom_or_plant(TaxonomicIDs_species)) 
Taxonomic_IDs_dataframe_species <- as.data.frame(TaxonomicIDs_species)
annotationrows_species <- bind_cols(Taxid_taxonomy_species,Taxonomic_IDs_dataframe_species,generate_species)
colnames(annotationrows_species) <- c("Domain", "Taxonomic_ID","Species")

setwd("C:/Users/user/Desktop/taxdmp")
generate_species <- getTaxonomy(TaxonomicIDs_species,desiredTaxa = "species")
generate_species = as.data.frame(generate_species)
generate_species$species = as.character(generate_species$species)
generate_species[902,1] <- print("unknown1")
generate_species[1455,1] <- print("special1")
generate_species[1790,1] <- print("trial3")
generate_species[1792,1] <- print("ever2")
generate_species[1791,1] <- print("mask1")
generate_species[521,1] <- print("another1")

annotationrows_species = annotationrows_species[,-2]
annotationrows_species <- tibble::column_to_rownames(annotationrows_species, var = "Species")


ann_colors = list(
  Domain = c(Bacteria = "#0D0887FF", Viruses = "#F0F921FF", Eukaryota = "#CC4678FF", Viridiplantae = "#73D055FF", Archaea = "#ad05f5"),
  Season = c(May = "#5cb300", July = "#ecf547", November = "#8f7c56"))

rownames(Hives_counts_vst_species) <- generate_species$species
rownames(annotationrows_species) <- generate_species$species

ordered_counts_species <- order(rowMeans(counts(Hives_dds_RLE_species,normalized=TRUE)), decreasing=TRUE)[1:30]
final_heatmap_species <- pheatmap(assay(Hives_counts_vst_species)[ordered_counts_species,],clustering_distance_cols = "correlation", cluster_rows = FALSE, show_rownames = TRUE, cluster_cols = TRUE, annotation_col = select(Hives_annotation,Season), labels_col = c("DirectSM H4","DirectSM H5","DirectSM H6" ,"DirectSM H7","SM H4", "SM H5","SM H6","SM H7"), annotation_row = annotationrows_species, annotation_colors = ann_colors[1:2], main = "Heatmap top 30 species", height = 100, fontsize = 20)

ggsave(filename="top30_heatmap_species.pdf",device="pdf", plot=final_heatmap_species, height = 13, width = 15, path = "./Figures/Figure_3")

# Significantly abundant species
dds_object <- DESeq(Hives_dds_RLE_species2)
res <- results(dds_object)
resultsNames(dds_object)

resLFC <- lfcShrink(dds_object, coef="Season_November_vs_July", type="apeglm")


plotMA(resLFC)

sum(res$padj < 0.05, na.rm=TRUE)

resOrdered <- res[order(res$pvalue),]

resSig <- subset(resOrdered, padj < 0.05)

# Classification function for the DESeq2 object
classification_deseq <- function(df) {
  setwd("C:/Users/user/Desktop/taxdmp")
  taxids <- df$Taxonomic_ID
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum"))
  Superkingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "superkingdom"))
  Kingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "kingdom"))
  Class <- as.character(getTaxonomy(taxids, desiredTaxa = "class"))
  Order <- as.character(getTaxonomy(taxids, desiredTaxa = "order"))
  Family <- as.character(getTaxonomy(taxids,desiredTaxa = "family"))
  Genus <- as.character(getTaxonomy(taxids,desiredTaxa = "genus"))
  Species <- as.character(getTaxonomy(taxids,desiredTaxa = "species"))
  cbind(df, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, Species)
}

# NOT NEEDED
taxids112 <- data.frame(Number=is.numeric)
taxids12 <- rownames(resSig) %>% as.data.frame() %>% dplyr::rename(Taxonomic_ID = ".")
taxids_with_tax <- classification_deseq(taxids12)
taxizetry <- ncbi_get_taxon_summary(id = taxids12$Taxonomic_ID)


# Plot differentially abundant species
filtering_sign <- rownames(resSig)
annotation_sign <- annotationrows_species %>% tibble::rownames_to_column(var = "Taxonomic_ID")
annotation_sign2 <- annotation_sign %>% filter(., Taxonomic_ID %in% filtering_sign) %>% column_to_rownames(var = "Taxonomic_ID")

plotUpDownSigGenes <- function(results, colNums, rld, title) {
  
  # make the lists
  upgenes <- rownames(head(resSig[ order( resSig$log2FoldChange ), ], n=40))
  downgenes <- rownames(head(resSig[ order( -resSig$log2FoldChange ), ], n=40))
  
  # this gives us the rows we want
  rows <- match(upgenes, row.names(rld))
  mat <- assay(rld)[rows,c(1:8)]
  mat <- mat - rowMeans(mat)
  
  TaxonomicIDs_species <- as.numeric(rownames((mat)))
  Taxid_taxonomy_species <- as.data.frame(get_super_kingdom_or_plant(TaxonomicIDs_species)) 
  Taxonomic_IDs_dataframe_species <- as.data.frame(TaxonomicIDs_species)
  
  setwd("C:/Users/user/Desktop/taxdmp")
  generate_species <- getTaxonomy(TaxonomicIDs_species,desiredTaxa = "species")
  generate_species = as.data.frame(generate_species)
  generate_species$species = as.character(generate_species$species)
  annotationrows_species1 <- bind_cols(Taxid_taxonomy_species,Taxonomic_IDs_dataframe_species,generate_species)
  colnames(annotationrows_species) <- c("Domain", "Taxonomic_ID","Species")
  
  annotationrows_species1 = annotationrows_species1[,-2]
  annotationrows_species1 <- tibble::column_to_rownames(annotationrows_species1, var = "Species")
  rownames(mat) <- generate_species$species
  
  
  
  # the labels are hard coded at the moment :(
  df <- as.data.frame(colData(rld)[c("Season")])
  pheatmap(mat, fontsize=15, annotation_colors = ann_colors, annotation_row = annotationrows_species1, cluster_rows = FALSE, show_rownames = TRUE, main = "Significantly differentially abundant species", height = 100, angle_col = "270", cluster_cols = TRUE, show_colnames = TRUE, border_color = NA, annotation_col = df, labels_col = c("DirectSM H4","DirectSM H5","DirectSM H6" ,"DirectSM H7","SM H4", "SM H5","SM H6","SM H7"))
}

plotUpDownSigGenes(
  resSig,
  colNums = c(1:8),
  rld)

contrastDEGenes <- subset(results(dds_try, contrast=c("Season","May","July")), padj < 0.05)
aCols <- c(2,3,4)


# get the log transforms again
rld <- rlog(dds_try, blind=F)

# call
sig_plot <- plotUpDownSigGenes(
  resSig,
  colNums = c(1:8),
  rld)
ggsave(plot= sig_plot,filename = "significantly_abundant_species4.pdf",device="pdf", height = 9, width = 14)