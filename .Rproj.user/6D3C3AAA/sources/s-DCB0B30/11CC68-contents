# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figures 3A
# This script is written in R version 3.6.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('taxonomizr')) install.packages('taxonomizr'); library('taxonomizr')
if (!require('here')) install.packages('here'): library('here')
# Taxonomizr will return the taxonomy for each species. However, this requires that a database is built locally (requires 60 GB of space).
# prepareDatabase('accessionTaxa.sql')
# This process will take over 3 hours on a regular laptop/PC. Othherwise, please consider dowloading the zipped file 
# from here: and simply unzip it in the cloned repository. Place it at the top level, honeyDSM-seq and not in the subfolders.
if (!require('DESeq2')) install.packages('DESeq2'); library('DESeq2')
if (!require('pheatmap')) install.packages('pheatmap'): library('pheatmap')
if (!require('taxize')) install.packages('taxize'): library('taxize')

# Import the normalised reads table
Hives_normalised_counts <- read.csv("./Figures/Figure_3/normalised_counts_species.csv")
colnames(Hives_normalised_counts)[1] <- c("Taxonomic_ID")

Rounded_counts <- Hives_normalised_counts %>% mutate_at(colnames(Hives_normalised_counts)[2:9], funs(round(.,digits = 0))) %>% column_to_rownames(var = "Taxonomic_ID")

classification_phyloseq <- function(df) {
  taxids <- rownames(df)
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum", 'accessionTaxa.sql'))
  Superkingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "superkingdom", 'accessionTaxa.sql'))
  Kingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "kingdom", 'accessionTaxa.sql'))
  Class <- as.character(getTaxonomy(taxids, desiredTaxa = "class", 'accessionTaxa.sql'))
  Order <- as.character(getTaxonomy(taxids, desiredTaxa = "order", 'accessionTaxa.sql'))
  Family <- as.character(getTaxonomy(taxids,desiredTaxa = "family", 'accessionTaxa.sql'))
  Genus <- as.character(getTaxonomy(taxids,desiredTaxa = "genus", 'accessionTaxa.sql'))
  Species <- as.character(getTaxonomy(taxids,desiredTaxa = "species", 'accessionTaxa.sql'))
  cbind(df, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, Species)
}


Rounded_counts_with_tax <- classification_phyloseq(Rounded_counts)


# Remove NAs from Superkingdom
Rounded_counts_with_tax <- Rounded_counts_with_tax %>% rownames_to_column(var="Taxonomic_ID") %>%  drop_na(Superkingdom)

# Calculate relative abundance
Relative_abundance_hives <- Rounded_counts_with_tax %>% mutate_at(colnames(Rounded_counts_with_tax)[2:9],funs((./sum(.))*100))



# We copy the file into a new one so we have a new varibale to adjust and alter without needing to complete everything 
# from the beginning
Relative_abundance_with_tax <- Relative_abundance_hives[,-1]

DSM_H4 <- Relative_abundance_with_tax[,c(1,9,10)] %>% filter(.,DirectSM_H4 != 0) 
DSM_H5 <- Relative_abundance_with_tax[,c(2,9,10)] %>% filter(.,DirectSM_H5 != 0)
DSM_H6 <- Relative_abundance_with_tax[,c(3,9,10)] %>% filter(.,DirectSM_H6 != 0)
DSM_H7 <- Relative_abundance_with_tax[,c(4,9,10)] %>% filter(.,DirectSM_H7 != 0)
SM_H4 <- Relative_abundance_with_tax[,c(5,9,10)] %>% filter(.,SM_H4 != 0)
SM_H5 <- Relative_abundance_with_tax[,c(6,9,10)] %>% filter(.,SM_H5 != 0)
SM_H6 <- Relative_abundance_with_tax[,c(7,9,10)] %>% filter(.,SM_H6 != 0)
SM_H7 <- Relative_abundance_with_tax[,c(8,9,10)] %>% filter(.,SM_H7 != 0)

# Add them to a list
new_lyst <- list(DSM_H4,DSM_H5,DSM_H6,DSM_H7,SM_H4,SM_H5,SM_H6,SM_H7)
names(new_lyst) <- c("DSM_H4","DSM_H5","DSM_H6","DSM_H7","SM_H4","SM_H5","SM_H6","SM_H7")

# Create a custom function to count the different domains
count_kingdoms <- function(df, name) {
  Bacteria <- sum(df$Superkingdom == "Bacteria")
  Viridiplantae <- sum(df$Kingdom == "Viridiplantae", na.rm = T)
  Eukaryota <- (sum(df$Superkingdom == "Eukaryota") - Viridiplantae)
  Viruses <- sum(df$Superkingdom == "Viruses")
  Archaea <- sum(df$Superkingdom == "Archaea")
  df_1 <- rbind(Bacteria,Viruses,Archaea,Viridiplantae,Eukaryota) %>% as.data.frame();
  return(df_1)
}


kingdoms_per_sample <- lapply(new_lyst,count_kingdoms) %>% do.call(cbind,.)
kingdoms_per_sample2 <- lapply(new_lyst,count_kingdoms) %>% do.call(rbind,.) %>% t() %>% as.data.frame()
rownames(kingdoms_per_sample2) <- NULL

colnames(kingdoms_per_sample) <- c("DirectSM_H4","DirectSM_H5","DirectSM_H6","DirectSM_H7","SM_H4","SM_H5","SM_H6","SM_H7")
kingdom_for_plot <- t(kingdoms_per_sample) %>% as.data.frame() %>% rownames_to_column(var = "Sample")

# Calculate percentage of each domain so we can add it to the plot later
new_kingdom_with_percent1 <- kingdom_for_plot
new_kingdom_with_percent2 = round((new_kingdom_with_percent1[,-1]/rowSums(new_kingdom_with_percent1[,-1]))*100, digits = 2)
new_kingdom_with_percent3 = new_kingdom_with_percent2[,-3]
names_samples <- as.data.frame(kingdom_for_plot[,1])
new_kingdom_with_percent4 = cbind(names_samples,new_kingdom_with_percent3)
colnames(new_kingdom_with_percent4) <- c("Sample","Bacteria","Viruses","Viridiplantae","Eukaryota")

kingdom_melt_percent5 <- reshape2::melt(new_kingdom_with_percent4, id = "Sample")
kingdom_melt_percent5$Sample <- factor(kingdom_melt_percent5$Sample, levels = c("SM_H4","SM_H5","SM_H6","SM_H7","DirectSM_H4","DirectSM_H5","DirectSM_H6","DirectSM_H7"))


# Reshapre the data so we can plot it
kingdom_melt <- reshape2::melt(kingdom_for_plot, id = "Sample") %>% filter(variable != "Archaea")
kingdom_melt$Sample <- factor(kingdom_melt$Sample, levels = c("SM_H4","SM_H5","SM_H6","SM_H7","DirectSM_H4","DirectSM_H5","DirectSM_H6","DirectSM_H7"))

# Set axis labels
x_axis_labels <- c("SM H4","SM H5","SM H6","SM H7","DirectSM H4","DirectSM H5","DirectSM H6","DirectSM H7")

# Plot the data. Annotate the parts that have more than 20 species.
barplot_species_absolute <- ggplot(kingdom_melt, aes(x=Sample, y= value, fill=variable)) + geom_bar(stat = "identity") + theme_bw() + scale_fill_manual(values=c("#0D0887FF","#F0F921FF","#73D055FF","#CC4678FF")) + labs(y = "Number of species", x="", fill="Domain") + theme(text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1, size = 10)) + scale_x_discrete(labels = x_axis_labels) + geom_text(aes(label=ifelse(value>20,paste0(kingdom_melt_percent5$value,"%"),"")), position=position_stack(vjust=0.5) , colour="white") 

ggsave(path = "./Figures/Figure_3" ,plot = barplot_species_absolute, filename = "barplot_species_absolute1.pdf", device = "pdf", height = 5, width = 9)
