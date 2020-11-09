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
devtools::install_github("hrbrmstr/waffle")

if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('pheatmap')) install.packages('pheatmap'): library('pheatmap')
if (!require('EnhancedVolcano')) BiocManager::install('EnhancedVolcano'): library('EnhancedVolcano')
if (!require("DEGreport")) BiocManager::install("DEGreport"): library('DEGreport')
# if (!require("lasso2")) install.packages("lasso2"): library('lasso2')
# Requires R 4.0 now. Install from archive.
normalised_species <- as.data.frame(read.csv("./Figures/Figure_3/normalised_methodseason_species.csv")[,-1])

normalised_species_plants <- normalised_species %>% filter(Phylum == "Streptophyta",) 
Relative_abundance_species_plants <- normalised_species_plants[,2:9] %>% mutate_all(.,funs((./sum(.))*100))
Relative_abundance_species_plants <- cbind(Relative_abundance_species_plants,normalised_species_plants$Taxonomic_ID,normalised_species_plants$Species)
colnames(Relative_abundance_species_plants)[c(9:10)] <- c("Taxonomic_ID","Species")
Relative_abundance_species_plants$Taxonomic_ID <- as.character(Relative_abundance_species_plants$Taxonomic_ID)


data_dsmh4 <- Relative_abundance_species_plants[,c(1,10)]
data_dsmh5 <- Relative_abundance_species_plants[,c(2,10)]
data_dsmh6 <- Relative_abundance_species_plants[,c(3,10)]
data_dsmh7 <- Relative_abundance_species_plants[,c(4,10)]
data_smh4 <- Relative_abundance_species_plants[,c(5,10)]
data_smh5 <- Relative_abundance_species_plants[,c(6,10)]
data_smh6 <- Relative_abundance_species_plants[,c(7,10)]
data_smh7 <- Relative_abundance_species_plants[,c(8,10)]

# Create a list
plants_list <- list(data_dsmh4,data_dsmh5,data_dsmh6,data_dsmh7,data_smh4,data_smh5,data_smh6,data_smh7)

get_pie_chart <- function(df, outputpath) {
  colnames(df) <- c("Abundance","Species")
  df_1 <- df[order(df[,1], decreasing = T),]
  df_1$Species <- as.character(df_1$Species)
  df_1$Species[c(4:nrow(df_1))] <- "Other"
  df_1 <- df_1 %>% group_by(Species) %>% summarise(Counts = sum(Abundance)) 
  df_1$Counts <- round(df_1$Counts, digits = 0)
  
  df_2 <- df_1 %>% 
    mutate(end = 2 * pi * cumsum(Counts)/sum(Counts),
           start = lag(end, default = 0),
           middle = 0.5 * (start + end),
           hjust = ifelse(middle > pi, 1, 0),
           vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
  
  library(ggforce) # for 'geom_arc_bar'
  plot <- ggplot(df_2) + 
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                     start = start, end = end, fill = Species)) +
    geom_text(aes(x = 1.05 * sin(middle), y = 1.05 * cos(middle), label = paste0(df_1$Species," ", df_1$Counts,"%"),
                  hjust = hjust, vjust = vjust)) +
    coord_fixed() +
    scale_x_continuous(limits = c(-3, 3),  # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) +
    scale_y_continuous(limits = c(-1, 1),      # Adjust so labels are not cut off
                       name = "", breaks = NULL, labels = NULL) + theme_void() + scale_fill_brewer(type = "div","Set1") + theme(legend.position = "none")
  ggsave(plot=plot, filename=outputpath,device = "pdf", width = 10)
}

get_pie_chart(data_smh4, outputpath = "./Figures/Figure_9/SM_H4.pdf")
get_pie_chart(data_smh5, outputpath = "./Figures/Figure_9/SM_H5.pdf")
get_pie_chart(data_smh6, outputpath = "./Figures/Figure_9/SM_H6.pdf")
get_pie_chart(data_smh7, outputpath = "./Figures/Figure_9/SM_H7.pdf")
get_pie_chart(data_dsmh4, outputpath = "./Figures/Figure_9/DSM_H4.pdf")
get_pie_chart(data_dsmh5, outputpath = "./Figures/Figure_9/DSM_H5.pdf")
get_pie_chart(data_dsmh6, outputpath = "./Figures/Figure_9/DSM_H6.pdf")
get_pie_chart(data_dsmh7, outputpath = "./Figures/Figure_9/DSM_H7.pdf")

