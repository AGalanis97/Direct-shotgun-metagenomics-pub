# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figure_X with the percent of species found in DSM, SM or Both
# This script is written in R version 3.6.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('here')) install.packages('here'): library('here')
if (!require('ggplot2')) install.packages('ggplot2'): library('ggplot2')
if (!require('reshape2')) install.packages('reshape2'): library('reshape2')
setwd(here::here())

Counts_with_taxonomy <- read.csv("./Figures/Figure_4/Data_fig_4/normalised_counts_species.csv")

# Counts_with_taxonomy <- Counts_with_taxonomy %>% tibble::rownames_to_column(var = "Taxonomic_ID") %>% filter(Phylum == "Streptophyta",) %>% tibble::column_to_rownames("Taxonomic_ID")


# Relative_abundance_plants_with_Hive7_filtered <- Counts_with_taxonomy %>% tibble::rownames_to_column(var = "Taxonomic_ID") %>% filter(Phylum == "Streptophyta",) %>% tibble::column_to_rownames("Taxonomic_ID")

Hive4_comparison <- Counts_with_taxonomy %>% select(c("Taxonomic_ID","DirectSM_H4","SM_H4")) %>% mutate(Presence = case_when(DirectSM_H4 == 0 & SM_H4 != 0 ~  "Only SM", DirectSM_H4 != 0 & SM_H4 == 0 ~  "Only DSM",DirectSM_H4 != 0 & SM_H4 != 0 ~  "Both", DirectSM_H4 == 0 & SM_H4 == 0 ~ "Both")) %>% filter(!is.na(Presence)) %>% tibble::column_to_rownames(var="Taxonomic_ID") 
Hive5_comparison <- Counts_with_taxonomy %>% select(c("Taxonomic_ID","DirectSM_H5","SM_H5")) %>% mutate(Presence = case_when(DirectSM_H5 == 0 & SM_H5 != 0 ~  "Only SM", DirectSM_H5 != 0 & SM_H5 == 0 ~  "Only DSM",DirectSM_H5 != 0 & SM_H5 != 0 ~  "Both",DirectSM_H5 == 0 & SM_H5 == 0 ~ "Both")) %>% filter(!is.na(Presence)) %>% tibble::column_to_rownames(var="Taxonomic_ID") 
Hive6_comparison <- Counts_with_taxonomy %>% select(c("Taxonomic_ID","DirectSM_H6","SM_H6")) %>% mutate(Presence = case_when(DirectSM_H6 == 0 & SM_H6 != 0 ~  "Only SM", DirectSM_H6 != 0 & SM_H6 == 0 ~  "Only DSM",DirectSM_H6 != 0 & SM_H6 != 0 ~  "Both",DirectSM_H6 == 0 & SM_H6 == 0 ~ "Both")) %>% filter(!is.na(Presence)) %>% tibble::column_to_rownames(var="Taxonomic_ID") 
Hive7_comparison <- Counts_with_taxonomy %>% select(c("Taxonomic_ID","DirectSM_H7","SM_H7")) %>% mutate(Presence = case_when(DirectSM_H7 == 0 & SM_H7 != 0 ~  "Only SM", DirectSM_H7 != 0 & SM_H7 == 0 ~  "Only DSM",DirectSM_H7 != 0 & SM_H7 != 0 ~  "Both",DirectSM_H7 == 0 & SM_H7 == 0 ~ "Both")) %>% filter(!is.na(Presence)) %>% tibble::column_to_rownames(var="Taxonomic_ID") 


Hive4_barplot <- c("Hive 4",sum(Hive4_comparison$Presence == "Only DSM"),sum(Hive4_comparison$Presence == "Only SM"),sum(Hive4_comparison$Presence == "Both"))
Hive5_barplot <- c("Hive 5",sum(Hive5_comparison$Presence == "Only DSM"),sum(Hive5_comparison$Presence == "Only SM"),sum(Hive5_comparison$Presence == "Both"))
Hive6_barplot <- c("Hive 6",sum(Hive6_comparison$Presence == "Only DSM"),sum(Hive6_comparison$Presence == "Only SM"),sum(Hive6_comparison$Presence == "Both"))
Hive7_barplot <- c("Hive 7",sum(Hive7_comparison$Presence == "Only DSM"),sum(Hive7_comparison$Presence == "Only SM"),sum(Hive7_comparison$Presence == "Both"))

data_method_to_reshape <- rbind(Hive4_barplot,Hive5_barplot,Hive6_barplot,Hive7_barplot) %>% as.data.frame() %>% remove_rownames() %>% dplyr::rename(Hive = V1) %>% dplyr::rename(Only_DSM = V2) %>% dplyr::rename(Only_SM = V3) %>% dplyr::rename(Both = V4)


data_method_per_hive <- melt(data_method_to_reshape, id.vars = "Hive") %>% dplyr::rename(Method = variable)
data_method_per_hive$value = as.numeric(data_method_per_hive$value)

barplot_comparison_method_per_hive <- ggplot(data_method_per_hive, aes(x=Hive, y=value, fill=Method)) + geom_bar(position = 'fill',stat = "identity", width = 0.5) + theme_bw() + labs(x="",y="Identified (%)") + scale_fill_viridis_d(labels = c("Only DirectSM","Only SM","Both")) + scale_y_continuous(labels = c(0,25,50,75,100)) + theme(text = element_text(size = 15), aspect.ratio = 0.3, legend.title = element_text(size = 13),
                                                                                                                                                                                                                                                                                                                                              legend.text = element_text(size = 10), legend.position = 'top') + coord_flip() 

ggsave(path = "./Figures/Figure_4/",plot = barplot_comparison_method_per_hive, filename = "barplot_comparison_method_per_hive.pdf")
