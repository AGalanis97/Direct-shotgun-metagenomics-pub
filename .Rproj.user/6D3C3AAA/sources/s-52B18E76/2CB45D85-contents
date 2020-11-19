# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figures 3A
# This script is written in R version 3.6.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('here')) install.packages('here'): library('here')
if (!require('plotly')) install.packages('plotly'): library('plotly')

normalised_families <- read.csv("./Figures/Figure_3/normalised_methodseason_family.csv")[,-1]
normalised_families$Taxonomic_ID <- as.character(normalised_families$Taxonomic_ID)

normalised_families_plants <- normalised_families %>% filter(Phylum == "Streptophyta",) 
Relative_abundance_families_plants <- normalised_families_plants[,2:9] %>% mutate_all(.,funs((./sum(.))*100))
Relative_abundance_families_plants <- cbind(Relative_abundance_families_plants,normalised_families_plants$Taxonomic_ID,normalised_families_plants$Family)
colnames(Relative_abundance_families_plants)[c(9:10)] <- c("Taxonomic_ID","Family")
Relative_abundance_families_plants$Taxonomic_ID <- as.character(Relative_abundance_families_plants$Taxonomic_ID)



plants <- Relative_abundance_families_plants %>% group_by(Family) %>% summarise_if(is.numeric, funs(sum))
plants$Sum <- rowSums(plants[,c(-1)])
plants$Name <- ifelse(plants$Sum > 0.5, plants$Sum, 0)
plants = plants %>% filter(.,Name != 0)

p <- plot_ly(plants, labels = ~Family, values = ~Name, type = 'pie',textposition = 'outside',textinfo = 'label+percent') %>%
  layout(title = 'Most abundant (>0.5%) families',
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
htmlwidgets::saveWidget(as_widget(p), "index.html")

# Count families, genera, species
normalised_species <- read.csv("./Figures/Figure_3/normalised_methodseason_species.csv")[,-1]
normalised_species <- normalised_species %>% filter(., Phylum == "Streptophyta")
normalised_species$Family <- as.character(normalised_species$Family)

family_number <- length(unique(na.omit(normalised_species$Family)))
genus_number <- length(unique(na.omit(normalised_species$Genus)))
species_number <- length(unique(na.omit(normalised_species$Species)))
