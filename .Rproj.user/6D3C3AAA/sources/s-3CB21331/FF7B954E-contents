# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figures X with the abundance of significant families across seasons
# This script is written in R version 3.6.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('taxonomizr')) install.packages('taxonomizr'); library('taxonomizr')
if (!require('here')) install.packages('here'): library('here')
if (!require('ggplot2')) install.packages('ggplot2'): library('ggplot2')
if (!require('reshape2')) install.packages('reshape2'): library('reshape2')
if (!require('viridis')) install.packages('viridis'): library('viridis')


# Import the significantly abundant families
significant_families <- as.data.frame(read.csv("./Figures/Figure_3/significant_families.csv")[,-1])
colnames(significant_families) <- c("Taxonomic_ID")
significant_families_vector <- as.character(significant_families$Taxonomic_ID)    

# Import the normalised read data
normalised_families <- read.csv("./Figures/Figure_3/normalised_methodseason_family.csv")[,-1]
normalised_families$Taxonomic_ID <- as.character(normalised_families$Taxonomic_ID)

normalised_families_plants <- normalised_families %>% filter(Phylum == "Streptophyta",) 
Relative_abundance_families_plants <- normalised_families_plants[,2:9] %>% mutate_all(.,funs((./sum(.))*100))
Relative_abundance_families_plants <- cbind(Relative_abundance_families_plants,normalised_families_plants$Taxonomic_ID,normalised_families_plants$Family)
colnames(Relative_abundance_families_plants)[c(9:10)] <- c("Taxonomic_ID","Family")
Relative_abundance_families_plants$Taxonomic_ID <- as.character(Relative_abundance_families_plants$Taxonomic_ID)


#Getting sum abundance for family
sum_abudance_sig <- Relative_abundance_families_plants %>% filter(.,Taxonomic_ID %in% significant_families_vector)

Sum_abundance_families <- sum_abudance_sig %>% dplyr::select(.,c("Family",(dplyr::select_if(.,is.numeric) %>% colnames()))) %>% dplyr::group_by(Family) %>%  dplyr::summarise_each(sum) %>% filter(rowSums(dplyr::select(., -Family)) > 0) %>% t(.) %>% as.data.frame(.) %>% janitor::row_to_names(.,row_number = 1) %>% tibble::rownames_to_column(var="Method") %>% reshape2::melt(., id = "Method")

Sum_abundance_families$value <- as.numeric(Sum_abundance_families$value)
Sum_abundance_families$Method <- str_replace(Sum_abundance_families$Method, pattern = "_", replacement = " ")
Sum_abundance_families$Method <- factor(Sum_abundance_families$Method, levels = c("DirectSM H5","SM H5","DirectSM H6","SM H6","DirectSM H7", "SM H7","DirectSM H4","SM H4"))
Sum_abundance_families$variable <- as.factor(Sum_abundance_families$variable)

x_axis_labels <- c("SM H4","SM H5","SM H6","SM H7","DirectSM H4","DirectSM H5","DirectSM H6","DirectSM H7")


Sum_abundance_families_plot = ggplot(Sum_abundance_families, aes(x = Method, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21, stroke = 0.5) +  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 15), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")  + scale_size_continuous(limits = c(0,100),range = c(0,15), breaks = c(1,10,50,75)) + scale_fill_viridis_d(guide = F) + scale_y_discrete(limits = rev(Sum_abundance_families$variable)) 

  
ggsave(path = "./Figures/Figure_3/" ,plot = Sum_abundance_families_plot, filename = "sum_abundance_families_sig.pdf", device = "pdf")

significant_genus <- as.data.frame(read.csv("./Figures/Figure_3/significant_genera.csv")[,-1])
colnames(significant_genus) <- c("Taxonomic_ID")
significant_genus_vector <- as.character(significant_genus$Taxonomic_ID)    

# Import the normalised read data
normalised_genus <- read.csv("./Figures/Figure_3/normalised_methodseason_genus.csv")[,-1]
normalised_genus$Taxonomic_ID <- as.character(normalised_genus$Taxonomic_ID)

normalised_genus_plants <- normalised_genus %>% filter(Phylum == "Streptophyta",) 
Relative_abundance_genus_plants <- normalised_genus_plants[,2:9] %>% mutate_all(.,funs((./sum(.))*100))
Relative_abundance_genus_plants <- cbind(Relative_abundance_genus_plants,normalised_genus_plants$Taxonomic_ID,normalised_genus_plants$Genus)
colnames(Relative_abundance_genus_plants)[c(9:10)] <- c("Taxonomic_ID","Genus")
Relative_abundance_genus_plants$Taxonomic_ID <- as.character(Relative_abundance_genus_plants$Taxonomic_ID)


#Getting sum abundance for family
sum_abudance_sig_genus <- Relative_abundance_genus_plants %>% filter(.,Taxonomic_ID %in% significant_genus_vector)

Sum_abundance_genus <- sum_abudance_sig_genus %>% dplyr::select(.,c("Genus",(dplyr::select_if(.,is.numeric) %>% colnames()))) %>% dplyr::group_by(Genus) %>%  dplyr::summarise_each(sum) %>% filter(rowSums(dplyr::select(., -Genus)) > 0) %>% t(.) %>% as.data.frame(.) %>% janitor::row_to_names(.,row_number = 1) %>% tibble::rownames_to_column(var="Method") %>% reshape2::melt(., id = "Method")

Sum_abundance_genus$value <- as.numeric(Sum_abundance_genus$value)
Sum_abundance_genus$Method <- str_replace(Sum_abundance_genus$Method, pattern = "_", replacement = " ")
Sum_abundance_genus$Method <- factor(Sum_abundance_genus$Method, levels = c("DirectSM H5","SM H5","DirectSM H6","SM H6","DirectSM H7", "SM H7","DirectSM H4","SM H4"))
Sum_abundance_genus$variable <- as.factor(Sum_abundance_genus$variable)

x_axis_labels <- c("SM H4","SM H5","SM H6","SM H7","DirectSM H4","DirectSM H5","DirectSM H6","DirectSM H7")


Sum_abundance_genus_plot = ggplot(Sum_abundance_genus, aes(x = Method, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21, stroke = 0.5) +  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 15), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")  + scale_size_continuous(limits = c(0,100),range = c(0,15), breaks = c(1,10,50,75)) + scale_fill_viridis_d(guide = F) + scale_y_discrete(limits = rev(Sum_abundance_genus$variable)) 


ggsave(path = "./Figures/Figure_3/" ,plot = Sum_abundance_genus_plot, filename = "sum_abundance_genus_sig.pdf", device = "pdf")
