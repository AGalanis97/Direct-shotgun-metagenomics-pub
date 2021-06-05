# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figure 4 with the abundance of significant and most abundant families across seasons
# This script is written in R version 4.0

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('taxonomizr')) install.packages('taxonomizr'); library('taxonomizr')
if (!require('here')) install.packages('here'): library('here')
if (!require('ggplot2')) install.packages('ggplot2'): library('ggplot2')
if (!require('reshape2')) install.packages('reshape2'): library('reshape2')
if (!require('viridis')) install.packages('viridis'): library('viridis')


# Import the normalised read data
normalised_genus <- read.csv("./Normalised_reads/normalised_methodseason_genus.csv")[,-1]
normalised_genus$Taxonomic_ID <- as.character(normalised_genus$Taxonomic_ID)

# Import the genera that were significantly differentially abundant
significant_genus <- as.data.frame(read.csv("./Figures/Figure_3/significant_genera.csv")[,-1])
significant_genus_vector <- as.character(significant_genus$Taxonomic_ID)    


#Getting sum abundance for genus
normalised_genus_plants <- normalised_genus %>% filter(Phylum == "Streptophyta",) 
Relative_abundance_genus_plants <- normalised_genus_plants[,2:9] %>% mutate_all(.,funs((./sum(.))*100))
Relative_abundance_genus_plants <- cbind(Relative_abundance_genus_plants,normalised_genus_plants$Taxonomic_ID,normalised_genus_plants$Genus)
colnames(Relative_abundance_genus_plants)[c(9:10)] <- c("Taxonomic_ID","Genus")
Relative_abundance_genus_plants$Taxonomic_ID <- as.character(Relative_abundance_genus_plants$Taxonomic_ID)

# Extract significant per sample so we can merge them together
sum_abudance_sig_genus <- Relative_abundance_genus_plants %>% filter(.,Taxonomic_ID %in% significant_genus_vector)

Hive4_significant <- sum_abudance_sig_genus[,c(1,5,10)]
Hive5_significant <- sum_abudance_sig_genus[,c(2,6,10)]
Hive6_significant <- sum_abudance_sig_genus[,c(3,7,10)]
Hive7_significant <- sum_abudance_sig_genus[,c(4,8,10)]

Hive4_plants_genus <- Relative_abundance_genus_plants[,c(1,5,10)]
Hive5_plants_genus <- Relative_abundance_genus_plants[,c(2,6,10)]
Hive6_plants_genus <- Relative_abundance_genus_plants[,c(3,7,10)]
Hive7_plants_genus <- Relative_abundance_genus_plants[,c(4,8,10)]

Hive4_plants_genus <- (Hive4_plants_genus[order(Hive4_plants_genus$DirectSM_H4, decreasing = T),])[1:10,]
Hive5_plants_genus <- (Hive5_plants_genus[order(Hive5_plants_genus$DirectSM_H5, decreasing = T),])[1:10,]
Hive6_plants_genus <- (Hive6_plants_genus[order(Hive6_plants_genus$DirectSM_H6, decreasing = T),])[1:10,]
Hive7_plants_genus <- (Hive7_plants_genus[order(Hive7_plants_genus$DirectSM_H7, decreasing = T),])[1:10,]

Hive4_total <- full_join(Hive4_plants_genus,Hive4_significant)
Hive5_total <- full_join(Hive5_plants_genus,Hive5_significant)
Hive6_total <- full_join(Hive6_plants_genus,Hive6_significant)
Hive7_total <- full_join(Hive7_plants_genus,Hive7_significant)

hive_list_genus <- list(Hive4_total,Hive5_total,Hive6_total,Hive7_total)
order_hive_plants_genus <- hive_list_genus %>% purrr::reduce(full_join,by="Genus")
order_hive_plants_genus <- order_hive_plants_genus %>% mutate_if(is.numeric,~replace(., is.na(.),0))

Sum_abundance_genus_top10 <- order_hive_plants_genus %>% dplyr::select(.,c("Genus",(dplyr::select_if(.,is.numeric) %>% colnames()))) %>% dplyr::group_by(Genus) %>%  dplyr::summarise_each(sum) %>% filter(rowSums(dplyr::select(., -Genus)) > 0) %>% t(.) %>% as.data.frame(.) %>% janitor::row_to_names(.,row_number = 1) %>% tibble::rownames_to_column(var="Method") %>% reshape2::melt(., id = "Method")

Sum_abundance_genus_top10$value <- as.numeric(Sum_abundance_genus_top10$value)
Sum_abundance_genus_top10$Method <- str_replace(Sum_abundance_genus_top10$Method, pattern = "_", replacement = " ")
Sum_abundance_genus_top10$Method <- factor(Sum_abundance_genus_top10$Method, levels = c("DirectSM H5","SM H5","DirectSM H6","SM H6","DirectSM H7", "SM H7","DirectSM H4","SM H4"))
Sum_abundance_genus_top10$variable <- as.factor(Sum_abundance_genus_top10$variable)

x_axis_labels <- c("SM H4","SM H5","SM H6","SM H7","DirectSM H4","DirectSM H5","DirectSM H6","DirectSM H7")

Sum_abundance_genus_plot_top10 = ggplot(Sum_abundance_genus_top10, aes(x = Method, y = variable)) + 
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21, stroke = 0.5) +  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 15), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right")  + scale_size_continuous(limits = c(0,100),range = c(0,15), breaks = c(1,10,50,75)) + scale_fill_viridis_d(guide = F) + scale_y_discrete(limits = rev(Sum_abundance_genus_top10$variable)) 

ggsave(path = "./Figures/Figure_4/" ,plot = Sum_abundance_genus_plot_top10, filename = "most_abundant_and_significant_plant_genus.pdf", device = "pdf", height = 7)