# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.
# This script is written in R version 4.0.3

# This script will reproduce Figure 2 and Supplementary Figures 3 B, C, D, and E

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggsci')) install.packages('ggsci'); library('ggsci')
if (!require('here')) install.packages('here'); library('here')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('hydroGOF')) install.packages('hydroGOF'); library('hydroGOF') # we use the RMSE function
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')



# Read the necessary classification resultfiles specific to Viridiplantae at Species level. Due to the here package we will simply refer to relative paths.
# The Rproj file is located at the top-level folder of this repository, so all paths are relative to it - or actually simply sub-directories.
# In the read.csv function we make use of the file.path to determine where the files are located. Otherwise read.csv is trying to read from the top level directory.

# We set the path where the variables are stored into a new variable

data_path <- "./Figures/Figure_2/Data_fig_2"

plant_species <- list.files(path = data_path, recursive = TRUE, pattern = "species_plants.csv")

comparison_species_plants <- plant_species %>%
  setNames(nm = .) %>%
  map_df(~read.csv(file.path(data_path,.)), .id = "file name")

comparison_species_plants$Software <- factor(comparison_species_plants$Software, levels = c("Mock","CCMetagen","DIAMOND","kraken2","minimap2"))

# This creates the stacked barplot
stacked_barplot_species_plants <- ggplot() + geom_bar(aes(y = Observed_Abundance, x = Software, fill= Species), data = comparison_species_plants, stat='identity', width=0.3, position="fill") + theme(legend.position = 'bottom', legend.text = element_text(size = 10), plot.title = element_text(hjust = 0.5)) + scale_fill_rickandmorty() + theme_bw(base_size =15) + scale_y_continuous(labels = function(x) paste0(x*100, "%"))

final_barplot_species_plants <- print(stacked_barplot_species_plants + labs(title= "Comparison Viridiplantae Species level",y="Observed Abundance (%)", x = ""))

# Export image as PDF
ggsave(path = "./Figures/Figure_2/", filename="final_barplot_species_plants.pdf",plot = final_barplot_species_plants, width = 9, device = "pdf")

# Same code as above, but for Genus level
plant_genus <- list.files(path = data_path, recursive = TRUE, pattern = "genus_plants.csv")

comparison_genus_plants <- plant_genus %>%
  setNames(nm = .) %>%
  map_df(~read.csv(file.path(data_path,.)), .id = "file name")

comparison_genus_plants$Software <- factor(comparison_genus_plants$Software, levels = c("Mock","CCMetagen","DIAMOND","kraken2","MG-RAST", "minimap2"))

stacked_barplot_genus_plants <- ggplot() + geom_bar(aes(y = Observed_Abundance, x = Software, fill= Genus), data = comparison_genus_plants, stat='identity', width=0.3, position = "fill") + theme(legend.position = 'bottom', legend.text = element_text(size = 8)) + scale_fill_rickandmorty() + theme_bw(base_size =15) + scale_y_continuous(labels = function(x) paste0(x*100, "%"))

final_barplot_genus_plants <- print(stacked_barplot_genus_plants + labs(title= "Comparison Viridiplantae Genus level",y="Observed Abundance (%)", x = ""))

ggsave(filename="final_barplot_genus_plants.pdf",plot = final_barplot_genus_plants, path = "./Figures/Figure_2/", device = "pdf", width = 9)

# Other dataset
others_dataset <- list.files(path = data_path, recursive = TRUE, pattern = "others.csv")

comparison_others <- others_dataset %>%
  setNames(nm = .) %>%
  map_df(~read.csv(file.path(data_path,.)), .id = "file name")

comparison_others$Software <- factor(comparison_others$Software, levels = c("Mock","CCMetagen","DIAMOND","kraken2","minimap2"))

stacked_barplot_others <- ggplot() + geom_bar(aes(y = Observed_Abundance, x = Software, fill= Species), data = comparison_others, stat='identity', width=0.3, position="fill") + theme(legend.position = 'bottom', legend.text = element_text(size = 8)) + scale_fill_rickandmorty() + theme_bw(base_size = 15) + scale_y_continuous(labels = function(x) paste0(x*100, "%"))

final_barplot_others <- print(stacked_barplot_others + labs(title= "Comparison for other organisms (non-Viridiplantae)",y="Observed Abundance (%)", x = ""))

ggsave(filename="final_barplot_others.pdf",plot = final_barplot_others, path = "./Figures/Figure_2/", device = "pdf", width = 9, height = 5)

# The code below reproduces the Community 2D figure

# First, load the data
CCMetagen_community <- read.csv("./Figures/Figure_2/Data_fig_2/CCMetagen_community.csv")
DIAMOND_community <- read.csv("./Figures/Figure_2/Data_fig_2/DIAMOND_community.csv")
kraken2_community <- read.csv("./Figures/Figure_2/Data_fig_2/kraken2_community.csv")
minimap2_community <- read.csv("./Figures/Figure_2/Data_fig_2/minimap2_community.csv")

# Create a function to make a scatter plot
scatterplot <- function(df, x, y) {
  ggscatter(df, x="Observed_Abundance", y= "Expected_Abundance", color = "Taxonomic_Domain", size = 3, palette = c("#E6AB02","#8a8a5c","#7570B3","#1B9E77","#eb53a1")) + stat_cor(size=7, label.x = 4, label.y = 18) + labs(x="Observed Abundance (%)", y="Expected Abundance (%)", color = "Taxonomic Domain", title = "") + theme_bw(base_size = 15) + geom_abline(intercept = 0, slope =1) + theme(legend.text = element_text(size=10))
}

# Apply the function to make scatter plots per tool

plot1 <- scatterplot(CCMetagen_community) + xlim(0,20) + ylim(0,20)
plot2 <- scatterplot(DIAMOND_community) + xlim(0,20) + ylim(0,20)
plot3 <- scatterplot(kraken2_community) + xlim(0,20) + ylim(0,20)
plot4 <- scatterplot(minimap2_community) + xlim(0,20) + ylim(0,20)

# Combine the scatter plots into one and edit legend before exporting the final one
final_scatterplot_no_export <- ggarrange(plot1,plot2,plot3,plot4, labels = c("A","B","C","D"),ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
final_scatterplot <- annotate_figure(final_scatterplot_no_export,  top = text_grob("Performance of tools on the community dataset", face = "bold", size = 14), bottom = text_grob("A: CCMetagen, B: DIAMOND, C: kraken2, D: minimap2"))

ggsave(path = "./Figures/Figure_2/", filename = "scatter_plot_performance.pdf", plot= final_scatterplot, height=8, width = 10, device = "pdf")

# RMSE calculation


temp_species_plants <- list.files(path = data_path, pattern = "species_plants.csv")
species_plants_list <- temp_species_plants %>%  map(~read.csv(file.path(data_path,.)))

temp_genus_plants <- list.files(path = data_path, pattern = "genus_plants.csv")
genus_plants_list <- temp_genus_plants[-c(4,6)] %>%  map(~read.csv(file.path(data_path,.))) # -c(4,6) to remove MG-RAST and the mock

temp_others <- list.files(path = data_path, pattern = "others.csv")
others_list <- temp_others[1:4] %>%  map(~read.csv(file.path(data_path,.))) # 1:4 in order to exclude the mock

temp_community <- list.files(path = data_path, pattern = "community.csv")
community_list <- temp_community %>%  map(~read.csv(file.path(data_path,.)))


# Create a custom function to add names to the list elements. 
# This custom function will grab the software name from the first column and assign it as the name of the element in the list. 
# We do this in order to know the results for each software separately.

getnames <- function(df){ 
  soft <- df[2,1] 
  return(soft)
}

list_naming <- unlist(lapply(species_plants_list,getnames), use.names = FALSE)
names(species_plants_list) <- list_naming

list_naming_genus <- unlist(lapply(species_plants_list,getnames), use.names = FALSE)
names(genus_plants_list) <- list_naming_genus

names(others_list) <- list_naming

names(community_list) <- list_naming

# This function reports the RMSE for each software. 

custom_rmse <- function(df){ 
  exp <- df[,"Expected_Abundance"] 
  obs <- df[,"Observed_Abundance"]
  rmse_value <- rmse(exp,obs)
  return(rmse_value)
}

RMSE_plant_species <- as.data.frame(do.call(rbind, lapply(species_plants_list,custom_rmse)))

RMSE_plant_genus <- as.data.frame(do.call(rbind, lapply(genus_plants_list,custom_rmse)))

RMSE_others <- as.data.frame(do.call(rbind, lapply(others_list,custom_rmse)))

RMSE_community <- as.data.frame(do.call(rbind, lapply(community_list,custom_rmse)))

RMSE_total <- cbind(RMSE_plant_species,RMSE_plant_genus,RMSE_others,RMSE_community)

colnames(RMSE_total) <- c("Plant_Species","Plant_Genus","Others","Community")
names_for_RMSE <- c("Viridiplantae: \n Species level","Viridiplantae: \n Genus level", "Others", "Community")

names(RMSE_total) <- names_for_RMSE
rounded_RMSE_total <- round(RMSE_total,digits = 3)
png("./Figures/Figure_2/RMSE_total.png",height = 50*nrow(rounded_RMSE_total), width = 100*ncol(rounded_RMSE_total))
grid.table(rounded_RMSE_total)
dev.off()

