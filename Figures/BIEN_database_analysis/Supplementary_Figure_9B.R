# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Supplementary Figure 9 B
# This script is written in R version 4.0.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('taxonomizr')) install.packages('taxonomizr'); library('taxonomizr')
if (!require('here')) install.packages('here'): library('here')
# Taxonomizr will return the taxonomy for each species. However, this requires that a database is built locally (requires 60 GB of space).
# prepareDatabase('accessionTaxa.sql')
# This process will take over 3 hours on a regular laptop/PC.  
# from here: and simply unzip it in the cloned repository. Place it at the top level, honeyDSM-seq and not in the subfolders.
if (!require('BIEN')) install.packages('BIEN'); library('BIEN')
if (!require('ggplot2')) install.packages('ggplot2'): library('ggplot2')

setwd(here::here())

# Calculate relative abundance
normalised_counts_methodseason <- read.csv("./Normalised_reads/normalised_methodseason_species_filtered.csv")[,-1]

Relative_abundance_species <- normalised_counts_methodseason[,-1] %>% select(1:8) %>% mutate_all(funs((./sum(.))*100))
rownames(Relative_abundance_species) <- normalised_counts_methodseason$Taxonomic_ID


Relative_abundance_species_with_taxonomy <- cbind(Relative_abundance_species,normalised_counts_methodseason[,c(10:17)])

# Extract the plants
Relative_abundance_plants <- Relative_abundance_species_with_taxonomy %>% tibble::rownames_to_column(var = "Taxonomic_ID") %>% filter(Phylum == "Streptophyta",) %>% tibble::column_to_rownames("Taxonomic_ID")
Relative_abundance_plants <- Relative_abundance_plants %>% mutate_if(., is.numeric,funs(ifelse(. <= 0.005, 0, .)))



# Now we will look for occurences of these plants in the mediterranean. 
# First we will look for occurence data in the specified location. It works with either species or genus, just specify!
occurence_data <- function(df,level=c("species","genus"), minlat,maxlat,maxlong,minlong){
  if(level=="species")
    return({
      plant_names <- df %>% dplyr::pull(.,Species) %>% as.character()  
      plant_names %>% BIEN::BIEN_occurrence_box(species = .,min.lat = minlat, max.lat =maxlat, min.long = minlong, max.long = maxlong, natives.only = FALSE, native.status = TRUE, cultivated = TRUE) %>% dplyr::select(c(scrubbed_species_binomial,country)) %>% dplyr::rename(Species = scrubbed_species_binomial)
    })
  if(level=="genus")
    return({
      plant_names <- df %>% dplyr::pull(.,Genus) %>% as.character()  
      plant_names %>% BIEN::BIEN_occurrence_box(genus = .,min.lat = minlat, max.lat =maxlat, min.long = minlong, max.long = maxlong, natives.only = FALSE, native.status = TRUE, cultivated = TRUE) %>% dplyr::select(c(scrubbed_species_binomial,country)) %>% dplyr::rename(Species = scrubbed_species_binomial)
    })
} 

# Returns dataframe with true/false depending on species presence in the mediterranean basin
species_in_box <- function(df, mn.lat,mx.lat,mn.long,mx.long) {
  plant_names <- df %>% dplyr::pull(.,"Species") %>% as.character()
  occurrence <- df %>% occurence_data(level="species", minlat = mn.lat,maxlat = mx.lat, minlong = mn.long, maxlong = mx.long) %>% dplyr::select(.,Species) %>% unique() %>% dplyr::pull()
  df_1 <- plant_names %in% occurrence %>% as.data.frame() %>% dplyr::rename(.,Occurs_in_box = .) %>% cbind(plant_names,.) 
}



DSMH4_plants <- as.data.frame(as.character(Relative_abundance_plants$Species[Relative_abundance_plants$DirectSM_H4 != 0]))
colnames(DSMH4_plants) <- c("Species")
SMH4_plants <- as.data.frame(as.character(Relative_abundance_plants$Species[Relative_abundance_plants$SM_H4 != 0]))
colnames(SMH4_plants) <- c("Species")

DSMH5_plants <- as.data.frame(as.character(Relative_abundance_plants$Species[Relative_abundance_plants$DirectSM_H5 != 0]))
colnames(DSMH5_plants) <- c("Species")
SMH5_plants <- as.data.frame(as.character(Relative_abundance_plants$Species[Relative_abundance_plants$SM_H5 != 0]))
colnames(SMH5_plants) <- c("Species")
DSMH6_plants <- as.data.frame(as.character(Relative_abundance_plants$Species[Relative_abundance_plants$DirectSM_H6 != 0]))
colnames(DSMH6_plants) <- c("Species")
SMH6_plants <- as.data.frame(as.character(Relative_abundance_plants$Species[Relative_abundance_plants$SM_H6 != 0]))
colnames(SMH6_plants) <- c("Species")
DSMH7_plants <- as.data.frame(as.character(Relative_abundance_plants$Species[Relative_abundance_plants$DirectSM_H7 != 0]))
colnames(DSMH7_plants) <- c("Species")
SMH7_plants <- as.data.frame(as.character(Relative_abundance_plants$Species[Relative_abundance_plants$SM_H7 != 0]))
colnames(SMH7_plants) <- c("Species")


# This will take a long long time so grab some coffee <3 I'm dead serious it'll take about 15 mins if not more
# First we use the parameters for Mediterranean box min lat = 33, max lat = 44, min long = -9, max long = 36.1
Hive4_med_DSM <- DSMH4_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -9,mx.long = 36.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H4")
Hive4_med_SM <- SMH4_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -9,mx.long = 36.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H4")
Hive5_med_DSM <- DSMH5_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -9,mx.long = 36.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H5")
Hive5_med_SM <- SMH5_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -9,mx.long = 36.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H5")
Hive6_med_DSM <- DSMH6_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -9,mx.long = 36.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H6")
Hive6_med_SM <- SMH6_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -9,mx.long = 36.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H6")
Hive7_med_DSM <- DSMH7_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -9,mx.long = 36.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H7")
Hive7_med_SM <- SMH7_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -9,mx.long = 36.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H7")

data_med_presence <- rbind(Hive4_med_DSM,Hive4_med_SM,Hive5_med_DSM,Hive5_med_SM,Hive6_med_DSM,Hive6_med_SM,Hive7_med_DSM,Hive7_med_SM)
data_med_presence$Method <- factor(data_med_presence$Method, levels = c("DirectSM_H4","SM_H4","DirectSM_H5","SM_H5","DirectSM_H6", "SM_H6","DirectSM_H7","SM_H7"))

x_axis_labels <- c("SM H4","DirectSM H4","DirectSM H5","SM H5","DirectSM H6","SM H6","DirectSM H7","SM H7")
plot_med_presence <- ggplot() + geom_bar(aes(y = n, x = Method, fill= Occurs_in_box), data = data_med_presence, stat='identity', width=0.3, position = "fill") + theme_bw() + labs(y = "Present in the Mediterranean (%)", x="",title = "Plant species in the Mediterranean") + theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1, size = 10)) + scale_fill_manual("Plant occurs in the Mediterranean", values=c("#ff9966","#3377ff")) + scale_y_continuous(labels = c(0,25,50,75,100)) + scale_x_discrete(labels = x_axis_labels)
ggsave(plot = plot_med_presence, filename = "./Figures/Figure_6/plants_med_filtered_abund1percent.pdf", height = 5, width = 6)

# Asia plot
Hive4_asia_DSM <- DSMH4_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = 90,mx.long = 135.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H4")
Hive4_asia_SM <- SMH4_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = 90,mx.long = 135.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H4")
Hive5_asia_DSM <- DSMH5_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = 90,mx.long = 135.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H5")
Hive5_asia_SM <- SMH5_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = 90,mx.long = 135.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H5")
Hive6_asia_DSM <- DSMH6_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = 90,mx.long = 135.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H6")
Hive6_asia_SM <- SMH6_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = 90,mx.long = 135.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H6")
Hive7_asia_DSM <- DSMH7_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = 90,mx.long = 135.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H7")
Hive7_asia_SM <- SMH7_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = 90,mx.long = 135.1) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H7")

data_asia_presence <- rbind(Hive4_asia_DSM,Hive4_asia_SM,Hive5_asia_DSM,Hive5_asia_SM,Hive6_asia_DSM,Hive6_asia_SM,Hive7_asia_DSM,Hive7_asia_SM)
data_asia_presence$Method <- factor(data_asia_presence$Method, levels = c("DirectSM_H4","SM_H4","DirectSM_H5","SM_H5","DirectSM_H6", "SM_H6","DirectSM_H7","SM_H7"))

x_axis_labels <- c("SM H4","DirectSM H4","DirectSM H5","SM H5","DirectSM H6","SM H6","DirectSM H7","SM H7")
plot_asia_presence <- ggplot() + geom_bar(aes(y = n, x = Method, fill= Occurs_in_box), data = data_asia_presence, stat='identity', width=0.3, position = "fill") + theme_bw() + labs(y = "Present in Asia (%)", x="",title = "Plant species in Asia") + theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1, size = 10)) + scale_fill_manual("Plant presence", values=c("#ff9966","#3377ff")) + scale_y_continuous(labels = c(0,25,50,75,100)) + scale_x_discrete(labels = x_axis_labels)

# America plot
Hive4_america_DSM <- DSMH4_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -121,mx.long = -75.9) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H4")
Hive4_america_SM <- SMH4_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -121,mx.long = -75.9) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H4")
Hive5_america_DSM <- DSMH5_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -121,mx.long = -75.9) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H5")
Hive5_america_SM <- SMH5_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -121,mx.long = -75.9) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H5")
Hive6_america_DSM <- DSMH6_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -121,mx.long = -75.9) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H6")
Hive6_america_SM <- SMH6_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -121,mx.long = -75.9) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H6")
Hive7_america_DSM <- DSMH7_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -121,mx.long = -75.9) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H7")
Hive7_america_SM <- SMH7_plants %>% species_in_box(.,mn.lat = 33,mx.lat = 44,mn.long = -121,mx.long = -75.9) %>% dplyr::group_by(Occurs_in_box) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H7")

data_america_presence <- rbind(Hive4_america_DSM,Hive4_america_SM,Hive5_america_DSM,Hive5_america_SM,Hive6_america_DSM,Hive6_america_SM,Hive7_america_DSM,Hive7_america_SM)
data_america_presence$Method <- factor(data_america_presence$Method, levels = c("DirectSM_H4","SM_H4","DirectSM_H5","SM_H5","DirectSM_H6", "SM_H6","DirectSM_H7","SM_H7"))

data_med_presence$Location <- c("Mediterranean")
data_america_presence$Location <- c("North_America")
data_asia_presence$Location <- c("Asia")

data_presence_combined <- rbind(data_med_presence,data_asia_presence,data_america_presence)
data_presence_combined$Location <- factor(data_presence_combined$Location)


x_axis_labels <- c("SM H4","DirectSM H4","DirectSM H5","SM H5","DirectSM H6","SM H6","DirectSM H7","SM H7")
plot_presence_combined <- ggplot() + geom_bar(aes(y = n, x = Method, fill= Occurs_in_box), data = data_presence_combined, stat='identity', width=0.3, position = "fill") + theme_bw() + facet_grid(. ~Location) + labs(y = "Present (%)", x="",title = "") + theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1, size = 10)) + scale_fill_manual("Plant presence", values=c("#ff9966","#3377ff"), labels = c("False","True")) + scale_y_continuous(labels = c(0,25,50,75,100)) + scale_x_discrete(labels = x_axis_labels)

ggsave(plot = plot_presence_combined, filename = "./Figures/Figure_6/plant_presence_locations.pdf", height = 4)
