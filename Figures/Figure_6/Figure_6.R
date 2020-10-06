# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Fig X seasonal fluctuation of plants
# This script is written in R version 3.6.3

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
if (!require('pheatmap')) install.packages('pheatmap'): library('pheatmap')

setwd(here::here())

# Calculate relative abundance
normalised_counts_methodseason <- read.csv("./Figures/Figure_6/Data_fig_6/normalised_counts_methodseason_species.csv")[,-1]
normalised_counts_methodseason <- normalised_counts_methodseason %>% round(.,digits = 0)

Relative_abundance_species <- normalised_counts_methodseason[,-1] %>% mutate_all(.,funs((./sum(.))*100))
rownames(Relative_abundance_species) <- normalised_counts_methodseason$Taxonomic_ID

# Add the taxonomy
classification_phyloseq <- function(df) {
  taxids <- rownames(df)
  Phylum <- as.character(getTaxonomy(taxids, desiredTaxa = "phylum", sqlFile = 'nameNode.sqlite'))
  Superkingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "superkingdom", sqlFile = 'nameNode.sqlite'))
  Kingdom <- as.character(getTaxonomy(taxids, desiredTaxa = "kingdom", sqlFile = 'nameNode.sqlite'))
  Class <- as.character(getTaxonomy(taxids, desiredTaxa = "class", sqlFile = 'nameNode.sqlite'))
  Order <- as.character(getTaxonomy(taxids, desiredTaxa = "order", sqlFile = 'nameNode.sqlite'))
  Family <- as.character(getTaxonomy(taxids,desiredTaxa = "family", sqlFile = 'nameNode.sqlite'))
  Genus <- as.character(getTaxonomy(taxids,desiredTaxa = "genus", sqlFile = 'nameNode.sqlite'))
  Species <- as.character(getTaxonomy(taxids,desiredTaxa = "species", sqlFile = 'nameNode.sqlite'))
  cbind(df, Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus, Species)
}

Relative_abundance_species_with_taxonomy <- classification_phyloseq(Relative_abundance_species)

# Extract the plants
Relative_abundance_plants <- Relative_abundance_species_with_taxonomy %>% tibble::rownames_to_column(var = "Taxonomic_ID") %>% filter(Phylum == "Streptophyta",) %>% tibble::column_to_rownames("Taxonomic_ID")
Relative_abundance_plants <- Relative_abundance_plants %>% mutate_if(., is.numeric,funs(ifelse(. <= 0.005, 0, .)))

# Get pairs per method
# First function removes the zeros across the samples 
remove_zero_across_samples <- function(df){
  df_x <- df %>% select_if(.,is.numeric)
  df_y <- df %>% select(.,"Species")
  df_z <- cbind(df_x,df_y) %>% remove_rownames(.) %>% dplyr::filter(rowSums(dplyr::select(., -Species)) > 0)
}

get_pair_method <- function(df,number=is.numeric){
  df_1 <- df %>% dplyr::select(.,contains(number)) %>% remove_rownames()
  df_2 <- df %>% dplyr::select(.,Species) %>% remove_rownames()
  df_3 <- cbind(df_1,df_2) %>% remove_zero_across_samples()
  return(df_3)
}

Hive4_plants <- Relative_abundance_plants %>% get_pair_method(number = "4") # 125 plants
Hive5_plants <- Relative_abundance_plants %>% get_pair_method(number = "5") # 173
Hive6_plants <- Relative_abundance_plants %>% get_pair_method(number = "6") # 107
Hive7_plants <- Relative_abundance_plants %>% get_pair_method(number = "7") # 67

# Now we look up some traits in the BIEN database. This will take a while, be patient !
Plants_BIEN <- function(df, Hive = is.numeric, method = c("DirectSM_H","SM_H")){
  x=paste(Hive)
  y=paste(method)
  plant_names <- df %>% dplyr::pull(.,Species) %>% as.character()
  flower_col <- plant_names %>% BIEN::BIEN_trait_traitbyspecies(species = ., trait = "flower color", all.taxonomy = FALSE, political.boundaries = FALSE, source.citation = FALSE) %>%  dplyr::select(c(scrubbed_species_binomial,trait_value))  
  flowering_begin <- plant_names %>% BIEN::BIEN_trait_traitbyspecies(species = ., trait = "plant flowering begin", all.taxonomy = FALSE, political.boundaries = FALSE, source.citation = FALSE) %>% dplyr::rename(flowering_begin_month = trait_value) %>% dplyr::filter(unit == "month") %>%  dplyr::select(c(scrubbed_species_binomial,flowering_begin_month))  
  flowering_duration <- plant_names %>% BIEN::BIEN_trait_traitbyspecies(species = ., trait = "plant flowering duration", all.taxonomy = FALSE, political.boundaries = FALSE, source.citation = FALSE) %>% dplyr::rename(flowering_duration_months = trait_value) %>% dplyr::select(c(scrubbed_species_binomial,flowering_duration_months))
  new_lst <- list(flower_col = flower_col, flowering_begin = flowering_begin,flowering_duration = flowering_duration)
  df_new <- new_lst %>% purrr::reduce(.,full_join, by = "scrubbed_species_binomial") %>% dplyr::rename(.,Species = scrubbed_species_binomial) %>% tibble::add_column(.,Method = factor(paste0(y,x)))
}

#Get plant data at species level
Hive4_DSM_plants <- Hive4_plants %>% Plants_BIEN(Hive = 4,method = "DirectSM_H")
Hive4_SM_plants <- Hive4_plants %>% Plants_BIEN(Hive = 4,method = "SM_H")
Hive5_DSM_plants <- Hive5_plants %>% Plants_BIEN(Hive = 5,method = "DirectSM_H")
Hive5_SM_plants <- Hive5_plants %>% Plants_BIEN(Hive = 5,method = "SM_H")
Hive6_DSM_plants <- Hive6_plants %>% Plants_BIEN(Hive = 6,method = "DirectSM_H")
Hive6_SM_plants <- Hive6_plants %>% Plants_BIEN(Hive = 6,method = "SM_H")
Hive7_DSM_plants <- Hive7_plants %>% Plants_BIEN(Hive = 7,method = "DirectSM_H")
Hive7_SM_plants <- Hive7_plants %>% Plants_BIEN(Hive = 7,method = "SM_H")


# Now we will look for occurences of these plants in the mediterranean. 
# First we will look for occurence data in the specified location. It works with either species or genus, just specify!
occurence_data <- function(df,level=c("species","genus")){
  if(level=="species")
    return({
      plant_names <- df %>% dplyr::pull(.,Species) %>% as.character()  
      plant_names %>% BIEN::BIEN_occurrence_box(species = .,min.lat = 33, max.lat =44, min.long = -9, max.long = 36.1, natives.only = FALSE, native.status = TRUE, cultivated = TRUE) %>% dplyr::select(c(scrubbed_species_binomial,country)) %>% dplyr::rename(Species = scrubbed_species_binomial)
    })
  if(level=="genus")
    return({
      plant_names <- df %>% dplyr::pull(.,Genus) %>% as.character()  
      plant_names %>% BIEN::BIEN_occurrence_box(genus = .,min.lat = 33, max.lat =44, min.long = -9, max.long = 36.1, natives.only = FALSE, native.status = TRUE, cultivated = TRUE) %>% dplyr::select(c(scrubbed_species_binomial,country)) %>% dplyr::rename(Species = scrubbed_species_binomial)
    })
} 

# Returns dataframe with true/false depending on species presence in the mediterranean basin
species_in_med <- function(df) {
  plant_names <- df %>% dplyr::pull(.,"Species") %>% as.character()
  occurrence <- df %>% occurence_data(level="species") %>% dplyr::select(.,Species) %>% unique() %>% dplyr::pull()
  df_1 <- plant_names %in% occurrence %>% as.data.frame() %>% dplyr::rename(.,Occurs_in_med = .) %>% cbind(plant_names,.) 
}

trial <- as.data.frame(as.character(Relative_abundance_plants$Species[Relative_abundance_plants$DirectSM_H4 != 0]))
colnames(trial) <- c("Species")
trial2 <- species_in_med(trial)

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
Hive4_med_DSM <- DSMH4_plants %>% species_in_med(.) %>% dplyr::group_by(Occurs_in_med) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H4")
Hive4_med_SM <- SMH4_plants %>% species_in_med(.) %>% dplyr::group_by(Occurs_in_med) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H4")
Hive5_med_DSM <- DSMH5_plants %>% species_in_med(.) %>% dplyr::group_by(Occurs_in_med) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H5")
Hive5_med_SM <- SMH5_plants %>% species_in_med(.) %>% dplyr::group_by(Occurs_in_med) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H5")
Hive6_med_DSM <- DSMH6_plants %>% species_in_med(.) %>% dplyr::group_by(Occurs_in_med) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H6")
Hive6_med_SM <- SMH6_plants %>% species_in_med(.) %>% dplyr::group_by(Occurs_in_med) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H6")
Hive7_med_DSM <- DSMH7_plants %>% species_in_med(.) %>% dplyr::group_by(Occurs_in_med) %>% dplyr::count() %>% tibble::add_column(.,Method = "DirectSM_H7")
Hive7_med_SM <- SMH7_plants %>% species_in_med(.) %>% dplyr::group_by(Occurs_in_med) %>% dplyr::count() %>% tibble::add_column(.,Method = "SM_H7")

data_med_presence <- rbind(Hive4_med_DSM,Hive4_med_SM,Hive5_med_DSM,Hive5_med_SM,Hive6_med_DSM,Hive6_med_SM,Hive7_med_DSM,Hive7_med_SM)
data_med_presence$Method <- factor(data_med_presence$Method, levels = c("DirectSM_H4","SM_H4","DirectSM_H5","SM_H5","DirectSM_H6", "SM_H6","DirectSM_H7","SM_H7"))

x_axis_labels <- c("SM H4","DirectSM H4","DirectSM H5","SM H5","DirectSM H6","SM H6","DirectSM H7","SM H7")
plot_med_presence <- ggplot() + geom_bar(aes(y = n, x = Method, fill= Occurs_in_med), data = data_med_presence, stat='identity', width=0.3, position = "fill") + theme_bw() + labs(y = "Present in the Mediterranean (%)", x="",title = "Plant species in the Mediterranean") + theme(text = element_text(size=12),axis.text.x = element_text(angle=45, hjust=1, size = 10)) + scale_fill_manual("Plant occurs in the Mediterranean", values=c("#ff9966","#3377ff")) + scale_y_continuous(labels = c(0,25,50,75,100)) + scale_x_discrete(labels = x_axis_labels)
ggsave(plot = plot_med_presence, filename = "./Figures/Figure_6/plants_med_filtered_abund005percent.pdf", height = 5, width = 6)
