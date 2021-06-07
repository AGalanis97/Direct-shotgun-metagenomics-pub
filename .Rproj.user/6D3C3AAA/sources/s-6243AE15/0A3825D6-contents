# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Supplementary Figures 4, 5, 6
# This script is written in R version 4.0.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('here')) install.packages('here'): library('here')
if (!require('hydroGOF')) install.packages('hydroGOF'); library('hydroGOF')

# Create RMSE evaluation functions
get_rmse_eval <- function(df) {
  df_1 <- data.frame()
  for (f in 2:((ncol(df))-1)) {
    x <- rmse_evaluation(df,df[ncol(df)],df[f])
    df_1 <- rbind(df_1,x)
  }
  return(df_1)
}

rmse_evaluation <- function(df, expected, observed){ 
  rmse_value <- rmse(expected,observed)
  return(rmse_value)
}

kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_A23 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A23",full.names = T, pattern = "\\.kraken2")

A23_mock <- lapply(import_A23, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_A23 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A23", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

A23_mock_files <- setNames(A23_mock, substring(naming_list_A23, first  = 1, last = nchar(naming_list_A23) -8))

A23_mock_files = lapply(A23_mock_files,setNames,kraken2_output_names)

A23_mock_files = lapply(A23_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
A23_mock_files = lapply(A23_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_A23 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_A23.csv") 
composition_A23$True_relative_abundance <- (composition_A23$True_reads/sum(composition_A23$True_reads)) * 100
composition_A23 <- composition_A23[,-c(2,3)]
A23_filtering <- composition_A23$Taxonomic_ID

A23_filtered <- A23_mock_files %>% lapply(filter, Taxonomic_ID %in% A23_filtering)

A23_filtered_2 <- lapply(A23_filtered, "[", c(4,6))

A23_mock_final <- lapply(names(A23_filtered_2), function(x){
  colnames(A23_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  A23_filtered_2[[x]]
})
names(A23_mock_final) <- names(A23_filtered_2) 


A23_confidence <- A23_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

A23_mock_with_confidence <- merge(A23_confidence,composition_A23, by="Taxonomic_ID")
A23_mock_with_confidence$True_reads <- as.numeric(A23_mock_with_confidence$True_reads)

rmse_A23 <- get_rmse_eval(A23_mock_with_confidence)




kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_B23 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B23",full.names = T, pattern = "\\.kraken2")

B23_mock <- lapply(import_B23, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_B23 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B23", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

B23_mock_files <- setNames(B23_mock, substring(naming_list_B23, first  = 1, last = nchar(naming_list_B23) -8))

B23_mock_files = lapply(B23_mock_files,setNames,kraken2_output_names)

B23_mock_files = lapply(B23_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
B23_mock_files = lapply(B23_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_B23 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_C23.csv") 
composition_B23$True_relative_abundance <- (composition_B23$True_reads/sum(composition_B23$True_reads)) * 100
composition_B23 <- composition_B23[,-c(2,3)]
B23_filtering <- composition_B23$Taxonomic_ID

B23_filtered <- B23_mock_files %>% lapply(filter, Taxonomic_ID %in% A23_filtering)

B23_filtered_2 <- lapply(B23_filtered, "[", c(4,6))

B23_mock_final <- lapply(names(B23_filtered_2), function(x){
  colnames(B23_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  B23_filtered_2[[x]]
})
names(B23_mock_final) <- names(B23_filtered_2) 


B23_confidence <- B23_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

B23_mock_with_confidence <- merge(B23_confidence,composition_B23, by="Taxonomic_ID")

rmse_B23 <- get_rmse_eval(B23_mock_with_confidence)


kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_C23 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C23",full.names = T, pattern = "\\.kraken2")

C23_mock <- lapply(import_C23, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_C23 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C23", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

C23_mock_files <- setNames(C23_mock, substring(naming_list_C23, first  = 1, last = nchar(naming_list_C23) -8))

C23_mock_files = lapply(C23_mock_files,setNames,kraken2_output_names)

C23_mock_files = lapply(C23_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
C23_mock_files = lapply(C23_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_C23 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_C23.csv") 
composition_C23$True_relative_abundance <- (composition_C23$True_reads/sum(composition_C23$True_reads)) * 100
composition_C23 <- composition_C23[,-c(2,3)]
C23_filtering <- composition_C23$Taxonomic_ID

C23_filtered <- C23_mock_files %>% lapply(filter, Taxonomic_ID %in% A23_filtering)

C23_filtered_2 <- lapply(C23_filtered, "[", c(4,6))

C23_mock_final <- lapply(names(C23_filtered_2), function(x){
  colnames(C23_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  C23_filtered_2[[x]]
})
names(C23_mock_final) <- names(C23_filtered_2) 


C23_confidence <- C23_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

C23_mock_with_confidence <- merge(C23_confidence,composition_C23, by="Taxonomic_ID")

rmse_C23 <- get_rmse_eval(C23_mock_with_confidence)


kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_A64 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A64",full.names = T, pattern = "\\.kraken2")

A64_mock <- lapply(import_A64, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_A64 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A64", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

A64_mock_files <- setNames(A64_mock, substring(naming_list_A64, first  = 1, last = nchar(naming_list_A64) -8))

A64_mock_files = lapply(A64_mock_files,setNames,kraken2_output_names)

A64_mock_files = lapply(A64_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
A64_mock_files = lapply(A64_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_A64 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_A64.csv") 
composition_A64$True_relative_abundance <- (composition_A64$True_reads/sum(composition_A64$True_reads)) * 100
composition_A64 <- composition_A64[,-c(2,3)]
A64_filtering <- composition_A64$Taxonomic_ID

A64_filtered <- A64_mock_files %>% lapply(filter, Taxonomic_ID %in% A64_filtering)

A64_filtered_2 <- lapply(A64_filtered, "[", c(4,6))

A64_mock_final <- lapply(names(A64_filtered_2), function(x){
  colnames(A64_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  A64_filtered_2[[x]]
})
names(A64_mock_final) <- names(A64_filtered_2) 


A64_confidence <- A64_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

A64_mock_with_confidence <- merge(A64_confidence,composition_A64, by="Taxonomic_ID")

rmse_A64 <- get_rmse_eval(A64_mock_with_confidence)



kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_B64 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B64",full.names = T, pattern = "\\.kraken2")

B64_mock <- lapply(import_B64, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_B64 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B64", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

B64_mock_files <- setNames(B64_mock, substring(naming_list_B64, first  = 1, last = nchar(naming_list_B64) -8))

B64_mock_files = lapply(B64_mock_files,setNames,kraken2_output_names)

B64_mock_files = lapply(B64_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
B64_mock_files = lapply(B64_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_B64 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_B64.csv") 
composition_B64$True_relative_abundance <- (composition_B64$True_reads/sum(composition_B64$True_reads)) * 100
composition_B64 <- composition_B64[,-c(2,3)]
B64_filtering <- composition_B64$Taxonomic_ID

B64_filtered <- B64_mock_files %>% lapply(filter, Taxonomic_ID %in% B64_filtering)

B64_filtered_2 <- lapply(B64_filtered, "[", c(4,6))

B64_mock_final <- lapply(names(B64_filtered_2), function(x){
  colnames(B64_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  B64_filtered_2[[x]]
})
names(B64_mock_final) <- names(B64_filtered_2) 


B64_confidence <- B64_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

B64_mock_with_confidence <- merge(B64_confidence,composition_B64, by="Taxonomic_ID")

rmse_B64 <- get_rmse_eval(B64_mock_with_confidence)



kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_C64 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C64",full.names = T, pattern = "\\.kraken2")

C64_mock <- lapply(import_C64, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_C64 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C64", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

C64_mock_files <- setNames(C64_mock, substring(naming_list_C64, first  = 1, last = nchar(naming_list_C64) -8))

C64_mock_files = lapply(C64_mock_files,setNames,kraken2_output_names)

C64_mock_files = lapply(C64_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
C64_mock_files = lapply(C64_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_C64 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_C64.csv") 
composition_C64$True_relative_abundance <- (composition_C64$True_reads/sum(composition_C64$True_reads)) * 100
composition_C64 <- composition_C64[,-c(2,3)]
C64_filtering <- composition_C64$Taxonomic_ID

C64_filtered <- C64_mock_files %>% lapply(filter, Taxonomic_ID %in% C64_filtering)

C64_filtered_2 <- lapply(C64_filtered, "[", c(4,6))

C64_mock_final <- lapply(names(C64_filtered_2), function(x){
  colnames(C64_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  C64_filtered_2[[x]]
})
names(C64_mock_final) <- names(C64_filtered_2) 


C64_confidence <- C64_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

C64_mock_with_confidence <- merge(C64_confidence,composition_C64, by="Taxonomic_ID")

rmse_C64 <- get_rmse_eval(C64_mock_with_confidence)

all_rmse_values <- cbind(rmse_A23,rmse_B23,rmse_C23,rmse_A64,rmse_B64,rmse_C64)
colnames(all_rmse_values) <- c("A23","B23","C23","A64","B64","C64")
all_rmse_values$Confidence <- c(seq(0,1, by=0.1))
for_plotting <- reshape2::melt(all_rmse_values, id = "Confidence")
for_plotting$Species_nr <- c(rep(c(23,64), each = 33))
for_plotting$Species_nr <- as.factor(for_plotting$Species_nr)
colnames(for_plotting) <- c("Confidence","Mock","RMSE","Species_nr")

boxplot2 <- ggboxplot(for_plotting, x = "Confidence", y = "RMSE", add = "jitter", color = "Species_nr", palette = get_palette("npg",2)) + ylim(0,7) + theme(text = element_text(size=12), legend.position = "top", axis.text.x = element_text(size=10)) + scale_x_discrete(labels = c(0,"",0.2,"","",0.5,"","",0.8,"",1)) + stat_compare_means(aes(group=Species_nr),method = "t.test", label = "p.signif")

boxplot3 <- ggpar(boxplot2, legend.title = "")

ggsave(plot = boxplot3, filename = "RMSE_all_new_mocks.pdf", width = 3,height = 3)


# Number of species
species_64 <- c(A64_mock_files,B64_mock_files,C64_mock_files)
species_23 <- c(A23_mock_files,B23_mock_files,C23_mock_files)

numberspecies64 <- species_64 %>% lapply(filter, Rank_code == "S") 
numberspecies641 <- as.data.frame(sapply(numberspecies64, nrow)) %>% tibble::rownames_to_column(var = "Sample")
numberspecies641$Mock <- rep(c("A64","B64","C64"),each=11)
numberspecies641$Confidence <- rep(seq(from=0.0,to=1,by=0.1))
numberspecies641_new <- numberspecies641 %>% rename(., Number_of_species = names(numberspecies641[2])) %>% select(., -Sample)
numberspecies641_new$Mock <- factor(numberspecies641$Mock)

numberspecies23 <- species_23 %>% lapply(filter, Rank_code == "S")
numberspecies231 <- as.data.frame(sapply(numberspecies23, nrow)) %>% tibble::rownames_to_column(var = "Sample")
numberspecies231$Mock <- rep(c("A23","B23","C23"),each=11)
numberspecies231$Confidence <- rep(seq(from=0.0,to=1,by=0.1))
numberspecies231_new <- numberspecies231 %>% rename(., Number_of_species = names(numberspecies231[2])) %>% select(., -Sample)
numberspecies231_new$Mock <- factor(numberspecies231$Mock)

all_mocks_species <- rbind(numberspecies641_new,numberspecies231_new)

all_mocks_species_plot <- ggline(all_mocks_species, "Confidence", "Number_of_species",
                                 linetype = "Mock", shape = "Mock",
                                 color = "Mock") + geom_hline(yintercept = 64, linetype = "dashed") + theme_bw() + theme(text = element_text(size=15), axis.text.x = element_text(size=12))
all_mocks_species_plot2 <- ggpar(all_mocks_species_plot, ylab = "Number of species identified")

ggsave(plot=all_mocks_species_plot2,filename="dataplotspecies_new_mocks.pdf", height = 6.5, width=4)


ggline(sens_plotting, x = "Confidence", y = "value", add = "mean_se",
       color = "Species_nr", palette = "npg")+
  stat_compare_means(label = "p.signif", ref.group = "0", method = "t.test")  + ylab("Sensitivity") + scale_x_discrete(labels=c(0,"",0.2,"","",0.5,"","",0.8,"",1))




# Sensitivity for 23 vs 64

counting_64_species <- function(df) {
  species_selection <- df %>% filter(., Rank_code == "S")
  non_zero_species <- species_selection %>% filter(., Sample_reads > 0)
  taxids <- non_zero_species %>% dplyr::pull(.,Taxonomic_ID) %>% as.character()
  df_1 <- taxids %in% A64_filtering %>% as.data.frame() %>% dplyr::rename(.,Presence = .) %>% cbind(taxids,.) 
}

counting_23_species <- function(df) {
  species_selection <- df %>% filter(., Rank_code == "S")
  non_zero_species <- species_selection %>% filter(., Sample_reads > 0)
  taxids <- non_zero_species %>% dplyr::pull(.,Taxonomic_ID) %>% as.character()
  df_1 <- taxids %in% A23_filtering %>% as.data.frame() %>% dplyr::rename(.,Presence = .) %>% cbind(taxids,.) 
}

count_true_false_species <- function(df) {
  vector1 <- df[,2];
  return(sum(vector1))
}

return_sensitivity_and_ppv_64 <- function(lyst) {
  A64_mock_files
  
  df_true_false <- lapply(lyst, counting_64_species)
  df_truefalse_count <- lapply(df_true_false, count_true_false_species) %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Origin") %>% rename(., True_species = V1)
  df_species_per_sample <- lapply(df_true_false, nrow) %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Origin") %>% rename(., Total_species = V1)
  
  df_truefalse_count_total <- cbind(df_truefalse_count,df_species_per_sample) %>% column_to_rownames(var = "Origin")
  df_truefalse_count_total <- df_truefalse_count_total[,-2]
  colnames(df_truefalse_count_total) <- c("Positive","Total")
  
  df_truefalse_count_total$Positivevalue = df_truefalse_count_total$Positive/df_truefalse_count_total$Total
  df_truefalse_count_total$Sensitivity = df_truefalse_count_total$Positive/(df_truefalse_count_total$Positive + (64 - df_truefalse_count_total$Positive));
  return(df_truefalse_count_total)
}

return_sensitivity_and_ppv_23 <- function(lyst) {
  df_true_false <- lapply(lyst, counting_23_species)
  df_truefalse_count <- lapply(df_true_false, count_true_false_species) %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Origin") %>% rename(., True_species = V1)
  df_species_per_sample <- lapply(df_true_false, nrow) %>% as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Origin") %>% rename(., Total_species = V1)
  
  df_truefalse_count_total <- cbind(df_truefalse_count,df_species_per_sample) %>% column_to_rownames(var = "Origin")
  df_truefalse_count_total <- df_truefalse_count_total[,-2]
  colnames(df_truefalse_count_total) <- c("Positive","Total")
  
  df_truefalse_count_total$Positivevalue = df_truefalse_count_total$Positive/df_truefalse_count_total$Total
  df_truefalse_count_total$Sensitivity = df_truefalse_count_total$Positive/(df_truefalse_count_total$Positive + (23 - df_truefalse_count_total$Positive));
  return(df_truefalse_count_total)
}


A64_sens_ppv <- return_sensitivity_and_ppv_64(A64_mock_files)
B64_sens_ppv <- return_sensitivity_and_ppv_64(B64_mock_files)
C64_sens_ppv <- return_sensitivity_and_ppv_64(C64_mock_files)

A23_sens_ppv <- return_sensitivity_and_ppv_23(A23_mock_files)
B23_sens_ppv <- return_sensitivity_and_ppv_23(B23_mock_files)
C23_sens_ppv <- return_sensitivity_and_ppv_23(C23_mock_files)

all_mocks_sens_ppv <- rbind(A64_sens_ppv,B64_sens_ppv,C64_sens_ppv,A23_sens_ppv,B23_sens_ppv,C23_sens_ppv)

all_mocks_sens_ppv$Confidence <- c(seq(0,1, by=0.1))
all_mocks_sens_ppv$Mock <- rep(c("A64","B64","C64","A23","B23","C23"), each=11)
all_mocks_sens_ppv$Species_nr <- c(rep(c(64,23), each = 33))
to_reshape_df <- all_mocks_sens_ppv[,-c(1,2,6)]
ppv_plotting <- reshape2::melt(to_reshape_df, id = c("Confidence","Species_nr"), measure.vars = "Positivevalue")
sens_plotting <- reshape2::melt(to_reshape_df, id = c("Confidence","Species_nr"), measure.vars = "Sensitivity")

ppv_plotting$Species_nr <- factor(ppv_plotting$Species_nr)
sens_plotting$Species_nr <- factor(sens_plotting$Species_nr)


ppv_all_mocks <- ggline(ppv_plotting, x = "Confidence", y = "value", add = "mean_se",
                        color = "Species_nr", palette = "npg")+ ylim(0,0.75) +
  stat_compare_means(label = "p.signif", ref.group = "0", method = "t.test", label.y = 0.7) + ylab("P.P.V.") + scale_x_discrete(labels=c(0,"",0.2,"","",0.5,"","",0.8,"",1)) 



sens_all_mocks <- ggline(sens_plotting, x = "Confidence", y = "value", add = "mean_se",
                     color = "Species_nr", palette = "npg")+
  stat_compare_means(label = "p.signif", ref.group = "0", method = "t.test")  + ylab("Sensitivity") + scale_x_discrete(labels=c(0,"",0.2,"","",0.5,"","",0.8,"",1))


ppv_all_mocks2 <- ggpar(ppv_all_mocks, legend.title = "")
sens_all_mocks2 <- ggpar(sens_all_mocks, legend.title = "")

ggsave(plot=ppv_all_mocks2,filename="ppv_all_mocks.pdf",width = 3, height = 3)
ggsave(plot=sens_all_mocks2,filename="sens_all_mocks.pdf",width = 3, height = 3)

# Get classified percentages 

classified_and_percentage <- function(df) {
  number_classified <- df[2,1]
  percent_classified <- (df[2,1]/(df[2,1]+df[1,1]))*100
  rbind(number_classified,percent_classified)
}


A64_class <- lapply(A64_mock_files,classified_and_percentage)
A64_df <- do.call(rbind,A64_class) %>% as.data.frame() %>% rename(., Value = V1)
A64_df$Type <- rep(c("Classified","Percent"))
A64_df$Confidence <- rep(seq(from=0.0, to = 1.0, by=0.1), each=2)
rownames(A64_df) <- NULL
A64_df$Species_nr <- 64
A64_df$Mock <- rep("A64")


total_class_percent <- function(y, species.nr = is.numeric) {
  class <- lapply(y,classified_and_percentage)
  y_df <- do.call(rbind,class) %>% as.data.frame() %>% rename(., Value = V1)
  y_df$Type <- rep(c("Classified","Percent"))
  y_df$Confidence <- rep(seq(from=0.0, to = 1.0, by=0.1), each=2)
  rownames(y_df) <- NULL
  y_df$Species_nr <- species.nr;
  return(y_df)
}

A64_total <- total_class_percent(A64_mock_files, species.nr = 64)
B64_total <- total_class_percent(B64_mock_files, species.nr = 64)
C64_total <- total_class_percent(C64_mock_files, species.nr = 64)

A23_total <- total_class_percent(A23_mock_files, species.nr = 23)
B23_total <- total_class_percent(B23_mock_files, species.nr = 23)
C23_total <- total_class_percent(C23_mock_files, species.nr = 23)

all_mocks_classified_percent <- rbind(A64_total,B64_total,C64_total,A23_total,B23_total,C23_total)
all_mocks_classified_percent$Species_nr <- as.factor(all_mocks_classified_percent$Species_nr)

all_mocks_classified_only <- all_mocks_classified_percent %>% filter(., Type == "Classified")
all_mocks_percent_only <- all_mocks_classified_percent %>% filter(., Type == "Percent")

classified_all_mocks <- ggscatter(all_mocks_classified_only, x = "Confidence", y = "Value",
                                  add = "loess", conf.int = TRUE, color = "Species_nr", palette = get_palette("npg",2)) + ylab("Classified reads") + scale_y_continuous(labels = scales::comma)
classified_all_mocks2 <- ggpar(classified_all_mocks,legend.title = "Number of species per sample")

ggsave(plot=classified_all_mocks2, filename="classified_reads_all_mocks.pdf")

percent_all_mocks <- ggscatter(all_mocks_percent_only, x = "Confidence", y = "Value",
                               add = "loess", conf.int = TRUE, color = "Species_nr", palette = get_palette("npg",2)) + ylab("Classified reads (%)") + scale_x_continuous(labels = c(0,"",0.2,"","",0.5,"","",0.8,"",1), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

percent_all_mocks2 <- ggpar(percent_all_mocks,legend.title = "", ylim = c(0,100))

ggsave(plot=percent_all_mocks2, filename="percent_classified_reads_all_mocks.pdf", width = 3, height = 3)

# library size comparison. We will use the previous library from A64/B64/C64
kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_A8M <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A8M",full.names = T, pattern = "\\.kraken2")

A8M_mock <- lapply(import_A8M, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_A8M <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A8M", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

A8M_mock_files <- setNames(A8M_mock, substring(naming_list_A8M, first  = 1, last = nchar(naming_list_A8M) -8))

A8M_mock_files = lapply(A8M_mock_files,setNames,kraken2_output_names)

A8M_mock_files = lapply(A8M_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
A8M_mock_files = lapply(A8M_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_A8M <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_A8M.csv") 
composition_A8M$True_relative_abundance <- (composition_A8M$True_reads/sum(composition_A8M$True_reads)) * 100
composition_A8M <- composition_A8M[,-c(2,3)]
A8M_filtering <- composition_A8M$Taxonomic_ID

A8M_filtered <- A8M_mock_files %>% lapply(filter, Taxonomic_ID %in% A8M_filtering)

A8M_filtered_2 <- lapply(A8M_filtered, "[", c(4,6))

A8M_mock_final <- lapply(names(A8M_filtered_2), function(x){
  colnames(A8M_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  A8M_filtered_2[[x]]
})
names(A8M_mock_final) <- names(A8M_filtered_2) 


A8M_confidence <- A8M_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

A8M_mock_with_confidence <- merge(A8M_confidence,composition_A8M, by="Taxonomic_ID")

rmse_A8M <- get_rmse_eval(A8M_mock_with_confidence)



kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_B8M <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B8M",full.names = T, pattern = "\\.kraken2")

B8M_mock <- lapply(import_B8M, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_B8M <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B8M", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

B8M_mock_files <- setNames(B8M_mock, substring(naming_list_B8M, first  = 1, last = nchar(naming_list_B8M) -8))

B8M_mock_files = lapply(B8M_mock_files,setNames,kraken2_output_names)

B8M_mock_files = lapply(B8M_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
B8M_mock_files = lapply(B8M_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_B8M <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_B8M.csv") 
composition_B8M$True_relative_abundance <- (composition_B8M$True_reads/sum(composition_B8M$True_reads)) * 100
composition_B8M <- composition_B8M[,-c(2,3)]
B8M_filtering <- composition_B8M$Taxonomic_ID

B8M_filtered <- B8M_mock_files %>% lapply(filter, Taxonomic_ID %in% B8M_filtering)

B8M_filtered_2 <- lapply(B8M_filtered, "[", c(4,6))

B8M_mock_final <- lapply(names(B8M_filtered_2), function(x){
  colnames(B8M_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  B8M_filtered_2[[x]]
})
names(B8M_mock_final) <- names(B8M_filtered_2) 


B8M_confidence <- B8M_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

B8M_mock_with_confidence <- merge(B8M_confidence,composition_B8M, by="Taxonomic_ID")
rmse_B8M <- get_rmse_eval(B8M_mock_with_confidence)


kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_C8M <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C8M",full.names = T, pattern = "\\.kraken2")

C8M_mock <- lapply(import_C8M, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_C8M <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C8M", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

C8M_mock_files <- setNames(C8M_mock, substring(naming_list_C8M, first  = 1, last = nchar(naming_list_C8M) -8))

C8M_mock_files = lapply(C8M_mock_files,setNames,kraken2_output_names)

C8M_mock_files = lapply(C8M_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
C8M_mock_files = lapply(C8M_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_C8M <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_C8M.csv") 
composition_C8M$True_relative_abundance <- (composition_C8M$True_reads/sum(composition_C8M$True_reads)) * 100
composition_C8M <- composition_C8M[,-c(2,3)]
C8M_filtering <- composition_C8M$Taxonomic_ID

C8M_filtered <- C8M_mock_files %>% lapply(filter, Taxonomic_ID %in% C8M_filtering)

C8M_filtered_2 <- lapply(C8M_filtered, "[", c(4,6))

C8M_mock_final <- lapply(names(C8M_filtered_2), function(x){
  colnames(C8M_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  C8M_filtered_2[[x]]
})
names(C8M_mock_final) <- names(C8M_filtered_2) 


C8M_confidence <- C8M_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

C8M_mock_with_confidence <- merge(C8M_confidence,composition_C8M, by="Taxonomic_ID")
rmse_C8M <- get_rmse_eval(C8M_mock_with_confidence)

all_rmse_values_library_size <- cbind(rmse_A64,rmse_B64,rmse_C64,rmse_A8M,rmse_B8M,rmse_C8M)
colnames(all_rmse_values_library_size) <- c("A64","B64","C64","A8M","B8M","C8M")
all_rmse_values_library_size$Confidence <- c(seq(0,1, by=0.1))
all_rmse_values_library_size_plotting <- reshape2::melt(all_rmse_values_library_size, id = "Confidence")
all_rmse_values_library_size_plotting$Library_size <- c(rep(c(1,8), each = 33))
all_rmse_values_library_size_plotting$Library_size <- as.factor(all_rmse_values_library_size_plotting$Library_size)
colnames(all_rmse_values_library_size_plotting) <- c("Confidence","Mock","RMSE","Library_size")

boxplot4 <- ggboxplot(all_rmse_values_library_size_plotting, x = "Confidence", y = "RMSE", add = "jitter", color = "Library_size", palette = get_palette("npg",11)) + ylim(0,4) + theme(text = element_text(size=15), axis.text.x = element_text(size=10)) + scale_x_discrete(labels = c(0,"",0.2,"","",0.5,"","",0.8,"",1))   + stat_compare_means(aes(group=Library_size),method = "t.test", label = "p.signif")


boxplot5 <- ggpar(boxplot4, legend.title = "")

ggsave(plot = boxplot5, filename = "RMSE_library_size.pdf", height = 3, width = 3)


# number of species

species_64 <- c(A64_mock_files,B64_mock_files,C64_mock_files)

numberspecies64 <- species_64 %>% lapply(filter, Rank_code == "S")
numberspecies641 <- as.data.frame(sapply(numberspecies64, nrow)) %>% tibble::rownames_to_column(var = "Sample")
numberspecies641$Mock <- rep(c("A64","B64","C64"),each=11)
numberspecies641$Confidence <- rep(seq(from=0.0,to=1,by=0.1))
numberspecies641_new <- numberspecies641 %>% rename(., Number_of_species = names(numberspecies641[2])) %>% select(., -Sample)
numberspecies641_new$Mock <- factor(numberspecies641$Mock)

species_8m <- c(A8M_mock_files,B8M_mock_files,C8M_mock_files)
numberspecies8m <- species_8m %>% lapply(filter, Rank_code == "S")
numberspecies8m1 <- as.data.frame(sapply(numberspecies8m, nrow)) %>% tibble::rownames_to_column(var = "Sample")
numberspecies8m1$Mock <- rep(c("A8M","B8M","C8M"),each=11)
numberspecies8m1$Confidence <- rep(seq(from=0.0,to=1,by=0.1))
numberspecies8m1_new <- numberspecies8m1 %>% rename(., Number_of_species = names(numberspecies8m1[2])) %>% select(., -Sample)
numberspecies8m1_new$Mock <- factor(numberspecies8m1$Mock)

all_mocks_species_library_size <- rbind(numberspecies641_new,numberspecies8m1_new)

all_mocks_species_plot_library_size <- ggline(all_mocks_species_library_size, "Confidence", "Number_of_species",
                                              linetype = "Mock", shape = "Mock",
                                              color = "Mock") + geom_hline(yintercept = 64, linetype = "dashed") + theme_bw() + theme(text = element_text(size=15))
all_mocks_species_plot_library_size2 <- ggpar(all_mocks_species_plot_library_size, ylab = "Number of species identified")

ggsave(plot=all_mocks_species_plot_library_size2,filename="dataplotspecies_new_mocks_library_size.pdf", height = 6.5, width = 4)


# Sensitivity for 1 million vs 8 million reads

counting_64_species <- function(df) {
  taxids <- df %>% dplyr::pull(.,Taxonomic_ID) %>% as.character()
  df_1 <- taxids %in% A64_filtering %>% as.data.frame() %>% dplyr::rename(.,Presence = .) %>% cbind(taxids,.) 
}

count_true_false_species <- function(df) {
  vector1 <- df[,2];
  return(sum(vector1))
}


A64_sens_ppv <- return_sensitivity_and_ppv_64(A64_mock_files)
B64_sens_ppv <- return_sensitivity_and_ppv_64(B64_mock_files)
C64_sens_ppv <- return_sensitivity_and_ppv_64(C64_mock_files)

A8M_sens_ppv <- return_sensitivity_and_ppv_64(A8M_mock_files)
B8M_sens_ppv <- return_sensitivity_and_ppv_64(B8M_mock_files)
C8M_sens_ppv <- return_sensitivity_and_ppv_64(C8M_mock_files)

all_mocks_sens_ppv_library_size <- rbind(A64_sens_ppv,B64_sens_ppv,C64_sens_ppv,A8M_sens_ppv,B8M_sens_ppv,C8M_sens_ppv)

all_mocks_sens_ppv_library_size$Confidence <- c(seq(0,1, by=0.1))
all_mocks_sens_ppv_library_size$Mock <- rep(c("A64","B64","C64","A8M","B8M","C8M"), each=11)
all_mocks_sens_ppv_library_size$Library_size <- c(rep(c(1,8), each = 33))
all_mocks_sens_ppv_library_size$Library_size <- as.factor(all_mocks_sens_ppv_library_size$Library_size)
to_reshape_df_library_size <- all_mocks_sens_ppv_library_size[,-c(1,2,6)]
ppv_plotting_library_size <- reshape2::melt(to_reshape_df_library_size, id = c("Confidence","Library_size"), measure.vars = "Positivevalue")
sens_plotting_library_size <- reshape2::melt(to_reshape_df_library_size, id = c("Confidence","Library_size"), measure.vars = "Sensitivity")


ppv_all_mocks_library_size <- ggline(ppv_plotting_library_size, x = "Confidence", y = "value", add = "mean_se",
                                     color = "Library_size", palette = "npg")+
  stat_compare_means(label = "p.signif", ref.group = "0", method = "t.test")  + ylab("P.P.V.") + scale_x_discrete(labels=c(0,"",0.2,"","",0.5,"","",0.8,"",1))

  
sens_all_mocks_library_size <- ggline(sens_plotting_library_size, x = "Confidence", y = "value", add = "mean_se",
                                        color = "Library_size", palette = "npg")+
    stat_compare_means(label = "p.signif", ref.group = "0", method = "t.test")  + ylab("P.P.V.") + scale_x_discrete(labels=c(0,"",0.2,"","",0.5,"","",0.8,"",1))
  
    
    
ppv_all_mocks_library_size2 <- ggpar(ppv_all_mocks_library_size, legend.title = "")
sens_all_mocks_library_size2 <- ggpar(sens_all_mocks_library_size, legend.title = "")

ggsave(plot=ppv_all_mocks_library_size2,filename="ppv_all_mocks_library_size.pdf",width = 3,height = 3)
ggsave(plot=sens_all_mocks_library_size2,filename="sens_all_mocks_library_size.pdf",width = 3,height = 3)

# classified library size



classified_and_percentage <- function(df) {
  number_classified <- df[2,1]
  percent_classified <- (df[2,1]/(df[2,1]+df[1,1]))*100
  rbind(number_classified,percent_classified)
}


total_class_percent_library_size <- function(y, species.nr = is.numeric, library.size = is.numeric) {
  class <- lapply(y,classified_and_percentage)
  y_df <- do.call(rbind,class) %>% as.data.frame() %>% rename(., Value = V1)
  y_df$Type <- rep(c("Classified","Percent"))
  y_df$Confidence <- rep(seq(from=0.0, to = 1.0, by=0.1), each=2)
  rownames(y_df) <- NULL
  y_df$Species_nr <- species.nr
  y_df$Library_size <- library.size;
  return(y_df)
}

A64_total_lib <- total_class_percent_library_size(A64_mock_files, species.nr = 64, library.size = 1)
B64_total_lib <- total_class_percent_library_size(B64_mock_files, species.nr = 64, library.size = 1)
C64_total_lib <- total_class_percent_library_size(C64_mock_files, species.nr = 64, library.size = 1)

A8M_total <- total_class_percent_library_size(A8M_mock_files, species.nr = 64, library.size = 8)
B8M_total <- total_class_percent_library_size(B8M_mock_files, species.nr = 64, library.size = 8)
C8M_total <- total_class_percent_library_size(C8M_mock_files, species.nr = 64, library.size = 8)

all_mocks_classified_percent_library_size <- rbind(A64_total_lib,B64_total_lib,C64_total_lib,A8M_total,B8M_total,C8M_total)
all_mocks_classified_percent_library_size$Library_size <- as.factor(all_mocks_classified_percent_library_size$Library_size)

all_mocks_classified_only_library_size <- all_mocks_classified_percent_library_size %>% filter(., Type == "Classified")
all_mocks_percent_only_library_size <- all_mocks_classified_percent_library_size %>% filter(., Type == "Percent")

classified_all_mocks_library_size <- ggscatter(all_mocks_classified_only_library_size, x = "Confidence", y = "Value",
                                               add = "loess", conf.int = TRUE, color = "Library_size", palette = get_palette("npg",2)) + ylab("Classified reads")
classified_all_mocks2_library_size <- ggpar(classified_all_mocks_library_size,legend.title = "Library size (Million reads)")

ggsave(plot=classified_all_mocks2_library_size, filename="classified_reads_all_mocks_library_size.png")

scattertry <- ggscatter(all_mocks_classified_only_library_size, x = "Confidence", y = "Value",
                        add = "loess", conf.int = TRUE, color = "Library_size", palette = get_palette("npg",2)) + ylab("Classified reads")+ stat_compare_means(label = "p.signif", method = "t.test", ref.group = ".all.")


ggsave(plot = boxplot_rmse_fragment, filename = "RMSE_fragment_size_with_stat.png")



percent_all_mocks_library_size <- ggscatter(all_mocks_percent_only_library_size, x = "Confidence", y = "Value",
                                            add = "loess", conf.int = TRUE, color = "Library_size", palette = get_palette("npg",2)) + ylab("Classified reads (%)")+ scale_x_continuous(labels = c(0,"",0.2,"","",0.5,"","",0.8,"",1), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

percent_all_mocks2_library_size <- ggpar(percent_all_mocks_library_size,legend.title = "", ylim = c(0,100))

ggsave(plot=percent_all_mocks2_library_size, filename="percent_classified_reads_all_mocks)library_size.pdf",height = 3,width = 3)



# Fragment size

kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_A200 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A200",full.names = T, pattern = "\\.kraken2")

A200_mock <- lapply(import_A200, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_A200 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A200", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

A200_mock_files <- setNames(A200_mock, substring(naming_list_A200, first  = 1, last = nchar(naming_list_A200) -8))

A200_mock_files = lapply(A200_mock_files,setNames,kraken2_output_names)

A200_mock_files = lapply(A200_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
A200_mock_files = lapply(A200_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_A200 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_A200.csv") 
composition_A200$True_relative_abundance <- (composition_A200$True_reads/sum(composition_A200$True_reads)) * 100
composition_A200 <- composition_A200[,-c(2,3)]
A200_filtering <- composition_A200$Taxonomic_ID

A200_filtered <- A200_mock_files %>% lapply(filter, Taxonomic_ID %in% A200_filtering)

A200_filtered_2 <- lapply(A200_filtered, "[", c(4,6))

A200_mock_final <- lapply(names(A200_filtered_2), function(x){
  colnames(A200_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  A200_filtered_2[[x]]
})
names(A200_mock_final) <- names(A200_filtered_2) 


A200_confidence <- A200_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

A200_mock_with_confidence <- merge(A200_confidence,composition_A200, by="Taxonomic_ID")

rmse_A200 <- get_rmse_eval(A200_mock_with_confidence)


kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_B200 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B200",full.names = T, pattern = "\\.kraken2")

B200_mock <- lapply(import_B200, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_B200 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B200", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

B200_mock_files <- setNames(B200_mock, substring(naming_list_B200, first  = 1, last = nchar(naming_list_B200) -8))

B200_mock_files = lapply(B200_mock_files,setNames,kraken2_output_names)

B200_mock_files = lapply(B200_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
B200_mock_files = lapply(B200_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_B200 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_B200.csv") 
composition_B200$True_relative_abundance <- (composition_B200$True_reads/sum(composition_B200$True_reads)) * 100
composition_B200 <- composition_B200[,-c(2,3)]
B200_filtering <- composition_B200$Taxonomic_ID

B200_filtered <- B200_mock_files %>% lapply(filter, Taxonomic_ID %in% B200_filtering)

B200_filtered_2 <- lapply(B200_filtered, "[", c(4,6))

B200_mock_final <- lapply(names(B200_filtered_2), function(x){
  colnames(B200_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  B200_filtered_2[[x]]
})
names(B200_mock_final) <- names(B200_filtered_2) 


B200_confidence <- B200_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

B200_mock_with_confidence <- merge(B200_confidence,composition_B200, by="Taxonomic_ID")

rmse_B200 <- get_rmse_eval(B200_mock_with_confidence)


kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_C200 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C200",full.names = T, pattern = "\\.kraken2")

C200_mock <- lapply(import_C200, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_C200 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C200", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

C200_mock_files <- setNames(C200_mock, substring(naming_list_C200, first  = 1, last = nchar(naming_list_C200) -8))

C200_mock_files = lapply(C200_mock_files,setNames,kraken2_output_names)

C200_mock_files = lapply(C200_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
C200_mock_files = lapply(C200_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_C200 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_C200.csv") 
composition_C200$True_relative_abundance <- (composition_C200$True_reads/sum(composition_C200$True_reads)) * 100
composition_C200 <- composition_C200[,-c(2,3)]
C200_filtering <- composition_C200$Taxonomic_ID

C200_filtered <- C200_mock_files %>% lapply(filter, Taxonomic_ID %in% C200_filtering)

C200_filtered_2 <- lapply(C200_filtered, "[", c(4,6))

C200_mock_final <- lapply(names(C200_filtered_2), function(x){
  colnames(C200_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  C200_filtered_2[[x]]
})
names(C200_mock_final) <- names(C200_filtered_2) 


C200_confidence <- C200_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

C200_mock_with_confidence <- merge(C200_confidence,composition_C200, by="Taxonomic_ID")

rmse_C200 <- get_rmse_eval(C200_mock_with_confidence)


kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_A110 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A110",full.names = T, pattern = "\\.kraken2")

A110_mock <- lapply(import_A110, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_A110 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/A110", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

A110_mock_files <- setNames(A110_mock, substring(naming_list_A110, first  = 1, last = nchar(naming_list_A110) -8))

A110_mock_files = lapply(A110_mock_files,setNames,kraken2_output_names)

A110_mock_files = lapply(A110_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
A110_mock_files = lapply(A110_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_A110 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_A110.csv") 
composition_A110$True_relative_abundance <- (composition_A110$True_reads/sum(composition_A110$True_reads)) * 100
composition_A110 <- composition_A110[,-c(2,3)]
A110_filtering <- composition_A110$Taxonomic_ID

A110_filtered <- A110_mock_files %>% lapply(filter, Taxonomic_ID %in% A110_filtering)

A110_filtered_2 <- lapply(A110_filtered, "[", c(4,6))

A110_mock_final <- lapply(names(A110_filtered_2), function(x){
  colnames(A110_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  A110_filtered_2[[x]]
})
names(A110_mock_final) <- names(A110_filtered_2) 


A110_confidence <- A110_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

A110_mock_with_confidence <- merge(A110_confidence,composition_A110, by="Taxonomic_ID")

rmse_A110 <- get_rmse_eval(A110_mock_with_confidence)


kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_B110 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B110",full.names = T, pattern = "\\.kraken2")

B110_mock <- lapply(import_B110, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_B110 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/B110", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

B110_mock_files <- setNames(B110_mock, substring(naming_list_B110, first  = 1, last = nchar(naming_list_B110) -8))

B110_mock_files = lapply(B110_mock_files,setNames,kraken2_output_names)

B110_mock_files = lapply(B110_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
B110_mock_files = lapply(B110_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_B110 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_B110.csv") 
composition_B110$True_relative_abundance <- (composition_B110$True_reads/sum(composition_B110$True_reads)) * 100
composition_B110 <- composition_B110[,-c(2,3)]
B110_filtering <- composition_B110$Taxonomic_ID

B110_filtered <- B110_mock_files %>% lapply(filter, Taxonomic_ID %in% B110_filtering)

B110_filtered_2 <- lapply(B110_filtered, "[", c(4,6))

B110_mock_final <- lapply(names(B110_filtered_2), function(x){
  colnames(B110_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  B110_filtered_2[[x]]
})
names(B110_mock_final) <- names(B110_filtered_2) 


B110_confidence <- B110_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

B110_mock_with_confidence <- merge(B110_confidence,composition_B110, by="Taxonomic_ID")

rmse_B110 <- get_rmse_eval(B110_mock_with_confidence)


kraken2_output_names <- c("Sample_rooted_reads", "Sample_reads","Rank_code","Taxonomic_ID","Name")

import_C110 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C110",full.names = T, pattern = "\\.kraken2")

C110_mock <- lapply(import_C110, function(tble) {
  read.delim(file = tble, 
             sep = '\t',
             header = FALSE,
             strip.white = TRUE)[,c(-1)]
})

naming_list_C110 <- list.files(path = "./Simulated_mocks/kraken2_threshold_evaluation/C110", pattern = "\\.kraken2")

# This will name each element (or dataframe) of the list
# according to the filename it originated from.

C110_mock_files <- setNames(C110_mock, substring(naming_list_C110, first  = 1, last = nchar(naming_list_C110) -8))

C110_mock_files = lapply(C110_mock_files,setNames,kraken2_output_names)

C110_mock_files = lapply(C110_mock_files,arrange, Taxonomic_ID)


# Optionally add the relative abundance as well
C110_mock_files = lapply(C110_mock_files, function(df) {
  df$Relative_abundance = df$Sample_reads / df[2,1] * 100;
  df$Relative_abundance_rooted = df$Sample_rooted_reads / df[1,2] * 100;
  return(df)
})

composition_C110 <- read.csv("./Simulated_mocks/kraken2_threshold_evaluation/Mock_C110.csv") 
composition_C110$True_relative_abundance <- (composition_C110$True_reads/sum(composition_C110$True_reads)) * 100
composition_C110 <- composition_C110[,-c(2,3)]
C110_filtering <- composition_C110$Taxonomic_ID

C110_filtered <- C110_mock_files %>% lapply(filter, Taxonomic_ID %in% C110_filtering)

C110_filtered_2 <- lapply(C110_filtered, "[", c(4,6))

C110_mock_final <- lapply(names(C110_filtered_2), function(x){
  colnames(C110_filtered_2[[x]]) <- c("Taxonomic_ID",x)
  C110_filtered_2[[x]]
})
names(C110_mock_final) <- names(C110_filtered_2) 


C110_confidence <- C110_mock_final %>% purrr::reduce(full_join, by = "Taxonomic_ID") %>% select("Taxonomic_ID", everything()) %>% mutate_all(~replace(., is.na(.), 0))

C110_mock_with_confidence <- merge(C110_confidence,composition_C110, by="Taxonomic_ID")

rmse_C110 <- get_rmse_eval(C110_mock_with_confidence)

all_rmse_values_fragment <- cbind(rmse_A110,rmse_B110,rmse_C110,rmse_A200,rmse_B200,rmse_C200)
colnames(all_rmse_values_fragment) <- c("A110","B110","C110","A200","B200","C200")
all_rmse_values_fragment$Confidence <- c(seq(0,1, by=0.1))
all_rmse_values_fragment_plotting <- reshape2::melt(all_rmse_values_fragment, id = "Confidence")
all_rmse_values_fragment_plotting$Fragment_size <- c(rep(c(110,200), each = 33))
all_rmse_values_fragment_plotting$Fragment_size <- as.factor(all_rmse_values_fragment_plotting$Fragment_size)
colnames(all_rmse_values_fragment_plotting) <- c("Confidence","Mock","RMSE","Fragment_size")

boxplot6 <- ggboxplot(all_rmse_values_fragment_plotting, x = "Confidence", y = "RMSE", add = "jitter", color = "Fragment_size", palette = get_palette("npg",2)) + ylim(0,5) + theme(text = element_text(size=12), axis.text.x = element_text(size=10), legend.position = "top") + scale_x_discrete(labels = c(0,"",0.2,"","",0.5,"","",0.8,"",1))  + stat_compare_means(aes(group=Fragment_size),method = "t.test", label = "p.signif")


boxplot7 <- ggpar(boxplot6, legend.title = "")

ggsave(plot = boxplot7, filename = "RMSE_fragment_size.pdf",height = 3,width = 3)


boxplot_rmse_fragment <- ggboxplot(all_rmse_values_fragment_plotting, x = "Confidence", y = "RMSE",
                                   color = "Fragment_size", palette = "jco") + ylim(0,4) +      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")
boxplot_rmse_fragment <- ggpar(boxplot_rmse_fragment, legend.title = "Fragment size (bp)")

ggsave(plot = boxplot_rmse_fragment, filename = "RMSE_fragment_size_with_stat.png")


# number of species fragment size

species_110 <- c(A110_mock_files,B110_mock_files,C110_mock_files)
species_200 <- c(A200_mock_files,B200_mock_files,C200_mock_files)

numberspecies110 <- species_110 %>% lapply(filter, Rank_code == "S")
numberspecies1101 <- as.data.frame(sapply(numberspecies110, nrow)) %>% tibble::rownames_to_column(var = "Sample")
numberspecies1101$Mock <- rep(c("A110","B110","C110"),each=11)
numberspecies1101$Confidence <- rep(seq(from=0.0,to=1,by=0.1))
numberspecies1101_new <- numberspecies1101 %>% rename(., Number_of_species = names(numberspecies1101[2])) %>% select(., -Sample)
numberspecies1101_new$Mock <- factor(numberspecies1101$Mock)

numberspecies200 <- species_200 %>% lapply(filter, Rank_code == "S")
numberspecies2001 <- as.data.frame(sapply(numberspecies200, nrow)) %>% tibble::rownames_to_column(var = "Sample")
numberspecies2001$Mock <- rep(c("A200","B200","C200"),each=11)
numberspecies2001$Confidence <- rep(seq(from=0.0,to=1,by=0.1))
numberspecies2001_new <- numberspecies2001 %>% rename(., Number_of_species = names(numberspecies2001[2])) %>% select(., -Sample)
numberspecies2001_new$Mock <- factor(numberspecies2001$Mock)

all_mocks_species_fragment_size <- rbind(numberspecies1101_new,numberspecies2001_new)

all_mocks_species_plot_fragment_size <- ggline(all_mocks_species_fragment_size, "Confidence", "Number_of_species",
                                               linetype = "Mock", shape = "Mock",
                                               color = "Mock") + geom_hline(yintercept = 64, linetype = "dashed") + theme_bw() + theme(text = element_text(size=15))
all_mocks_species_plot_fragment_size2 <- ggpar(all_mocks_species_plot_fragment_size, ylab = "Number of species identified")

ggsave(plot=all_mocks_species_plot_fragment_size2,filename="dataplotspecies_new_mocks_fragment_size.pdf",height = 6.5,width = 4)


# Sensitivity for fragment size

counting_64_species <- function(df) {
  taxids <- df %>% dplyr::pull(.,Taxonomic_ID) %>% as.character()
  df_1 <- taxids %in% A64_filtering %>% as.data.frame() %>% dplyr::rename(.,Presence = .) %>% cbind(taxids,.) 
}

count_true_false_species <- function(df) {
  vector1 <- df[,2];
  return(sum(vector1))
}


A110_sens_ppv <- return_sensitivity_and_ppv_64(A110_mock_files)
B110_sens_ppv <- return_sensitivity_and_ppv_64(B110_mock_files)
C110_sens_ppv <- return_sensitivity_and_ppv_64(C110_mock_files)

A200_sens_ppv <- return_sensitivity_and_ppv_64(A200_mock_files)
B200_sens_ppv <- return_sensitivity_and_ppv_64(B200_mock_files)
C200_sens_ppv <- return_sensitivity_and_ppv_64(C200_mock_files)

all_mocks_sens_ppv_fragment_size <- rbind(A110_sens_ppv,B110_sens_ppv,C110_sens_ppv,A200_sens_ppv,B200_sens_ppv,C200_sens_ppv)

all_mocks_sens_ppv_fragment_size$Confidence <- c(seq(0,1, by=0.1))
all_mocks_sens_ppv_fragment_size$Mock <- rep(c("A110","B110","C110","A200","B200","C200"), each=11)
all_mocks_sens_ppv_fragment_size$Fragment_size <- c(rep(c(110,200), each = 33))
all_mocks_sens_ppv_fragment_size$Fragment_size <- as.factor(all_mocks_sens_ppv_fragment_size$Fragment_size)
to_reshape_df_fragment_size <- all_mocks_sens_ppv_fragment_size[,-c(1,2,6)]
ppv_plotting_fragment_size <- reshape2::melt(to_reshape_df_fragment_size, id = c("Confidence","Fragment_size"), measure.vars = "Positivevalue")
sens_plotting_fragment_size <- reshape2::melt(to_reshape_df_fragment_size, id = c("Confidence","Fragment_size"), measure.vars = "Sensitivity")

ppv_all_mocks_fragment_size <- ggline(ppv_plotting_fragment_size, x = "Confidence", y = "value", add = "mean_se",
                                      color = "Fragment_size", palette = "npg")+
  stat_compare_means(label = "p.signif", ref.group = "0", method = "t.test")  + ylab("P.P.V.") + scale_x_discrete(labels=c(0,"",0.2,"","",0.5,"","",0.8,"",1))


sens_all_mocks_fragment_size <- ggline(sens_plotting_fragment_size, x = "Confidence", y = "value", add = "mean_se",
                                       color = "Fragment_size", palette = "npg")+
  stat_compare_means(label = "p.signif", ref.group = "0", method = "t.test")  + ylab("Sensitivity") + scale_x_discrete(labels=c(0,"",0.2,"","",0.5,"","",0.8,"",1))

ppv_all_mocks_fragment_size2 <- ggpar(ppv_all_mocks_fragment_size, legend.title = "Fragment size (bp)")
sens_all_mocks_fragment_size2 <- ggpar(sens_all_mocks_fragment_size, legend.title = "Fragment size (bp)")

ggsave(plot=ppv_all_mocks_fragment_size2,filename="ppv_all_mocks_fragment_size.pdf",width = 3,height = 3)
ggsave(plot=sens_all_mocks_fragment_size2,filename="sens_all_mocks_fragment_size.pdf",width = 3,height = 3)


classified_and_percentage <- function(df) {
  number_classified <- df[2,1]
  percent_classified <- (df[2,1]/(df[2,1]+df[1,1]))*100
  rbind(number_classified,percent_classified)
}


total_class_percent_fragment_size <- function(y, species.nr = is.numeric, fragment.size = is.numeric) {
  class <- lapply(y,classified_and_percentage)
  y_df <- do.call(rbind,class) %>% as.data.frame() %>% rename(., Value = V1)
  y_df$Type <- rep(c("Classified","Percent"))
  y_df$Confidence <- rep(seq(from=0.0, to = 1.0, by=0.1), each=2)
  rownames(y_df) <- NULL
  y_df$Species_nr <- species.nr
  y_df$Fragment_size <- fragment.size;
  return(y_df)
}

A110_total <- total_class_percent_fragment_size(A110_mock_files, species.nr = 64, fragment.size = 110)
B110_total <- total_class_percent_fragment_size(B110_mock_files, species.nr = 64, fragment.size = 110)
C110_total <- total_class_percent_fragment_size(C110_mock_files, species.nr = 64, fragment.size = 110)

A200_total <- total_class_percent_fragment_size(A200_mock_files, species.nr = 64, fragment.size = 200)
B200_total <- total_class_percent_fragment_size(B200_mock_files, species.nr = 64, fragment.size = 200)
C200_total <- total_class_percent_fragment_size(C200_mock_files, species.nr = 64, fragment.size = 200)


all_mocks_classified_percent_fragment_size <- rbind(A110_total,B110_total,A200_total,B200_total)
all_mocks_classified_percent_fragment_size$Fragment_size <- as.factor(all_mocks_classified_percent_fragment_size$Fragment_size)

all_mocks_classified_only_fragment_size <- all_mocks_classified_percent_fragment_size %>% filter(., Type == "Classified")
all_mocks_percent_only_fragment_size <- all_mocks_classified_percent_fragment_size %>% filter(., Type == "Percent")

classified_all_mocks_fragment_size <- ggscatter(all_mocks_classified_only_fragment_size, x = "Confidence", y = "Value",
                                                add = "loess", conf.int = TRUE, color = "Fragment_size", palette = get_palette("npg",2)) + ylab("Classified reads")
classified_all_mocks2_fragment_size <- ggpar(classified_all_mocks_fragment_size,legend.title = "Fragment size (bp)")

ggsave(plot=classified_all_mocks2_fragment_size, filename="classified_reads_all_mocks_fragment_size.png")

percent_all_mocks_fragment_size <- ggscatter(all_mocks_percent_only_fragment_size, x = "Confidence", y = "Value",
                                             add = "loess", conf.int = TRUE, color = "Fragment_size", palette = get_palette("npg",2)) + ylab("Classified reads (%)")+ scale_x_continuous(labels = c(0,"",0.2,"","",0.5,"","",0.8,"",1), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))

percent_all_mocks2_fragment_size <- ggpar(percent_all_mocks_fragment_size,legend.title = "Fragment size (bp)", ylim = c(0,100))

ggsave(plot=percent_all_mocks2_fragment_size, filename="percent_classified_reads_all_mocks_fragment_size.pdf",height = 3,width = 3)
