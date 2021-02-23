# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Solenn Patalano, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figures 6 A and B
# This script is written in R version 4.0.3

# Load or install (if not present) the required packages
if (!require('readxl')) install.packages('readxl'); library(readxl) 
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library(RColorBrewer)
if (!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if (!require('dplyr')) install.packages('dplyr');library(dplyr)
if (!require('tidyr')) install.packages('tidyr');library(tidyr)
if (!require('purrr')) install.packages('purrr');library(purrr)
if (!require('data.table')) install.packages('data.table');library(data.table)
if (!require('ggpubr')) install.packages('ggpubr');library(ggpubr)
if (!require('car')) install.packages('car'); library(car) #test for normality
if (!require('rcompanion')) install.packages('rcompanion'); library(rcompanion) #for transformation cube square
if (!require('here')) install.packages('here'): library('here')


#Import data
# Varroa_genus is the metagenomic data at genus level
varroa2<-read_excel("Varroa_genus.xlsx")
varroa2<-varroa2[,-c(2:3,6:7,13)]

# clean matrix
varroa2_mat <- dplyr::select_if(varroa2, is.numeric)


# log2 transform CPMs
# I was reading that usually people log2 transform cpm values this way.
# If needed you can apply the same function for Deseq data, I have only log2 transformed them.
log_reads <- function(x) {
  log2(x+1)
}

varroa2_mat<-log_reads(varroa2_mat)

#### OVERALL
# calculate the correlations
r <- cor(varroa2_mat)
round(r,2)

#correlation matrix
ggcorrplot(r, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE)

#####SM
varroa2_mat_SM<-varroa2_mat[,-c(1:2)]
r_SM <- cor(varroa2_mat_SM)
round(r_SM,2)
ggcorrplot(r_SM, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE)
####DSM
varroa2_mat_DSM<-varroa2_mat[,-c(3:4)]
r_DSM <- cor(varroa2_mat_DSM)
round(r_DSM,2)
ggcorrplot(r_DSM, 
           hc.order = TRUE, 
           type = "lower",
           lab = TRUE)


##### FINAL PLOTS 
#import organised  data 
Varroa_final<-read_excel("Varroa_genus.xlsx", sheet="Final", col_types = c("text","text","numeric", "text","numeric"))

Varroa_final<- Varroa_final%>%
  filter(!Monitoring=="average 3 months")

#transform read count in log
Varroa_final$Reads<-log_reads(Varroa_final$Reads)
Varroa_final$Natural_Fall<-log_reads(Varroa_final$Natural_Fall)

#plot
ggscatter( Varroa_final, x = "Reads", y = "Natural_Fall", color = "Monitoring", palette = "jco", add = "reg.line", shape="Hive") +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  facet_wrap(~ Method,scales = "free") +
  theme(plot.title = element_text(size=14), 
        axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))+
  ggtitle("Metagenomic correlation with Natural Varroa fall" ) +
  ylab("Varroa count/day (log)") +
  xlab("Normalised reads count (log)") +
  stat_cor(aes(color = Monitoring))
