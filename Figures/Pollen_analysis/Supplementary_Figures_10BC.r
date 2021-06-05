# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Supplementary Figure 10 B and C
# This script is written in R version 4.0.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('here')) install.packages('here'): library('here')
if (!require('ggpubr')) install.packages('ggpubr'): library('ggpubr')
if (!require('reshape2')) install.packages('reshape2'): library('reshape2')
if (!require('RColorBrewer')) install.packages('RColorBrewer'): library('RColorBrewer')
if (!require('eulerr')) install.packages('eulerr'): library('eulerr')

# Import reads
normalised_species <- as.data.frame(read.csv("./Normalised_reads/normalised_methodseason_family.csv")[,-1])

# Extract plants and calculate relative abundance
normalised_species_plants <- normalised_species %>% filter(Phylum == "Streptophyta",) 
Relative_abundance_species_plants <- normalised_species_plants[,2:9] %>% mutate_all(.,funs((./sum(.))*100))
Relative_abundance_species_plants <- cbind(Relative_abundance_species_plants,normalised_species_plants$Taxonomic_ID,normalised_species_plants$Family)
colnames(Relative_abundance_species_plants)[c(9:10)] <- c("Taxonomic_ID","Family")
Relative_abundance_species_plants$Taxonomic_ID <- as.character(Relative_abundance_species_plants$Taxonomic_ID)
Relative_abundance_species_plants$Family <- as.character(Relative_abundance_species_plants$Family)

# Extract families and summarise at family level (data is already in that format but it's just as QC and because initially i was
# making this with the species data)
data_dsmh4 <- Relative_abundance_species_plants[,c(1,10)] %>% group_by(Family) %>% summarise_at(.vars = "DirectSM_H4", .funs = sum)
data_dsmh5 <- Relative_abundance_species_plants[,c(2,10)] %>% group_by(Family) %>% summarise_at(.vars = "DirectSM_H5", .funs = sum)
data_dsmh6 <- Relative_abundance_species_plants[,c(3,10)] %>% group_by(Family) %>% summarise_at(.vars = "DirectSM_H6", .funs = sum)
data_dsmh7 <- Relative_abundance_species_plants[,c(4,10)] %>% group_by(Family) %>% summarise_at(.vars = "DirectSM_H7", .funs = sum)
data_smh4 <- Relative_abundance_species_plants[,c(5,10)] %>% group_by(Family) %>% summarise_at(.vars = "SM_H4", .funs = sum)
data_smh5 <- Relative_abundance_species_plants[,c(6,10)] %>% group_by(Family) %>% summarise_at(.vars = "SM_H5", .funs = sum)
data_smh6 <- Relative_abundance_species_plants[,c(7,10)] %>% group_by(Family) %>% summarise_at(.vars = "SM_H6", .funs = sum)
data_smh7 <- Relative_abundance_species_plants[,c(8,10)] %>% group_by(Family) %>% summarise_at(.vars = "SM_H7", .funs = sum)


# Import pollen data
hive4_pollen <- as.data.frame(read.csv("./Figures/Figure_11/Hive4_pollen_analysis.csv", fileEncoding="UTF-8-BOM"))
hive5_pollen <- as.data.frame(read.csv("./Figures/Figure_11/Hive5_pollen_analysis.csv", fileEncoding="UTF-8-BOM"))
hive6_pollen <- as.data.frame(read.csv("./Figures/Figure_11/Hive6_pollen_analysis.csv", fileEncoding="UTF-8-BOM"))
hive7_pollen <- as.data.frame(read.csv("./Figures/Figure_11/Hive7_pollen_analysis.csv", fileEncoding="UTF-8-BOM"))

hive4_pollen$Hive <- "Hive 4"
hive5_pollen$Hive <- "Hive 5"
hive6_pollen$Hive <- "Hive 6"
hive7_pollen$Hive <- "Hive 7"

# Combine pollen data
pollen_analysis <- rbind(hive4_pollen,hive5_pollen,hive6_pollen,hive7_pollen)

mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(21)
pollen_analysis$Family <- as.character(pollen_analysis$Family)
pollen_analysis$Hive <- factor(pollen_analysis$Hive, levels = c("Hive 5","Hive 6","Hive 7","Hive 4"))

# Plot the pollen data per hive
all_pollen <- ggplot(pollen_analysis, aes(Hive, Abundance,fill=Family)) + 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values = mycolors)+
  labs(title = "All plant families detected through pollen analysis", y = "Abundance") + scale_y_continuous(labels = scales::percent)

ggsave(filename = "./Figures/Figure_11/all_plants_pollen_analysis.pdf")

# Select the top10 families per sample
top10_pollen_analysis <- pollen_analysis %>% group_by(Hive) %>% arrange(Abundance, .by_group = TRUE) %>% top_n(10, Abundance)

top10_pollen_analysis$Hive <- factor(top10_pollen_analysis$Hive, levels = c("Hive 5","Hive 6","Hive 7","Hive 4"))

top10_pollen_analysis$Family <- as.character(top10_pollen_analysis$Family)

# Plot the top10
top10_pollen <-  ggplot(top10_pollen_analysis, aes(Hive, Abundance,fill=Family)) + 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_manual(values = mycolors)+
  labs(title = "Top10 plant families per sample pollen analysis", y = "Abundance") + scale_y_continuous(labels = scales::percent)

ggsave(filename = "./Figures/Figure_11/top10_pollen_analysis.pdf")


# Extract families detected by each method, which means that the family shouldn't have 0 counts.
data_dsmh4$Method <- "DirectSM"
data_dsmh4$Hive <- "Hive 4"
data_dsmh4 <- data_dsmh4[data_dsmh4$DirectSM_H4 != 0,]
colnames(data_dsmh4)[2] <- "Abundance metagenomic"

data_dsmh5$Method <- "DirectSM"
data_dsmh5$Hive <- "Hive 5"
data_dsmh5 <- data_dsmh5[data_dsmh5$DirectSM_H5 != 0,]
colnames(data_dsmh5)[2] <- "Abundance metagenomic"

data_dsmh6$Method <- "DirectSM"
data_dsmh6$Hive <- "Hive 6"
data_dsmh6 <- data_dsmh6[data_dsmh6$DirectSM_H6 != 0,]
colnames(data_dsmh6)[2] <- "Abundance metagenomic"

data_dsmh7$Method <- "DirectSM"
data_dsmh7$Hive <- "Hive 7"
data_dsmh7 <- data_dsmh7[data_dsmh7$DirectSM_H7 != 0,]
colnames(data_dsmh7)[2] <- "Abundance metagenomic"

data_smh4$Method <- "SM"
data_smh4$Hive <- "Hive 4"
data_smh4 <- data_smh4[data_smh4$SM_H4 != 0,]
colnames(data_smh4)[2] <- "Abundance metagenomic"

data_smh5$Method <- "SM"
data_smh5$Hive <- "Hive 5"
data_smh5 <- data_smh5[data_smh5$SM_H5 != 0,]
colnames(data_smh5)[2] <- "Abundance metagenomic"

data_smh6$Method <- "SM"
data_smh6$Hive <- "Hive 6"
data_smh6 <- data_smh6[data_smh6$SM_H6 != 0,]
colnames(data_smh6)[2] <- "Abundance metagenomic"

data_smh7$Method <- "SM"
data_smh7$Hive <- "Hive 7"
data_smh7 <- data_smh7[data_smh7$SM_H7 != 0,]
colnames(data_smh7)[2] <- "Abundance metagenomic"

# Combine the data
all_metagenomic_plants <- rbind(data_dsmh4,data_dsmh5,data_dsmh6,data_dsmh7,data_smh4,data_smh5,data_smh6,data_smh7)

# Combine pollen data with metagenomic data
combined_data_meta_pollen <- full_join(all_metagenomic_plants,pollen_analysis, by = c("Family","Hive"))

colnames(combined_data_meta_pollen)[2] <- "Abundance_metagenomic"

# Several options here

# This will select rows that were detected by metagenomics
combined_data_meta_pollen <- combined_data_meta_pollen[combined_data_meta_pollen$Abundance_metagenomic > 0,]

# This will select rows detected by pollen analysis
combined_data_meta_pollen <- combined_data_meta_pollen[combined_data_meta_pollen$Abundance > 0,]

# This will select rows detected by both
combined_data_meta_pollen <- na.omit(combined_data_meta_pollen)

# Only needed if no selection
# combined_data_meta_pollen$Method <- replace_na(combined_data_meta_pollen$Method, "Pollen_analysis_only")

combined_data_meta_pollen$Method <- as.factor(combined_data_meta_pollen$Method)

meta_pollen_scatterplot <- ggscatter(combined_data_meta_pollen, x = "Abundance_metagenomic", y = "Abundance", palette = "jco", add = "reg.line", color = "Method", conf.int = F) +
 # theme_minimal() +
  facet_wrap("Hive", scales = "free_x") +
  theme(plot.title = element_text(size=14), 
        axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12)) +  stat_cor(aes(color=Method)) + labs(x="Metagenomic abundance (%)", y="Pollen analysis abundance (%)")


ggsave(filename = "./Figures/Figure_11/scatterplot_meta_pollen_only_those_detected_by_gyro.pdf")


# Euleur plots
plot_euler_pollen <- function(smplant, dsmplant, pollen) {
  hive <- smplant[1,4]
  sm <- smplant$Family
  dsm <- dsmplant$Family
  pol <- pollen$Family
  smpollen <- ifelse(sum(sm %in% pol) == length(Reduce(intersect, list(dsm, sm, pol))), 0, 
                      ifelse(sum(sm %in% pol) > length(Reduce(intersect, list(dsm, sm, pol))),sum(sm %in% pol)-length(Reduce(intersect, list(dsm, sm, pol))),length(Reduce(intersect, list(dsm, sm, pol)))-sum(sm %in% pol)))
  dsmpollen <- ifelse(sum(dsm %in% pol) == length(Reduce(intersect, list(dsm, sm, pol))), 0, 
                      ifelse(sum(dsm %in% pol) > length(Reduce(intersect, list(dsm, sm, pol))),sum(dsm %in% pol)-length(Reduce(intersect, list(dsm, sm, pol))),length(Reduce(intersect, list(dsm, sm, pol)))-sum(dsm %in% pol)))
  
  plants <- euler(c("DirectSM" = length(setdiff(dsm,sm)), "SM" = length(setdiff(sm,dsm)), "Pollen analysis" = (length(pol)-(length(Reduce(intersect, list(dsm, sm, pol)))+length(sum(sm %in% pol)))),
                "DirectSM&SM" = sum(dsm %in% sm), "DirectSM&Pollen analysis" = dsmpollen, "SM&Pollen analysis" = smpollen, "DirectSM&SM&Pollen analysis" = length(Reduce(intersect, list(dsm, sm, pol)))))
  eulcols <- brewer.pal(3, "Accent")
  
  euplot <- plot(plants,
       quantities = F, fills= eulcols, main = paste(hive), edges = F)

}


hive_4_pollen_plot <- plot_euler_pollen(smplant = data_smh4, dsmplant = data_dsmh4, pollen = hive4_pollen)
hive_5_pollen_plot <- plot_euler_pollen(smplant = data_smh5, dsmplant = data_dsmh5, pollen = hive5_pollen)
hive_6_pollen_plot <- plot_euler_pollen(smplant = data_smh6, dsmplant = data_dsmh6, pollen = hive6_pollen)
hive_7_pollen_plot <- plot_euler_pollen(smplant = data_smh7, dsmplant = data_dsmh7, pollen = hive7_pollen)

ggsave(filename = "hive_4_euler.pdf", plot = hive_4_pollen_plot)
ggsave(filename = "hive_5_euler.pdf", plot = hive_5_pollen_plot)
ggsave(filename = "hive_6_euler.pdf", plot = hive_6_pollen_plot)
ggsave(filename = "hive_7_euler.pdf", plot = hive_7_pollen_plot)


# Previous venn diagrams

venn.diagram(
  x = list(data_dsmh4$Family, hive4_pollen$Family, data_smh4$Family),
  category.names = c("DirectSM", "Pollen analysis" , "SM" ), main = "Hive 4", main.fontfamily = "sans",
  filename = "hive4_venn_pollen_metagenomic_filt_euler.png", output=T,
  print.mode = c("raw","percent"),
  # Circles
  lwd = c("1","1.5","0.5"),
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 0.9,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, -1),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.col = myCol,
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = list(data_dsmh5$Family, hive5_pollen$Family, data_smh5$Family),
  category.names = c("DirectSM", "Pollen analysis" , "SM" ), main = "Hive 5", main.fontfamily = "sans",
  filename = "hive5_venn_pollen_metagenomic_filt.png", output=T,
  print.mode = c("raw","percent"),
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 0.9,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, -1),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.col = myCol,
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = list(data_dsmh6$Family, hive6_pollen$Family, data_smh6$Family),
  category.names = c("DirectSM", "Pollen analysis" , "SM" ), main = "Hive 6", main.fontfamily = "sans",
  filename = "hive6_venn_pollen_metagenomic_filt.png", output=T,
  print.mode = c("raw","percent"),
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 0.9,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, -1),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.col = myCol,
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = list(data_dsmh7$Family, hive7_pollen$Family, data_smh7$Family),
  category.names = c("DirectSM", "Pollen analysis" , "SM" ), main = "Hive 7", main.fontfamily = "sans",
  filename = "hive7_venn_pollen_metagenomic_filt.png", output=T,
  print.mode = c("raw","percent"),
  # Circles
  lwd = 1,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 0.9,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, -1),
  cat.dist = c(0.055, 0.055, 0.055),
  cat.col = myCol,
  cat.fontfamily = "sans",
  rotation = 1
)



# Extract top 3 per method
top3_dsmh4 <- data_dsmh4[order(data_dsmh4$`Abundance metagenomic`,decreasing = T),][1:4,]
top3_dsmh5 <- data_dsmh5[order(data_dsmh5$`Abundance metagenomic`,decreasing = T),][1:4,]
top3_dsmh6 <- data_dsmh6[order(data_dsmh6$`Abundance metagenomic`,decreasing = T),][1:4,]
top3_dsmh7 <- data_dsmh7[order(data_dsmh7$`Abundance metagenomic`,decreasing = T),][1:4,]
top3_smh4 <- data_smh4[order(data_smh4$`Abundance metagenomic`,decreasing = T),][1:4,]
top3_smh5 <- data_smh5[order(data_smh5$`Abundance metagenomic`,decreasing = T),][1:4,]
top3_smh6 <- data_smh6[order(data_smh6$`Abundance metagenomic`,decreasing = T),][1:4,]
top3_smh7 <- data_smh7[order(data_smh7$`Abundance metagenomic`,decreasing = T),][1:4,]

top3_meta <- rbind(top3_dsmh4,top3_dsmh5,top3_dsmh6,top3_dsmh7,top3_smh4,top3_smh5,top3_smh6,top3_smh7)
colnames(top3_meta)[2] <- "Abundance"


top_pollen <- pollen_analysis %>% group_by(Hive) %>% arrange(Abundance, .by_group = TRUE) %>% top_n(4,Abundance) %>% ungroup()
top_pollen$Method <- "Pollen analysis"

top_3_methods <- rbind(top3_meta,top_pollen)

#colnames(top_3_methods)[2] <- "Abundance_metagenomic"

mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(15)


top_most_abundandant_pollen_meta <- ggbarplot(top_3_methods,y= "Abundance",x= "Method",
          fill = "Family", color = "Family", palette = mycolors,
          label = "Family",
          position = position_dodge(0.9), rotate=T, lab.vjust =0.4, lab.hjust = -0.3, legend = "none", main = "Top 4 families per method") + scale_y_continuous(limits = c(0,100))  + facet_grid("Hive") + labs(y="Abundance (%)",x="")


ggsave(filename = "./Figures/Figure_11/most_abundant_pollen_meta_filt.pdf")
