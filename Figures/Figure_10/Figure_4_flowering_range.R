# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figure 4 (flowering range)
# This script is written in R version 4.0.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('here')) install.packages('here'): library('here')

if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')


flowering_period <- read.csv("./Figures/Figure_10/plant_flowering_table.csv")

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(18)

flowering_plot <- ggplot(flowering_period, aes(x=flowering_period$ï..Identified_species_metagenomic, color=ï..Identified_species_metagenomic))+
  geom_linerange(aes(ymin=Month_start,ymax=Month_end),linetype=1, size=2, alpha=0.7)+
  geom_point(aes(y=Month_start),size=4)+
  geom_point(aes(y=Month_end),size=4)+ theme_bw() + theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 0.5), text = element_text(size=12)) + 
   scale_color_manual(values = mycolors) + coord_flip() + scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), limits = c(1,12), labels =  c("January","February","March","April","May","June","July","August","September","October","November","December")) + labs(x="",y="Flowering range")


ggsave(plot = flowering_plot,filename = "./Figures/Figure_10/flowering_plot_plants.pdf", height = 8)

