# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Supplementary Figure 9 A
# This script is written in R version 4.0.3

# Load or install (if not present) the required packages
if (!require('maps')) install.packages('maps'); library('maps')
if (!require('sf')) install.packages('sf'); library('sf')
if (!require('rnaturalearth')) install.packages('rnaturalearth'); library('rnaturalearth')
if (!require('here')) install.packages('here'): library('here')
if (!require('rnaturalearthdata')) install.packages('rnaturalearthdata'); library('rnaturalearthdata')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggspatial')) install.packages('ggspatial'); library('ggspatial')
if (!require('rgeos')) install.packages('rgeos'); library('rgeos')


# load data
world <- ne_countries(scale = "medium", returnclass = "sf")
# gene world map
worldmaps <- ggplot(data = world) +
  geom_sf(fill="grey90")+ theme(panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25),panel.background = element_rect(fill = "aliceblue")) + 
  ggtitle("Sampled areas", subtitle = "BIEN database search") + annotate("rect",size=1.3,xmin=-9,xmax=36.1,ymin=33,ymax=44, colour="darkorange", alpha=0) + annotate("rect",size=1.3,xmin=90,xmax=135.1,ymin=33,ymax=44, colour="cyan3", alpha=0) + annotate("rect",size=1.3,xmin=-121,xmax=-75.9,ymin=33,ymax=44, colour="deeppink", alpha=0)

mediterranean <- ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-9, 36.1), ylim = c(33, 44), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw() + theme(panel.border = element_rect(colour="darkorange", size=2),panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25),panel.background = element_rect(fill = "aliceblue"))

asia <- ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(90, 135.1), ylim = c(33, 44), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw() + theme(panel.border = element_rect(colour="cyan3", size=2),panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25),panel.background = element_rect(fill = "aliceblue"))

northamerica <- ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitude", y = "Latitude") +
  coord_sf(xlim = c(-121,-75.9), ylim = c(33, 44), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw() + theme(panel.border = element_rect(colour="deeppink", size=2),panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25),panel.background = element_rect(fill = "aliceblue"))


ggsave(filename = "med_map.pdf", plot = mediterranean)
ggsave(filename = "asia_map.pdf",plot = asia)
ggsave(filename = "usa_map.pdf",plot = northamerica)
ggsave(filename = "world_map.pdf",plot = worldmaps)
