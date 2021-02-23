# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figures X with the environmental data
# This script is written in R version 3.6.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('here')) install.packages('here'): library('here')
setwd(here::here())

weather <- read.csv("./Figures/Figure_5/weather.csv")

weather$DAY <- as.Date(weather$DAY, format = "%d/%m/%Y")

temperature_plot <- ggplot(weather, aes(x = DAY, y = MEAN..TEMP)) + geom_line() + geom_ribbon(aes(x = DAY, ymax = HIGH, ymin = LOW), alpha = 0.6, fill = "tan2") + labs(x = "", y = "Temperature (°C)") + theme_bw() + geom_vline(xintercept = weather$DAY[121], size = 1.5, linetype = "dotted", color = "black") + geom_vline(xintercept = weather$DAY[182], size = 1.5, linetype = "dotted", color = "black") + geom_vline(xintercept = weather$DAY[305], size = 1.5, linetype = "dotted", color = "black") + annotate("text", x =  weather$DAY[111], y = 3, label = "Hive 5", angle = 90) + annotate("text", x =  weather$DAY[172], y = 3, label = "Hive 6 and 7", angle = 90) + annotate("text", x =  weather$DAY[295], y = 3, label = "Hive 4", angle = 90)

ggsave(plot=temperature_plot, filename="temperature_plot.pdf",path="./Figures/Figure_5/", width = 9)

rain_plot <- ggplot(weather, aes(x = DAY, y = RAIN)) + geom_line(color="steelblue1") + theme_bw() + labs(x = "", y = "") + geom_vline(xintercept = weather$DAY[121], size = 1.5, linetype = "dotted", color = "black") + geom_vline(xintercept = weather$DAY[182], size = 1.5, linetype = "dotted", color = "black") + geom_vline(xintercept = weather$DAY[305], size = 1.5, linetype = "dotted", color = "black") + annotate("text", x =  weather$DAY[111], y = 60, label = "Hive 5", angle = 90) + annotate("text", x =  weather$DAY[172], y = 60, label = "Hive 6 and 7", angle = 90) + annotate("text", x =  weather$DAY[295], y = 60, label = "Hive 4", angle = 90) + labs(y="Rainfall (mm)")
ggsave(plot=rain_plot, filename="rain_plot.pdf",path="./Figures/Figure_5/", width = 9)


wind_plot <- ggplot(weather, aes(x = DAY, y = AVG.WIIND.SPEED)) + geom_line(color="orange4") + theme_bw() + labs(x = "", y = "Speed (km/h)") + geom_vline(xintercept = weather$DAY[121], size = 1.5, linetype = "dotted", color = "black") + geom_vline(xintercept = weather$DAY[182], size = 1.5, linetype = "dotted", color = "black") + geom_vline(xintercept = weather$DAY[305], size = 1.5, linetype = "dotted", color = "black")  + annotate("text", x =  weather$DAY[111], y = 25, label = "Hive 5", angle = 90) + annotate("text", x =  weather$DAY[172], y = 25, label = "Hive 6 and 7", angle = 90) + annotate("text", x =  weather$DAY[295], y = 25, label = "Hive 4", angle = 90)
ggsave(plot=wind_plot, filename="wind_plot.pdf",path="./Figures/Figure_5/", width = 9)