# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figures 3A, 3B (30 most abundant species) and 3C but will also create similar Figures to 3B for Family and Genus levels
# This script is written in R version 3.6.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggtree')) install.packages('ggtree'); library('ggtree')
if (!require('here')) install.packages('here'): library('here')
setwd(here::here())


simulations_tree <- read.nexus("./Figures/Figure_5/Data_fig_5/simulations_tree_data.nexus")


annotated_tree_simulations <- ggtree(simulations_tree) + geom_cladelabel(node=67, label="Eukaryota", align=T, fontsize=8, offset = 16) + geom_cladelabel(node=87, label="Bacteria", align=T, fontsize=8, offset = 16) + geom_cladelabel(node=66, label="Viruses", align=T, fontsize=8, offset = 16) + geom_tiplab(align = TRUE, offset = 1, label = gsub("_", " ", simulations_tree@phylo[["tip.label"]]), fontface="italic") + xlim(c(0,50)) + geom_balance(node=76, fill='green', color='white', alpha=0.3, extend=1) + geom_hilight(node=87, fill="brown", alpha = 0.2, extend = 5) + geom_hilight(node=67, fill="lightblue", alpha = 0.3, extend = 1) + geom_hilight(node=66, fill="violet", alpha = 0.2, extend = 21)

ggsave(plot=annotated_tree_simulations,filename="annotatedtree.pdf", height = 15, width = 11)


organisms_simulation <- read.csv("./Figures/Figure_5/Data_fig_5/organisms.csv")
organisms_simulation$Genome_size <- as.numeric(levels(organisms_simulation$Genome_size))[organisms_simulation$Genome_size]
neworganisms_simulation <- organisms_simulation[-1]
neworganisms_simulation$Genome_size[46] <- 22103.6
neworganisms_simulation$Genome_size <- log10(neworganisms_simulation$Genome_size)
neworganisms_simulation$Species <- gsub(" ", "_", neworganisms_simulation$Species, fixed = TRUE)
neworganisms_simulation <- neworganisms_simulation[,c(1:3)]
neworganisms_simulation$Species <- simulations_tree$tip.label

facet_tree_simulation <- ggtree(simulations_tree)

first_plot <- facet_plot(facet_tree_simulation, panel = "Genome size (Mb; log10)", data = neworganisms_simulation, mapping = aes(x=Genome_size), geom = geom_point) + xlim_expand(c(-3, 5), panel = "Genome size (Mb; log10)") + geom_tiplab(align = TRUE, parse = FALSE, size=3, label = gsub("_", " ", simulations_tree[["tip.label"]]), fontface="italic") + xlim_tree(c(0,350)) + theme_tree2() 
final_facet_tree_simulation <- first_plot + geom_facet(panel = "GC content (%)", data = neworganisms_simulation, aes(x=GC_count), geom = geom_point)

facet_widths(final_facet_tree_simulation, c(Tree = 0.5))
ggsave(plot=final_facet_tree_simulation, filename="tree_plot_genome_size.pdf", height = 10, width = 8)
