# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Anastasios Galanis, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figures X with the tree and differences in families
# This script is written in R version 3.6.3

# Load or install (if not present) the required packages
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('purrr')) install.packages('purrr'); library('purrr')
if (!require('taxonomizr')) install.packages('taxonomizr'); library('taxonomizr')
if (!require('here')) install.packages('here'): library('here')
# Taxonomizr will return the taxonomy for each species. However, this requires that a database is built locally (requires 60 GB of space).
# prepareDatabase('nameNode.sqlite')
# This process will take over 3 hours on a regular laptop/PC. Othherwise, please consider dowloading the zipped file 
# from here: and simply unzip it in the cloned repository. Place it at the top level, honeyDSM-seq and not in the subfolders.
if (!require('DESeq2')) install.packages('DESeq2'); library('DESeq2')
if (!require('pheatmap')) install.packages('pheatmap'): library('pheatmap')
if (!require('ggtree')) install.packages('ggtree'): library('ggtree')


data_path <- "./Figures/Figure_3/Data_fig_3"



# Now we will build a figure that shows whether Families are overrepresented in DSM or SM
Hives_normalised_counts <- read.csv("./Figures/Figure_3/normalised_methodseason_family.csv")
Hives_normalised_counts <- Hives_normalised_counts[,-1] %>% column_to_rownames(var = "Taxonomic_ID")

# Remove NAs from Superkingdom
Rounded_counts_with_tax <- Hives_normalised_counts %>% rownames_to_column(var="Taxonomic_ID") %>%  drop_na(Superkingdom)

# Subtract the abundance of DSM and SM to calculate their difference
Family_counts_DSM <- Rounded_counts_with_tax[,c(2:5)]
Family_counts_SM <- Rounded_counts_with_tax[,c(6:9)]
Family_counts_SM$Sum_SM <- rowSums(Family_counts_SM)
Family_counts_DSM$Sum_DSM <- rowSums(Family_counts_DSM)


differences <- Family_counts_DSM$Sum_DSM - Family_counts_SM$Sum_SM
differences <- as.data.frame(differences)

# Exoirt the names so that we can retrieve the tree.
# The csv file is uploaded in https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
# The resulting tree is downloaded as phylip format and then loaded on https://icytree.org
# Then it is downloaded as newick format. 
# The newick tree may need some editing to remove occurrences of this symbol -> ' (yes a single one)
ids <- as.data.frame(Rounded_counts_with_tax$Family)
colnames(ids) <- "ids"
write.table(ids,"./Figures/Figure_3/taxids_fam.csv",row.names = FALSE, col.names = FALSE)


# Import the tree
tree1 <- read.tree("./Figures/Figure_3/treefam.newick")
tree1

# generate a circular tree
circ <- ggtree(tree1, branch.length = "none", layout = "circular", ladderize = TRUE) + geom_tiplab(size=1.2)
circ

# Create a dataframe with node numbers and their labels so that we can colour the domains later
node <- circ[["data"]][["node"]]
label <- circ[["data"]][["label"]]
node_table <- cbind(node,label)

# We create a dataframe that colours all nodes grey. Custom colour of gheatmap is black and I don't like it.
colouring <- data.frame(nodes=label,color = "grey")

# Attach the data to the tree
circ2 <- circ %<+% colouring + aes(color = I(color))
circ2

# Import the log2fc data
family_lfc <- read.csv("./Output_data/Figure_3_output/family_lfc.csv")[,-1]
direction <- as.data.frame(direction)
colnames(direction) <- c("Direction")
family_lfc$Family <- as.character(family_lfc$Family)
family_lfc$Family <- str_remove(family_lfc$Family, pattern = " ")


# Extract tip label names so we can reorder
reordering <- as.data.frame(tree1$tip.label)
colnames(reordering) <- c("treetip")
reordering$treetip <- as.character(reordering$treetip)

# For some reason it won't accept me to do it any other way so it has to be done twice..
reordering$treetip <- str_remove(reordering$treetip, "\"")
reordering$treetip <- str_remove(reordering$treetip, "\"")

# Reordered dataframe with data
family_lfc2 <- family_lfc[match(reordering$treetip, family_lfc$Family),]

rownames(family_lfc) <- family_lfc$Family

family_lfc_plot <- as.data.frame(family_lfc$log2FoldChange)
colnames(family_lfc_plot) <- c("Log2FoldChange")
rownames(family_lfc_plot) <- tree1$tip.label

# Extract the direction of the log2fc. For the data I have used the contrast DirectSM vs SM.
# Thereofre, positive values mean SM > DSM and negative values DSM > SM
direction <- as.data.frame(ifelse(family_lfc$log2FoldChange < 0, "Higher DirectSM","Higher SM"))
colnames(direction) <- c("Direction")
rownames(direction) <- tree1$tip.label

# Get the adjusted pvalues
padj_values_family <- as.data.frame(family_lfc$padj)
colnames(padj_values_family) <- c("padjvalue")
rownames(padj_values_family) <- tree1$tip.label

# We can now build the heatmap
p2 <- gheatmap(circ2, padj_values_family, colnames_angle=0, offset = 0, width = 0.2, colnames = F) + scale_fill_distiller(palette = "RdGy", limits=c(0,1), direction = 1, breaks = c(0, 1)) + theme(legend.position = "bottom", legend.title.align = 0.5) + guides(fill = guide_colourbar(barwidth = 5, ticks = FALSE, title = "p-value (adjusted)", title.position = "top"))


# Add annotation to the Domains as we've done so far 
p3 <- p2  + geom_cladelabel(node = 486, color ="#7570B3", align = T, label = "", barsize = 3) + geom_cladelabel(node = 543, color = "#1B9E77", align = T, label = "", barsize = 3, offset = 0) + geom_cladelabel(node = 578, color = "#E6AB02", align = T, label = "", barsize = 3) + geom_cladelabel(node = 479, color = "#eb53a1", align = T, label = "", barsize = 3) + geom_cladelabel(node = 575, color = "#D95F02", align = T, label = "", barsize = 3)

p4 <- p3 + new_scale_fill()
p5 <- gheatmap(p4, family_lfc_plot, colnames_angle=0, offset = 1, width = 0.15, colnames = F) + scale_fill_distiller(palette = "PiYG", limits = c(-2,2), labels = c("DirectSM","","SM"), breaks = c(-2,0,2)) + theme(legend.position = "bottom", legend.title.align = 0.5) + guides(fill = guide_colourbar(barwidth = 5, ticks = FALSE, title = "Log2FC", title.position = "top"))

ggsave(p5, filename="circlularplot_families.pdf", device = "pdf", width = 15, height = 15, path = "./Figures/Figure_3/")


family_lfc_replace_plant <- family_lfc
family_lfc_replace_plant$Kingdom <- as.character(family_lfc_replace_plant$Kingdom)
family_lfc_replace_plant$Kingdom <- replace_na(family_lfc_replace_plant$Kingdom, "Unknown")
family_lfc_replace_plant$Superkingdom <- ifelse(family_lfc_replace_plant$Kingdom == "Viridiplantae", 
                                                as.character(family_lfc_replace_plant$Kingdom), as.character(family_lfc_replace_plant$Superkingdom))

family_lfc_threshold <- subset(family_lfc_replace_plant, abs(family_lfc_replace_plant$log2FoldChange) >= 1)
family_lfc_threshold_total <- family_lfc_threshold %>% group_by(Superkingdom) %>% summarise(Number = n())
family_lfc_threshold_total <- as.data.frame(family_lfc_threshold_total)

lfc_bars <- ggplot(family_lfc_threshold_total, aes(x=Superkingdom, y= Number, fill=Superkingdom)) + geom_bar(stat="identity") + coord_flip() + scale_fill_manual(values = c("#E6AB02","#7570B3","#1B9E77","#eb53a1"))+ theme_minimal() + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.text.y = element_text(size=15), plot.title = element_text(size=18)) + geom_text(aes(label=Number), vjust = 0, hjust=-2, size = 10) + ylim(c(0,11)) + labs(title="Distribution of most variable Families", x = "")

ggsave(plot = lfc_bars, filename = "./Figures/Figure_3/distribution_families.pdf")
