BiocManager::install("Rsamtools")


import_and_count_bam <- function(file.loc) {
  x <- scanBam(file = file.loc)
  qualreads <- x[(x[[1]][["mapq"]] > 40)]
  length(qualreads)
#  reads <- nrow(as.data.frame(x[[1]]["flag"]))
}


species <- c("Apis_mellifera","Artemisia_annua","Asparagus_officinalis","Cajanus_cajan","Helianthus_annuus","Nosema_ceranae","Olea_europaea","Papaver_somniferum","Varroa_destructor","Vitis_vinifera")


full_import <- function(x, sample_name,f) {
  y <- list.files(paste0("./wgs/",x), pattern = "\\.bam$",full.names = T)
  z <- lapply(y, import_and_count_bam)
  a <- as.data.frame(t(as.data.frame(z)))
  rownames(a) <- species
  a <- (a/f)*1000000
  colnames(a) <- sample_name;
  return(a)
}
SMH5 <- full_import("spd3", sample_name = "SMH5", f = 8160732)
SMH4 <- full_import("spd5", sample_name = "SMH4", f = 4058520)
SMH6 <- full_import("spd6", sample_name = "SMH6", f = 3619844)
SMH7 <- full_import("spd7", sample_name = "SMH7", f = 2038988)
DSMH5 <- full_import("spd4", sample_name = "DirectSMH5", f = 1444924)
DSMH4 <- full_import("spd8", sample_name = "DirectSMH4", f = 3282782)
DSMH6 <- full_import("spd9", sample_name = "DirectSMH6", f = 2379244)
DSMH7 <- full_import("spd10", sample_name = "DirectSMH7", f = 2280627)

all_alignments <- cbind(DSMH4,DSMH5,DSMH6,DSMH7,SMH4,SMH5,SMH6,SMH7)
all_alignments <- log2(all_alignments)
all_alignments <- all_alignments %>% rownames_to_column(var="Species")


all_alignments_melt <- reshape2::melt(all_alignments, .id_vars = "Species")

samples <- c("DirectSMH4","DirectSMH5", "DirectSMH6","DirectSMH7","SMH4","SMH5","SMH6","SMH7")


species_level <- read.csv("species_level.csv", header = T, encoding = "UTF-8")
colnames(species_level) <- samples
#species_level <- log2(species_level)
species_level$Species <- species

dsmh4_cpm <- (dsmh4/3282781)*1000000
smh4_cpm <- (smh4/4058520) * 1000000
dsmh5 <- (dsmh5/1444924)*1000000
smh5_cpm <- (smh5/8160732)*1000000
dsmh6_cpm <- (dsmh6/2379244)*1000000
smh6 <- (smh6/3619844)*1000000
dsmh7_cpm <- (dsmh7/2280627)*1000000
smh7_cpm <- (smh7/2038988)*1000000

species_level$DirectSMH4 <- species_level$DirectSMH4/480566
species_level$DirectSMH5 <- species_level$DirectSMH5/154577
species_level$DirectSMH6 <- species_level$DirectSMH6/676415
species_level$DirectSMH7 <- species_level$DirectSMH7/827329
species_level$SMH4 <- species_level$SMH4/730475
species_level$SMH5 <- species_level$SMH5/235515
species_level$SMH6 <- species_level$SMH6/1247102
species_level$SMH7 <- species_level$SMH7/1353429



species_level$DirectSMH4 <- species_level$DirectSMH4/3282781
species_level$DirectSMH5 <- species_level$DirectSMH5/1444924
species_level$DirectSMH6 <- species_level$DirectSMH6/2379244
species_level$DirectSMH7 <- species_level$DirectSMH7/2280627
species_level$SMH4 <- species_level$SMH4/4058520
species_level$SMH5 <- species_level$SMH5/8160732
species_level$SMH6 <- species_level$SMH6/3619844
species_level$SMH7 <- species_level$SMH7/2038988

species_level_melt <- species_level %>% reshape2::melt(., .id_vars="Species")

combined_data <- merge(all_alignments_melt,species_level_melt, by=c("Species","variable"))
colnames(combined_data) <- c("Species","Sample","CPM","DESeq2")

combined_data$Method <- rep(c("DirectSM","DirectSM","DirectSM","DirectSM","SM","SM","SM","SM"), 10)

combined_data <- combined_data[!(combined_data$Species %in% c("Artemisia_annua","Cajanus_cajan","Helianthus_annuus","Nosema_ceranae")),]
combined_data$CPM <- log2(combined_data$CPM)
combined_data$DESeq2 <- log2(combined_data$DESeq2)
combined_data$Species <- str_replace(combined_data$Species, pattern = "_", replacement = " ")

cpm_dsm <- combined_data[combined_data$Method == "DirectSM",]
colnames(cpm_dsm) <- c("Species","Sample","CPM_DSM","DESeq2_DSM","Method")
cpm_dsm <- cpm_dsm[,c(1,3,4)]


cpm_sm <- combined_data[combined_data$Method == "SM",]
colnames(cpm_sm) <- c("Species","Sample","CPM_SM","DESeq2_SM","Method")
cpm_sm <- cpm_sm[,c(1,3,4)]


all_cpm_sm <- cbind(cpm_sm,cpm_dsm)
all_cpm_sm <- all_cpm_sm[,-4]

all_cpm_sm <- all_cpm_sm %>% mutate_all(function(x) ifelse(is.infinite(x), 0, x))  



t <- ggscatter(all_cpm_sm, x = "DESeq2_SM", y = "CPM_SM", palette = "npg",
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "coral",
                                 fill = "lightgray")
) + labs(x="DSM DESeq2 abundance (log2)", y="DSM CPM abundance (log2)", title = "Species level \nCorrelation between alignment DSM and DESeq2 DSM") + theme(legend.position = "right")


s <- facet(t, facet.by = c("Species")) + stat_cor(method= "pearson",label.sep = "\n", label.y = 10)

ggsave(s, filename = "alignment_deseq_dsm_species_level_correlation.pdf")



overall_align_krak <- ggscatter(combined_data, x = "DESeq2", y = "CPM", color = "Species", palette = "jco",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "coral",
                            fill = "lightgray")
)+ stat_cor(method = "pearson", label.y = 13)  + labs(x="kraken2-DESeq2 abundance (log2)", y="CPM (log2)", title = "Genus level \nOverall correlation between alignment and kraken2-DESeq2") + theme(legend.position = "right")
  # Add correlation coefficient

ggsave(overall_align_krak, filename = "overall_correlation_align_kraken_genus.pdf")

# Genus level






combined_data$Species <- str_wrap(combined_data$Species, width = 0, exdent = 2)


ggscatter(combined_data, x = "DESeq2", y = "CPM", palette = "jco", add = "reg.line", shape = "Sample", color = "Species") +
  theme_classic() +
  facet_grid(Method~Species) +
  theme(plot.title = element_text(size=14), 
        axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))+  stat_cor() + labs(x="kraken2-DESeq2 abundance (log2)", y="CPM (log2)")

p <- ggscatter(combined_data, x = "DESeq2", y = "CPM", palette = "jco", add = "reg.line", shape = "Sample", color = "Species") +
     theme_light() + theme(plot.title = element_text(size=14), 
        axis.text=element_text(size=12,colour = "black"),
        axis.title=element_text(size=12),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        panel.grid = element_blank(),
        text = element_text(size=12)) + labs(x="kraken2-DESeq2 abundance (log2)", y="CPM (log2)", title = "Species level \nOverall correlation between alignment and kraken2-DESeq2")

g <- facet(p, facet.by = c("Method","Species")) + stat_cor(method= "pearson",label.sep = "\n", label.y = 12)

ggsave(plot = g, filename = "facet_correlation_align_kraken2_per_species.pdf", width = 12)


gp <- ggplotGrob(g)   
      
for(i in 1:2){
  grob.i <- grep("strip-r", gp$layout$name)[i]
  gp$grobs[[grob.i]]$grobs[[1]]$children[[1]]$gp$fill <- "brown"
}

grid::grid.draw(gp)




















edsmh4_cpm <- (dsmh4/3282781)*1000000
smh4_cpm <- (smh4/4058520) * 1000000
dsmh5 <- (dsmh5/1444924)*1000000
smh5_cpm <- (smh5/8160732)*1000000
dsmh6_cpm <- (dsmh6/2379244)*1000000
smh6 <- (smh6/3619844)*1000000
dsmh7_cpm <- (dsmh7/2280627)*1000000
smh7_cpm <- (smh7/2038988)*1000000


log2(smh4_cpm)





SMH5_files <- list.files("./varroa 2/wgs/spd3", pattern = "\\.bam$",full.names = T)
SMH5_files_import <- lapply(SMH5_files, import_and_count_bam)
SMH5 <- as.data.frame(SMH5_files_import)

rownames(SMH5) <- species



names(all_bams_import) <- species







smh4 <- import_and_count_bam("./varroabam/smh4.bam")
dsmh4 <- import_and_count_bam("./varroabam/dsmh4.bam")
smh5 <- import_and_count_bam("./varroabam/smh5.bam")
dsmh5 <- import_and_count_bam("./varroabam/dsmh5.bam")
smh6 <- import_and_count_bam("./varroabam/smh6.bam")
dsmh6 <- import_and_count_bam("./varroabam/dsmh6.bam")
smh7 <- import_and_count_bam("./varroabam/smh7.bam")
dsmh7 <- import_and_count_bam("./varroabam/dsmh7.bam")

trialbam <- scanBam("./varroabam/smh6.bam")
trialbam_filtered <- trialbam_filtered[(trialbam[[1]][["mapq"]] > 20)] #(trialbam[[1]][["mapq"]] > 20)

length(trialbam_filtered)

plot(dens)


varroa_data1 <- round(as.data.frame(cbind(dsmh4,dsmh5,dsmh6,dsmh7,smh4,smh5,smh6,smh7)))

varroa_data <- do.call("rbind", replicate(5, varroa_data1, simplify = FALSE))




dsmh4_cpm <- (dsmh4/3282781)*1000000
smh4_cpm <- (smh4/4058520) * 1000000
dsmh5 <- (dsmh5/1444924)*1000000
smh5_cpm <- (smh5/8160732)*1000000
dsmh6_cpm <- (dsmh6/2379244)*1000000
smh6 <- (smh6/3619844)*1000000
dsmh7_cpm <- (dsmh7/2280627)*1000000
smh7_cpm <- (smh7/2038988)*1000000
