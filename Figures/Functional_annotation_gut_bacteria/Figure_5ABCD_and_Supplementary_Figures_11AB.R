# Direct shotgun metagenomics captures species abundance of honey samples
# R script author: Solenn Patalano, 2020
# This script will INSTALL packages that have not been previously installed. 
# If you do not wish to install new packages, please avoid running this script.

# This script will reproduce Figures 5 A, B, C, D and Supplementary Figures 11 A and B
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


myPalette <- brewer.pal(6, "Dark2")

#All non plant species 
nonplants<-read.csv("./Normalised_reads/normalised_methodseason_species.csv", header = TRUE)
nonplants<-nonplants[!grepl("Viridiplantae", nonplants$Kingdom),]

#filter species with less that 100 reads and export for manual and NCBI-based classification
nonplants$Sum<-rowSums(nonplants[3:10])
nonplantsfilt<-nonplants[!(rowSums(nonplants[3:10])<=100), ] # this threshold give the same distribution that with all read
write.table(nonplantsfilt, file = "nonplantsfilt.txt", sep = "\t",row.names = TRUE, col.names = NA)

#merge with Pavlopoulos analysis (NCBI search)
NCBI<-read_excel("NCBI.xlsx")
NCBI<-NCBI[!grepl("Apis mellifera filamentous virus", NCBI$Species),]

#Differential non plant species across seasons (from Cluster analysis)
SeasonalDifCl1<-read.csv("./Figures/Figure_3/Cluster/Species_cluster_1.csv", header = TRUE)
SeasonalDifCl3<-read.csv("./Figures/Figure_3/Cluster/Species_cluster_3.csv", header = TRUE)
SeasonalDifCl4<-read.csv("./Figures/Figure_3/Cluster/Species_cluster_4.csv", header = TRUE)
SeasonalDifAll<-rbind(SeasonalDifCl1,SeasonalDifCl3,SeasonalDifCl4)
SeasonalDifAll<-SeasonalDifAll[!grepl("Viridiplantae", SeasonalDifAll$Kingdom),]
SeasonalDifAll<-SeasonalDifAll[,-(2:11)]
SeasonalDifAll<-merge(SeasonalDifAll,NCBI,by="Species")


#plot distribution
All<-nonplants[!grepl("Apis mellifera filamentous virus", nonplants$Species),] #remove Am f virus
All<-data.frame(tapply(All$Sum, All$Superkingdom, FUN=sum))
pie (t(All), labels = c("Archaea","Bacteria","Eucaryote","Virus"))

Filter<-data.frame(tapply(NCBI$Sum, NCBI$Superkingdom, FUN=sum))
pie (t(Filter), labels = c("Bacteria","Eucaryote","NA","Virus"),col=myPalette,main="Non-plants superkindgoms distribution in honey")

Relationship <- data.frame(tapply(NCBI$Sum, NCBI$`Relation category`  , FUN=sum))
pie (t(Relationship), labels = c("Bacterial gut community","Host","Human cross contamination","Others","Pathogen","Unknown"),border="white", col=myPalette, main="Relationship of non-plant species with bees ")

Differential <- data.frame(tapply(SeasonalDifAll$Sum, SeasonalDifAll$`Relation category`  , FUN=sum))
pie (t(Differential), labels = c("Bacterial gut community","Host","Human cross contamination","Others","Pathogen","Unknown"),border="white", col=myPalette, main="Relationship of the seasonal specific non-plant species with bees ")



#Extract data of the known microbiota (from Kesnerová et al. 2019 and doi: 10.1186/s12864-015-1476-6)

Bartonella_sp <- nonplants[grepl("Bartonella", nonplants$Genus), ]
Bartonella_sp$Category<- c("Bartonella sp.")

Bifidobacterium_sp <- nonplants[grepl("Bifidobacterium", nonplants$Genus), ]
Bifidobacterium_sp$Category<- c("Bifidobacterium sp.")

Bombella_sp <-  nonplants[grepl("Bombella", nonplants$Genus), ]
Bombella_sp$Category<- c("Bombella sp.")

Frischellaperrara <- nonplants[grepl("Frischella perrara", nonplants$Species), ]
Frischellaperrara$Category<- c("Frischella perrara")

Gilliamellaapicola <- nonplants[grepl("Gilliamella apicola", nonplants$Species), ]
Gilliamellaapicola$Category<- c("Gilliamella apicola")

Lactobacilluskunkeei <- nonplants[grepl("Lactobacillus kunkeei", nonplants$Species), ]
Lactobacilluskunkeei$Category<- c("Lactobacillus kunkeei")

Lactobacillus_Firm4 <- nonplants[grepl("Lactobacillus mellis", nonplants$Species) | grepl("Lactobacillus mellifer", nonplants$Species) , ] #doi: 10.1186/s12864-015-1476-6
Lactobacillus_Firm4$Category<- c("Lactobacillus-Firm4") #not plot because of too lower reads

Lactobacillus_Firm5 <- nonplants[grepl("Lactobacillus apis", nonplants$Species) | grepl("Lactobacillus helsingborgensis", nonplants$Species)| grepl("Lactobacillus melliventris", nonplants$Species) | grepl("Lactobacillus kullabergensis", nonplants$Species)| grepl("Lactobacillus kimbladii", nonplants$Species) , ] #doi: 10.1186/s12864-015-1476-6
Lactobacillus_Firm5$Category<- c("Lactobacillus-Firm5")

Lactobacillus_Others <- nonplants[grepl("Lactobacillus sp. Fhon2N", nonplants$Species) | grepl("Lactobacillus apinorum", nonplants$Species)| grepl("Lactobacillus delbrueckii", nonplants$Species) | grepl("Liquorilactobacillus nageliis", nonplants$Species), ]  #Others Lactobacillus with more that 100 reads coverage
Lactobacillus_Others$Category<- c("Lactobacillus-Others")

Lonsdaleabritannica <- nonplants[grepl("Lonsdalea britannica", nonplants$Species), ]
Lonsdaleabritannica$Category<- c("Lonsdalea britannica")

Snodgrassellaalvi <- nonplants[grepl("Snodgrassella alvi", nonplants$Species), ]
Snodgrassellaalvi $Category<- c("Snodgrassella alvi ") #not plot because of too lower reads

microbiota<-rbind(Bartonella_sp,Bifidobacterium_sp,Bombella_sp,Frischellaperrara,Gilliamellaapicola,Lactobacilluskunkeei,Lactobacillus_Firm5,Lactobacillus_Others,Lonsdaleabritannica) 
#Because L. kunkeii dominate the DNA from microbiota, I remove it to clarify the barplots
microbiota_nokunkeei<-rbind(Bartonella_sp,Bifidobacterium_sp,Bombella_sp,Frischellaperrara,Gilliamellaapicola,Lactobacillus_Firm5,Lactobacillus_Others,Lonsdaleabritannica) 


#Abundance barplots
### ALL
p<-microbiota %>% 
  group_by(Category) %>% 
  summarise_at(vars(DirectSM_H5,SM_H5,DirectSM_H7,SM_H7,DirectSM_H6,SM_H6,DirectSM_H4,SM_H4  ), funs(sum))
pp <- melt(p, id.vars="Category")
ppp <-pp %>% 
  group_by(Category) %>% 
  summarise(mean=mean(value), sd=sd(value))

ggplot(ppp, aes(reorder(Category,- mean), y=mean, fill=Category)) + 
  geom_bar(stat="identity", color="black", position=position_dodge(), fill="grey") +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),position = position_dodge(0.9), width = .3) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45,hjust=1), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  xlab("Species") +
  ylab("Attributed reads") 


#### NO KUNKEII
pk<-microbiota_nokunkeei %>% 
  group_by(Category) %>% 
  summarise_at(vars(DirectSM_H5,SM_H5,DirectSM_H7,SM_H7,DirectSM_H6,SM_H6,DirectSM_H4,SM_H4  ), funs(sum))
ppk <- melt(pk, id.vars="Category")
ggplot(ppk, aes(variable, value,fill=Category)) + 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_brewer(palette="Paired")+
  labs(title = "CORE and NON-CORE microbiota detected in Honey",
       subtitle = "(Excluding L. kunkeii)")



#Because others abundant bacteria were attributed as gut microbiota from the NCBI file and the differential analyses I decide to include them also in the microbiota database

Sodalis <- nonplants[grepl("Sodalis glossinidius", nonplants$Species) | grepl("Sodalis praecaptivus", nonplants$Species) , ] 
Sodalis$Category<- c("Sodalis")

Sodalisglossinidius <- nonplants[grepl("Sodalis glossinidius", nonplants$Species)  , ] 
Sodalisglossinidius$Category<- c("Sodalis glossinidius")

Sodalispraecaptivus <- nonplants[grepl("Sodalis praecaptivus", nonplants$Species)  , ] 
Sodalispraecaptivus$Category<- c("Sodalis praecaptivus")

Pantoeaagglomerans<-nonplants[grepl("Pantoea agglomerans", nonplants$Species) , ]
Pantoeaagglomerans$Category<- c("Pantoea agglomerans")

Leuconostocpseudomesenteroides<-nonplants[grepl("Leuconostoc pseudomesenteroides", nonplants$Species) , ]
Leuconostocpseudomesenteroides$Category<- c("Leuconostoc pseudomesenteroides")

Enterobacter_sp<-nonplants[grepl("Enterobacter sp. SA187", nonplants$Species) , ]
Enterobacter_sp$Category<- c("Enterobacter sp.")

Arsenophonusnasoniae<-nonplants[grepl("Arsenophonus nasoniae", nonplants$Species) , ]
Arsenophonusnasoniae$Category<- c("Arsenophonus nasoniae")

Parasaccharibacterapium<-nonplants[grepl("Parasaccharibacter apium", nonplants$Species) , ]
Parasaccharibacterapium$Category<- c("Parasaccharibacter apium")

Morganellamorganii<-nonplants[grepl("Morganella morganii", nonplants$Species) , ]
Morganellamorganii$Category<- c("Morganella morganii")

Klebsiellaoxytoca<-nonplants[grepl("Klebsiella oxytoca", nonplants$Species) , ]
Klebsiellaoxytoca$Category<- c("Klebsiella oxytoca")

microbiota_extended<-rbind(microbiota_nokunkeei,Sodalis, Pantoeaagglomerans,Leuconostocpseudomesenteroides,Enterobacter_sp,Arsenophonusnasoniae,Parasaccharibacterapium,Morganellamorganii,Klebsiellaoxytoca)
other_microbiota<-rbind(Sodalis,Pantoeaagglomerans,Leuconostocpseudomesenteroides,Enterobacter_sp,Arsenophonusnasoniae,Parasaccharibacterapium,Morganellamorganii,Klebsiellaoxytoca)

#### EXTENDED
pkext<-microbiota_extended %>% 
  group_by(Category) %>% 
  summarise_at(vars(DirectSM_H5,SM_H5,DirectSM_H7,SM_H7,DirectSM_H6,SM_H6,DirectSM_H4,SM_H4  ), funs(sum))
ppkext <- melt(pkext, id.vars="Category")
ggplot(ppkext, aes(variable, value,fill=Category)) + 
  geom_bar(position="fill", stat="identity")+
  theme_minimal() +
  labs(title = "EXTENTED GUT MICROBIOTA")

#### OTHERS MICROBIOTA
pkother<-other_microbiota %>% 
  group_by(Category) %>% 
  summarise_at(vars(DirectSM_H5,SM_H5,DirectSM_H7,SM_H7,DirectSM_H6,SM_H6,DirectSM_H4,SM_H4  ), funs(sum))
ppkother <- melt(pkother, id.vars="Category")
ggplot(ppkother, aes(variable, value,fill=Category)) + 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  scale_fill_brewer(palette="Paired")+
  labs(title = "OTHER MICROBIOTA")

### SEASONAL DIFF MICROBIOTA (all except frishella are autonm specific so I removed it for plotting purpose  )
Diff_microbiota<-rbind(Enterobacter_sp,Klebsiellaoxytoca,Leuconostocpseudomesenteroides,Sodalisglossinidius,Sodalispraecaptivus)#,Frischellaperrara
pkdiff<-Diff_microbiota %>% 
  group_by(Category) %>% 
  summarise_at(vars(DirectSM_H5,SM_H5,DirectSM_H7,SM_H7,DirectSM_H6,SM_H6,DirectSM_H4,SM_H4  ), funs(sum))
ppkdiff <- melt(pkdiff, id.vars="Category")
ggplot(ppkdiff, aes(variable, value,fill=Category)) + 
  geom_bar(position="stack", stat="identity")+
  theme_minimal()+
  scale_fill_brewer(palette="Oranges")+
  facet_grid(Category~variable, scales="free", space="free_x")

functional<-read.csv("functional_results_significant_bacteria.csv", header = TRUE)
Func_Sum<-functional %>% 
  group_by(Species) %>% 
  summarise_at(vars(DirectSM_H5,SM_H5,DirectSM_H7,SM_H7,DirectSM_H6,SM_H6,DirectSM_H4,SM_H4  ), funs(sum))

Functional_Reduced<- functional[grepl("Enterobacter sp. SA187", functional$Species) | grepl("Klebsiella oxytoca", functional$Species) | grepl("Leuconostoc pseudomesenteroides", functional$Species)| grepl("Sodalis glossinidius", functional$Species)| grepl("Sodalis praecaptivus", functional$Species) , ] 

Functional_Reduced<-Functional_Reduced[,-c(1:2,4:6,8:11,14)]

Func<- reshape2::melt(Functional_Reduced, id=c("GO_description","Species"))
colnames(Func)[3] <- "Method"
colnames(Func)[4] <- "Frequency"


Func_10 <- Func%>%
  mutate_if(sapply(Func, is.character), as.factor)%>%
  filter(!GO_description=="NO_NAME")%>%
  group_by(GO_description) %>%
  filter(n() == 10) #only GO terms present in the 5 species and both methods

summary(Func_10)

b<-Func_10 %>% 
  nest(data = c(Species, Method,Frequency)) %>% 
  mutate(model = map(data, ~anova(lm(Frequency ~ Species, .))), 
         tidy = map(model, broom::tidy)) %>% 
  select(GO_description, tidy) %>% 
  unnest(tidy)

b$significance<-ifelse(b$p.value<0.005,"significant","not sign.")  

#################################

#In order to perform statistical test across ALL GO, I need to attribute a zero value to the missing species
fill<-reshape2::dcast(Func, GO_description + Method ~ Species,value.var = "Frequency", fun.aggregate = mean)
fill[is.na(fill)] <- 0 #The missing value for certain GO corresponds to true zeros (no hits during alignments)

Func_all<- reshape2::melt(fill, id=c("GO_description","Method"))
colnames(Func_all)[3] <- "Species"
colnames(Func_all)[4] <- "Frequency"


Func_all_m <- Func_all %>%
  mutate_if(sapply(Func_all, is.character), as.factor)%>%
  filter(!GO_description=="NO_NAME")%>%
  group_by(GO_description) %>%
  filter(n() == 10) 

summary(Func_all_m)

#cube root transformation
#Func_all_m$Frequency<-sign(Func_all_m$Frequency) * abs(Func_all_m$Frequency)^(1/3)
#plotNormalHistogram(Func_all_m$Frequency)

#####CHOOSE#######

###### ANOVA
c<-Func_all_m %>% 
  nest(data = c(Species, Method,Frequency)) %>% 
  mutate(model = map(data, ~anova(lm(Frequency ~ Species, .))), 
         tidy = map(model, broom::tidy)) %>% 
  select(GO_description, tidy) %>% 
  unnest(tidy)

###### NON PARAMETRIC
#c<-Func_all_m %>% 
#nest(data = c(Species, Method,Frequency)) %>% 
#mutate(model = map(data, ~kruskal.test(Frequency ~ Species, .)), 
#tidy = map(model, broom::tidy)) %>% 
#select(GO_description, tidy) %>% 
#unnest(tidy)

#####CHOOSE#######

c$significance<-ifelse(c$p.value<0.005,"significant","not sign.")  


###### DATA CHECKS ######
anova <- aov(Frequency ~ Species , data = Func_all_m)
# Test for Homogeneity of variances (OK IF around red line)
plot(anova, 1)
# Second Test for Homogeneity of variances (OK if p value significatif)
leveneTest(Frequency ~ Species, data = Func_all_m)
# Test for Normality distribution (OK If most point follow the line)
plot(anova, 2)
# Second test for normality Extract the residuals
anova_residuals <- residuals(object = anova )
# Run Shapiro-Wilk test (OK if p value significatif)
#shapiro.test(x = anova_residuals ) #This test can't be applied to so many data


#### FINAL STAT, CLEAN and EXPORT TABLE ####

significant<-c%>% 
  filter(!term=="Residuals")

stat<-Func_all_m %>%
  group_by(GO_description, Species) %>%
  summarise(
    n=n(),
    mean = mean(Frequency, na.rm = TRUE),
    sd = sd(Frequency, na.rm = TRUE)
  )

stat1<-reshape2::dcast(stat, GO_description ~ Species,value.var = "mean")
colnames(stat1)[2] <- "Enterobacter sp. SA187 - Mean SM/DSM"
colnames(stat1)[3] <- "Klebsiella oxytoca - Mean SM/DSM"
colnames(stat1)[4] <- "Leuconostoc pseudomesenteroides - Mean SM/DSM"
colnames(stat1)[5] <- "Sodalis glossinidius - Mean SM/DSM"
colnames(stat1)[6] <- "Sodalis praecaptivus - Mean SM/DSM"
stat2<-reshape2::dcast(stat, GO_description ~ Species,value.var = "sd")
colnames(stat2)[2] <- "Enterobacter sp. SA187 - SD SM/DSM"
colnames(stat2)[3] <- "Klebsiella oxytoca - SD SM/DSM"
colnames(stat2)[4] <- "Leuconostoc pseudomesenteroides - SD SM/DSM"
colnames(stat2)[5] <- "Sodalis glossinidius - SD SM/DSM"
colnames(stat2)[6] <- "Sodalis praecaptivus - SD SM/DSM"

GO_sign<- merge(stat1,stat2,by="GO_description")
GO_sign<-merge(GO_sign,significant,by="GO_description")
write.table(GO_sign, "GO_significant.txt", sep="\t")


# Export the columns 
GO_sign_columns <- GO_sign[,2:6]

# Rename the columns
colnames(GO_sign_columns) <- c("Enterobacter sp. SA187","Klebsiella oxytoca","Leuconostoc pseudomesenteroides","Sodalis glossinidius","Sodalis praecaptivus")

# Extract the maximum number (or Frequency) per GO (this is a choice as we assumed that the highest mean is bringing the significance)
maximum_column_names <- colnames(GO_sign_columns)[max.col(GO_sign_columns,ties.method = "first")] %>% as.data.frame()

# Combine it with the pvalue data
maximums_pvalues_and_names <- cbind(maximum_column_names, GO_sign$GO_description,GO_sign$p.value)
colnames(maximums_pvalues_and_names) <- c("Species","GO_description","pvalue")
maximums_pvalues_and_names <- maximums_pvalues_and_names %>% group_by(Species)

# Create a list with the GOs per species. Each dataframe contains the data for each species. 
maximum_split <- maximums_pvalues_and_names %>% group_split()
top20_per_species <- lapply(maximum_split, function(df) {
  df <- df[order(df$pvalue),]
  df <- df[1:20,]
})

# Manually annotate the functions (copy paste to the GO database)
enterobacter <- top20_per_species[[1]]
enterobacter$functions <- c("MF","MF","MF","MF","MF","MF","MF","MF","CC","MF","CC","MF","BP","BP","MF","BP","MF","MF","MF","BP")
enterobacter$Topfunctions<-c( "Carbohydrate process", "Others","Others","Others","Others","Carbohydrate process","Others","Others","Others","Carbohydrate process","Others","Carbohydrate process","Others","Others","Others","Others","Others","Others","Others","Carbohydrate process" )
enterobacter$pvalue<--log10(enterobacter$pvalue)
enterobacter<-merge(enterobacter,GO_sign,by=c("GO_description"))
enterobacter<-enterobacter[,-c(7:10,12:22)]
colnames(enterobacter)[6] <- "Frequency"
colnames(enterobacter)[7] <- "SD"

klebsiella <- top20_per_species[[2]]
klebsiella$functions <- c("MF","MF","MF","BP","MF","MF","MF","MF","MF","MF","BP","BP","BP","MF","MF","MF","MF","MF","CC","MF")
klebsiella$Topfunctions<-c("Others","C-C lyase activity","tRNA modif","tRNA modif","Phosphorilation","C-C lyase activity","C-C lyase activity","C-C lyase activity","Phosphorilation","Others","Phosphorilation","Others","Others","Others","Others","Others","tRNA modif","Others","Others","Others")
klebsiella$pvalue<--log10(klebsiella$pvalue)
klebsiella<-merge(klebsiella,GO_sign,by=c("GO_description"))
klebsiella<-klebsiella[,-c(6,8:11,13:22)]
colnames(klebsiella)[6] <- "Frequency"
colnames(klebsiella)[7] <- "SD"

leuconostoc <- top20_per_species[[3]]
leuconostoc$functions <- c("MF","BP","BP","MF","CC","MF","CC","BP","BP","MF","MF","MF","MF","MF","BP","MF","MF","MF","BP","BP")
leuconostoc$Topfunctions<-c("tRNA modif","tRNA modif","DNA repair","Glycolysis","Others","Glycolysis","Others","Others","Others","Others","Others","tRNA modif","Others","Others","tRNA modif","DNA repair","Others","DNA repair","DNA repair","DNA repair")
leuconostoc$pvalue<--log10(leuconostoc$pvalue)
leuconostoc<-merge(leuconostoc,GO_sign,by=c("GO_description"))
leuconostoc<-leuconostoc[,-c(6:7,9:12,14:22)]
colnames(leuconostoc)[6] <- "Frequency"
colnames(leuconostoc)[7] <- "SD"

sodgloss <- top20_per_species[[4]]
sodgloss <- sodgloss %>% na.omit()
sodgloss$functions <- c("MF","BP","MF","BP","BP","MF","BP","BP","MF","MF","BP","MF")
sodgloss$Topfunctions<- c("Transposition","Transposition","Others","Transposition","Others","Others","Transposition","Transposition","Others","Others","Others","Others")
sodgloss$pvalue<--log10(sodgloss$pvalue) 
sodgloss<-merge(sodgloss,GO_sign,by=c("GO_description"))
sodgloss<-sodgloss[,-c(6:8,10:13,15:22)]
colnames(sodgloss)[6] <- "Frequency"
colnames(sodgloss)[7] <- "SD"

sodpra <- top20_per_species[[5]]
sodpra$functions <- c("MF","MF","CC","MF","MF","MF","MF","MF","BP","MF","MF","MF","MF","MF","MF","MF","CC","BP","MF","BP")
sodpra$Topfunctions<- c("Amino-acid process","Amino-acid process","Others","Others","Others","Others","Amino-acid process","Amino-acid process","Amino-acid process","Amino-acid process","Others","Others","Others","Others","Others","Others","Others","Others","Others","Others")
sodpra$pvalue<--log10(sodpra$pvalue) 
sodpra<-merge(sodpra,GO_sign,by=c("GO_description"))
sodpra<-sodpra[,-c(6:9,11:14,16:22)]
colnames(sodpra)[6] <- "Frequency"
colnames(sodpra)[7] <- "SD"

# Plot function
go_plot <- function(df, plotfilename) {
  df2 <- df %>% 
    # 1. Remove any grouping
    ungroup() %>% 
    # 2. Arrange by
    #   i.  facet group (period)
    #   ii. value (val)
    arrange(functions, pvalue) %>%
    # 3. Add order column of row numbers
    mutate(order = row_number())
  
  
  plot_df <- ggplot(df2, aes(order, pvalue)) +
    geom_col(aes(fill = functions), position = "dodge", width = 0.5) +
    scale_x_continuous(
      breaks = df2$order,
      labels = df2$GO_description) +
    # scale_y_continuous(expand = c(0, 0)) +
    facet_grid(functions ~ ., scales = "free", space = "free") +
    coord_flip() +
    theme_minimal() +
    scale_fill_brewer(palette = "Set2", name = "functions") +
    theme(panel.grid.major.y = element_blank()) + 
    labs(x="", y="p-value", title = df2$Species[1]) +
    theme(text = element_text(size=12),legend.position = "none",
          axis.ticks.y = element_blank()) + scale_y_log10("pvalue",
                                                          breaks = trans_breaks("log10", function(x) 10^x),
                                                          labels = trans_format("log10", math_format(10^.x)))
  ggsave(filename = plotfilename, plot = plot_df, width = 10, height = 15)
}

go_plot(sodgloss, plotfilename = "sodalis_gloss.pdf")
go_plot(klebsiella, plotfilename = "klesbiella.pdf")
go_plot(sodpra, plotfilename = "sodalis_prae.pdf")
go_plot(enterobacter, plotfilename = "enterobacter.pdf")
go_plot(leuconostoc, plotfilename = "leuconostoc.pdf")

#plot all top functions (SOL)
automn_final<-rbind(enterobacter,klebsiella,leuconostoc,sodgloss,sodpra)
automn_final<-automn_final[- grep("Others", automn_final$Topfunctions),]
#automn_final$Frequency<-log10(automn_final$Frequency)

automn_final<- automn_final %>% 
  ungroup() %>% 
  arrange(Species,Topfunctions, pvalue) %>%
  mutate(order = row_number())


automn_final$Species <- str_wrap(automn_final$Species, width = 15, exdent = 2)


ggplot(automn_final, aes(x = pvalue, y = reorder(GO_description,-pvalue))) + 
  geom_point(aes( size = Frequency, fill=Species), shape = 21, colour="black") + 
  theme_minimal()+
  facet_grid(vars(Species), scales="free", space="free_x") +
  xlab("-log10(p.value)") +
  ylab("GO description") + scale_fill_brewer(palette = "Oranges") + scale_size_continuous(limits = c(0,200),range = c(0,10), breaks = c(1,50,100,200))

