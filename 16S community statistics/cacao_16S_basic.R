#Set working directory
setwd("~/Library/CloudStorage/GoogleDrive-nn33@hawaii.edu/.shortcut-targets-by-id/1-4rD_zMa5xHszdH4D_WkYU9f6pwkT0Kq/cacao_combined/cacao_physeq_full_analysis/cacao_16S_ps")

#load libraries etc.
library(tidyverse)
library(phyloseq)
library(vegan)
library(gridExtra)
library(multcompView)

#Getting ready
theme_set(theme_bw())
set.seed(1943)
cacao_colors = c("Urban Garden Center"="#FD8D3C", "Kualoa Ranch" = "#46AEA0")

#import unrarefied contingency table into PhyloSeq Object
cacao_16S = import_biom("cacao_prok.biom", treefilename = "cacao_rooted-tree-filtered2-prok.nwk", refseqfilename = "cacao_prok-only.fasta")

#rename columns
colnames(tax_table(cacao_16S)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#load metadata
cacao_16S_metadata <- read.csv("cacao_metadata_combined.csv")   

#Add metadata to main dataset
cacao_16S_metadata$group <- paste(cacao_16S_metadata$Location,cacao_16S_metadata$selection_days_ferm)
# str(cacao_16S_metadata)
# View(cacao_16S_metadata)

#Create taxa table
tax_table_cacao_16S_0 <- as.data.frame(tax_table(cacao_16S))
# View(tax_table_cacao_16S_0)
# write.csv(tax_table_cacao_16S_0, "tax_table_cacao_16S_0.csv")

tax_table_cacao_16S <- tax_table_cacao_16S_0

#Rename taxa
tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Acidobacteriota", "Acidobacteria", x)}))
tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Actinobacteriota", "Actinobacteria", x)}))
tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Armatimonadota", "Armatimonadetes", x)}))
tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Bacteroidota", "Bacteroidetes", x)}))
tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Gemmatimonadota", "Gemmatimonadetes", x)}))
tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Halobacterota", "Halobacteria", x)}))
tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Planctomycetota", "Planctomycetes", x)}))
tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Verrucomicrobiota", "Verrucomicrobia", x)}))

#Remove taxonomic designations
tax_table_cacao_16S$Kingdom<- gsub("d__", "", tax_table_cacao_16S$Kingdom)
tax_table_cacao_16S$Phylum <- gsub("p__", "", tax_table_cacao_16S$Phylum)
tax_table_cacao_16S$Class <- gsub("c__", "", tax_table_cacao_16S$Class)
tax_table_cacao_16S$Order <- gsub("o__", "", tax_table_cacao_16S$Order)
tax_table_cacao_16S$Family <- gsub("f__", "", tax_table_cacao_16S$Family)
tax_table_cacao_16S$Genus <- gsub("g__", "", tax_table_cacao_16S$Genus)
tax_table_cacao_16S$Species <- gsub("s__", "", tax_table_cacao_16S$Species)
tax_table_cacao_16S$Fam_Gen <- paste(tax_table_cacao_16S$Family," | ", tax_table_cacao_16S$Genus,sep="")

row.names(tax_table_cacao_16S) <- row.names(tax_table_cacao_16S_0)
tax_table(cacao_16S) <- as.matrix(tax_table_cacao_16S)

row.names(cacao_16S_metadata) <- cacao_16S_metadata$X.SampleID

cacao_16S_metadata$SampleID <- cacao_16S_metadata$X.SampleID

sample_data(cacao_16S) <- cacao_16S_metadata

# View(data.frame(sample_data(cacao_16S)))

tax_table_cacao_16S <- as.data.frame(tax_table(cacao_16S))
# write.csv(tax_table_cacao_16S, "tax_table_cacao_16S.csv")

# not necessary: #subset_taxa to include only those identified as bacteria
# not necessary: cacao_16S_sub0 = subset_taxa(cacao_16S, Kingdom=="Bacteria")
cacao_16S = subset_samples(cacao_16S, Amplicon== "16S")
cacao_16S = subset_samples(cacao_16S, Description != "16SMOCK")
cacao_16S = subset_samples(cacao_16S, SampleID != "02.2.5.1.16S.a")
cacao_16S = subset_samples(cacao_16S, days_ferment != "3" & days_ferment !="6" & days_ferment != "1")
# View(data.frame(sample_data(cacao_16S)))

#rarecurve(t(otu_table(cacao_16S)), step=50, cex=0.5)

total_depth_0 <- colSums((data.frame(otu_table(cacao_16S))))
# View(total_depth_0)

cacao_16S_sub = prune_samples(sample_sums(cacao_16S)>= 1231, cacao_16S)

total_depth_0_1 <- colSums((data.frame(otu_table(cacao_16S))))
# View(total_depth_0_1)

# richnessp <- plot_richness(cacao_16S, x="selection_days_ferm", measures=c("Shannon", "Simpson", "Observed"), color="Location")
# richnessp

#----Rarefy samples----
set.seed(1943)
cacao_16S_rare = rarefy_even_depth((cacao_16S),rngseed = TRUE)
#rarecurve(t(otu_table(cacao_16S_rare)), step=50, cex=0.5)

#sum to make sure rarefy was successful
total_depth_cacao_16S_rare <- colSums((data.frame(otu_table(cacao_16S_rare))))
# View(total_depth_cacao_16S_rare)

#calculate relative abundance
cacao_16S_relabund = transform_sample_counts(cacao_16S_rare, function(x) x / sum(x) )
write.csv(otu_table(cacao_16S_relabund), "cacao_16S_relabund.csv")

#Build rarefaction curves
# rel_cacao_16S_0 <- prune_samples(sample_sums(cacao_16S_sub2) >= 50, cacao_16S_sub)
# rarecurve(t(otu_table(rel_cacao_16S_0)), step=50, cex=0.5)
# View(data.frame(otu_table(rel_cacao_16S_0)))

#Bray-Curtis distance matrix
#cacao_16S_rareBray <- phyloseq::distance(cacao_16S_rare, method = "bray")
#cacao_16S_rare_data = data.frame(sample_data(cacao_16S_rare))

#Weighted Unifrac distance matrix
cacao_16S_rareUnifrac <- phyloseq::distance(cacao_16S_rare, method = "wunifrac")

#change metadata table into data frame
cacao_16S_rare_data = data.frame(sample_data(cacao_16S_rare))
# View(cacao_16S_rare_data)

########################################################
#Richness
Observed_otus <- estimate_richness(cacao_16S_rare, measures="Observed")

#check to see if rownames available
has_rownames(Observed_otus)#rownames available

#remove the 'X' in front of each rownames
rownames(Observed_otus) <- substr(rownames(Observed_otus),2,nchar(rownames(Observed_otus)))

#convert rownames to actual column
Observed_otus <- rownames_to_column(Observed_otus, var="SampleID")

#join the observed OTU data to main metadata table
cacao_16S_rare_data <- left_join(Observed_otus, cacao_16S_rare_data, by="SampleID", row.names="SampleID")

##===================================================================================================
##                                  Build alpha diversity Plots (Supplemental Figure)
##===================================================================================================
observed_asv_location_plot <- ggplot(cacao_16S_rare_data, aes(x=Location, y=Observed)) +
  #stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(outlier.color="gray", show.legend = FALSE) +
  theme_classic() +
  scale_fill_manual(values=c("#8DD3C7", "#FDB462", "#B3DE68", "#FB8072"), name="Sample type") +
  scale_x_discrete(limit=c("Kualoa Ranch", "Urban Garden Center"), 
                   labels = c("Kualoa Ranch \n phyllosphere", "Urban Garden Center \n phyllosphere")) +
  labs(x="", y = "Observed bacterial ASVs") +
  ggtitle("C) Site") +
  theme(plot.tag = element_text(size=18))

observed_asv_selection_plot <- ggplot(cacao_16S_rare_data, aes(x=tree_selection, y=Observed)) +
  #stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(outlier.color="gray", show.legend = FALSE) +
  theme_classic() +
  scale_fill_manual(values=c("#8DD3C7", "#FDB462", "#B3DE68", "#FB8072"), name="Sample type") +
  scale_x_discrete(limit=c("3", "5"), 
                   labels = c("HSCT3 phyllosphere", "HSCT5 phyllosphere")) +
  labs(x="", y = "Observed bacterial ASVs") +
  ggtitle("D) Cacao variety") +
  theme(plot.tag = element_text(size=18))


#Plot together
observed_asv_location_plot + observed_asv_selection_plot

#Save the plot
ggsave("FigureSX-alphadiv.pdf", device="pdf", width=10.5, height=5)

##===================================================================================================
##                                  Testing ANOVA assumptions; perform ANOVA
##===================================================================================================
# Using the Shapiro Test -- if p-value not significant = data is normal = OK for ANOVA
shapiro.test(cacao_16S_rare_data$Observed)#NS
#shapiro.test(metadata$faith_pd)
#shapiro.test(metadata$shannon_entropy)

# Test for homogeneity of variances -- if p-value not significant = data is homogeneous = OK for ANOVA
bartlett.test(cacao_16S_rare_data$Observed~cacao_16S_rare_data$Location)#NS
bartlett.test(cacao_16S_rare_data$Observed~cacao_16S_rare_data$tree_selection)#S

#bartlett.test(metadata$faith_pd~metadata$HostSampleType)#NS
#bartlett.test(metadata$shannon_entropy~metadata$HostSampleType)#S

anova_observed_asv <- aov(cacao_16S_rare_data$Observed ~ cacao_16S_rare_data$Location)
summary(anova_observed_asv)#S
tukey_observed_asv <- TukeyHSD(anova_observed_asv, ordered="TRUE")
tukey_observed_asv

anova_observed_asv <- aov(cacao_16S_rare_data$Observed ~ cacao_16S_rare_data$tree_selection)
summary(anova_observed_asv)#NS
tukey_observed_asv <- TukeyHSD(anova_observed_asv, ordered="TRUE")
tukey_observed_asv

##===================================================================================================
##                                  Community Composition Comparisons
##===================================================================================================
#Location differences
set.seed(1943)
cacao_bac_adonis <- adonis2(cacao_16S_rareUnifrac~Location, data=cacao_16S_rare_data, strata=cacao_16S_rare_data$tree_selection, permutations=999, method ="unifrac")
cacao_bac_adonis

#Variety differences
set.seed(1943)
cacao_bac_adonis <- adonis2(cacao_16S_rareUnifrac~tree_selection, data=cacao_16S_rare_data, strata=cacao_16S_rare_data$Location, permutations=999, method ="unifrac")
cacao_bac_adonis

#Testing interactions
set.seed(1943)
cacao_bac_adonis <- adonis2(cacao_16S_rareUnifrac~Location*days_ferment*tree_selection, data = cacao_16S_rare_data, permutations=999, method ="bray")
cacao_bac_adonis
#plot(cacao_fun_adonis$aov.tab)

###################################
#FIGURE 1A
###################################
bac_uni_Wp +
  theme_bw() +
  geom_point(size=3, show.legend=FALSE) +
  ggtitle("A) Bacteria") +
  theme(text=element_text(size=12)) +
  scale_color_manual(values=cacao_colors) +
  scale_shape_manual(values=c(16, 15))

#Save the plot
ggsave("Figure1A-ordination.pdf", device="pdf", width=5, height=4)

# svg("bac_uni_Wp.svg", width = 8, height = 4, pointsize=12)
# plot(bac_uni_Wp)
# dev.off()

###############################################
#Some exploratory analyses
# plot_heatmap(physeq = rel_cacao_16S, method = "MDS", distance = "wUniFrac", 
#              title = "weighted-UniFrac", taxa.label = FALSE, sample.label="Location")

#Relative abundance is kinda working....
# top40 <- names(sort(taxa_sums(rel_cacao_16S), decreasing=TRUE))[1:40]
# rel_cacao_16S.top40 <- prune_taxa(top40, rel_cacao_16S)
# rel_cacao_16S.top40_2 <- merge_samples(rel_cacao_16S.top40, "group")
# rel_cacao_16S.top40_3 <- transform_sample_counts(rel_cacao_16S.top40_2, function(x) x / sum(x))
# 
# View(data.frame(sample_data(rel_cacao_16S.top40_3)))
# 
# plot_bar(rel_cacao_16S.top40_3, x="tree_selection", fill="Family") + facet_wrap(~Location, scales="free_x")

# 
# cacao_bac_dbrda <- dbrda(rel_cacao_16S.dis ~ tree_selection + Location, data.env, scannf = TRUE, distance = "bray")
# cacao_bac_dbrdaplot <- plot(cacao_bac_dbrda, type = "p", col = data.env$Col)
# cacao_bac_dbrdaplot
# anova(cacao_bac_dbrda, by="terms", permutations = 999)
# cacao_bac_dbrda
# 
# cacao_bac_dbrdadf <- as.data.frame(cacao_bac_dbrdaplot$sites)
# cacao_bac_dbrdavec <- as.data.frame(cacao_bac_dbrdaplot$biplot)
# cacao_bac_dbrdavec$env <- rownames(cacao_bac_dbrdavec) 
# 
# 
# 
# cacao_bac_dbrda_p1 <- ggplot(data=cacao_bac_dbrdadf, aes(dbRDA1, dbRDA2)) + geom_point(aes(colour = data.env$Location, shape = data.env$selection_days_ferm), stroke=1, size=3)
# cacao_bac_dbrda_p1


# svg("cacao_bac_dbrda.svg", width = 8, height = 4, pointsize=12)
# plot(cacao_bac_dbrda_p1)
# dev.off()



# cacao_bac_dbrda_p2 <- cacao_bac_dbrda_p1 + geom_segment(data=cacao_bac_dbrdavec,aes(x=0,xend=dbRDA1,y=0,yend=dbRDA2), colour="grey") + geom_text(data=cacao_bac_dbrdavec,aes(x=dbRDA1,y=dbRDA2,label=env,size=4))
# cacao_bac_dbrda_p2




#constrained ordinations-we'll probably go with dbRDA for now
# rel_cacao_cca <- cca(rel_cacao_16S_vegan ~tree_selection + Location, data.env, scannf = TRUE)
# plot(rel_cacao_cca)
# 
# data.env$tree_selection <- as.character(data.env$tree_selection)
# 
# cacao_bac_dbrda <- dbrda(rel_cacao_16S.dis ~tree_selection + Location, data.env, scannf = TRUE, distance = "bray")
# plot(cacao_bac_dbrda)
# anova(cacao_bac_dbrda, by="terms")
# 
# cacao_bac_dbrdadf <- as.data.frame(cacao_bac_dbrda$sites)
# 
# cacao_bac_dbrdavec <- as.data.frame(cacao_bac_dbrda$biplot) 
# 
# cacao_bac_dbrdavec$env <- rownames(cacao_bac_dbrda) 
# 
# 
# ggcacao_bac1 <- ggplot(data=cacao_bac_dbrdadf, aes(dbRDA1, dbRDA2)) + geom_point(aes(colour = data.env$Location), stroke=1, size=3)
# ggcacao_bac1
# 
# View(cacao_bac_dbrdadf)
# 
# ggcacao_bac2 <- ggcacao_bac1 + geom_segment(data=cacao_bac_dbrdavec,aes(x=0,xend=dbRDA1,y=0,yend=dbRDA2), colour="grey") + geom_text(data=cacao_bac_dbrdavec,aes(x=dbRDA1,y=dbRDA2,label=env,size=4))
# ggcacao_bac2
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # conord3 <- capscale(rel_cacao_16S_vegan ~Chlor, data.env, scannf = TRUE, distance = "bray")
# # plot(conord3)
# 
# 
# 
# conord2a <- dbrda(con.dis ~ Chlor + Magnesium, data.env2, scannf = TRUE, distance = "bray")
# conord2aplot <- plot(conord2a, type = "p", col = data.env2$Col)
# conord2aplot
# anova(conord2a, by="terms", permutations = 999)
# conord2a
# 
# conord2adf <- as.data.frame(conord2aplot$sites)
# conord2avec <- as.data.frame(conord2aplot$biplot)
# conord2avec$env <- rownames(conord2avec) 
# 
# 
# 
# ggconord1 <- ggplot(data=conord2adf, aes(dbRDA1, dbRDA2)) + geom_point(aes(colour = data.env2$Chlor), stroke=1, size=3)
# 
# ggconord2 <- ggconord1 + geom_segment(data=conord2avec,aes(x=0,xend=dbRDA1,y=0,yend=dbRDA2), colour="grey") + geom_text(data=conord2avec,aes(x=dbRDA1,y=dbRDA2,label=env,size=4))
# ggconord2
# 
# fig <- ordiplot(conord2a, type = "none")
# points(fig, "sites", pch=21, col="red", bg="yellow")
# text(fig, "species", col="blue", cex=0.9)
# fig
# 
# fig2 <- ordisurf(conord2a~Chlor, data.env2, method = "GCV.Cp", isotropic=FALSE, bs="tp", knots =c(3,4), fx=FALSE,select=FALSE)
# fig2
# 
# summary(fig2)
# conord2adf$chlorTF <-datachlor2$ChlorTF 
# 
# qplot1 <- qplot(data=conord2adf, x=dbRDA1, y=dbRDA2, colour=chlorTF) + stat_ellipse()
# 
# qplot1
# 
# adonischlortf <- adonis(con.dis~chlorTF, data=conord2adf, permutations=999)
# adonischlortf
# plot(adonischlortf$aov.tab)
# 
# qplot2 <- qplot1  + geom_point(aes(colour = data.env2$Chlor), stroke=1, size=3)
# qplot2
# datachlor2 <- as.data.frame(data.env2$Chlor==2)
# View(datachlor2)
# row.names(datachlor2) <- data.env2$SampleID
# datachlor2$Chlor <- data.env2$Chlor
# datachlor2$ChlorTF <-datachlor2$`data.env2$Chlor == 2`
# View(data.env2)
# View(datachlor2)
# datachlor2$chlor2code <- c(0,	0,	1,	1,	0,	1,	1,	0,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	1,	0,	1,	1,	0,	0,	0,	1,	0)
# 
# str(datachlor2)
# 
# 
# 
# 
# conord3 <- capscale(con.dis ~Chlor, data.env, scannf = TRUE, distance = "bray")
# plot(conord3)
# 
