#set working directory
setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1--pvx1vLSBijmxHMCZPrxKlyFBLp2jqF/Rick/cacao_combined/cacao_physeq_full_analysis/cacao_Fungal_ITS_ps")

#load libraries
library(tidyverse)
library(phyloseq)
library(vegan)
library(gridExtra)

#Getting ready
theme_set(theme_bw())
set.seed(1943)
cacao_colors = c("Urban Garden Center"="#FD8D3C", "Kualoa Ranch" = "#46AEA0")

#import unrarefied contingency table
cacao_ITS = import_biom( "cacao_fungi.biom", refseqfilename = "cacao_fungi-only.fasta")

#rename columns 
colnames(tax_table(cacao_ITS)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

otu_ITS_df <- as.data.frame(otu_table(cacao_ITS))
# write.csv(otu_ITS_df, "otu_ITS_df.csv")

#import metadata table
cacao_ITS_metadata <- read.csv("cacao_metadata_combined.csv")

#add new column
cacao_ITS_metadata$location_days_ferm <- paste(cacao_ITS_metadata$Location,cacao_ITS_metadata$days_ferment)

# str(cacao_ITS_metadata)
# View(cacao_ITS_metadata)

row.names(cacao_ITS_metadata) <- cacao_ITS_metadata$X.SampleID

cacao_ITS_metadata$SampleID <- cacao_ITS_metadata$X.SampleID


sample_data(cacao_ITS) <- cacao_ITS_metadata

# View(data.frame(sample_data(cacao_ITS)))

tax_table_cacao_ITS_0 <- as.data.frame(tax_table(cacao_ITS))
tax_table_cacao_ITS <- tax_table_cacao_ITS_0

tax_table_cacao_ITS$Kingdom<- gsub("k__", "", tax_table_cacao_ITS$Kingdom)
tax_table_cacao_ITS$Phylum <- gsub("p__", "", tax_table_cacao_ITS$Phylum)
tax_table_cacao_ITS$Class <- gsub("c__", "", tax_table_cacao_ITS$Class)
tax_table_cacao_ITS$Order <- gsub("o__", "", tax_table_cacao_ITS$Order)
tax_table_cacao_ITS$Family <- gsub("f__", "", tax_table_cacao_ITS$Family)
tax_table_cacao_ITS$Genus <- gsub("g__", "", tax_table_cacao_ITS$Genus)
tax_table_cacao_ITS$Species <- gsub("s__", "", tax_table_cacao_ITS$Species)
tax_table_cacao_ITS$Fam_Gen <- paste(tax_table_cacao_ITS$Family,"_",tax_table_cacao_ITS$Genus,sep="")

row.names(tax_table_cacao_ITS) <- row.names(tax_table_cacao_ITS_0)
tax_table(cacao_ITS) <- as.matrix(tax_table_cacao_ITS)

row.names(cacao_ITS_metadata) <- cacao_ITS_metadata$X.SampleID

cacao_ITS_metadata$SampleID <- cacao_ITS_metadata$X.SampleID


sample_data(cacao_ITS) <- cacao_ITS_metadata
#View(cacao_ITS_metadata)

tax_table_cacao_ITS <- as.data.frame(tax_table(cacao_ITS))
# View(tax_table_cacao_ITS)

# write.csv(tax_table_cacao_ITS, "tax_table_cacao_ITS.csv")


###remove OTUs from negative controls. after looking at the original otu table, ITSNEG.ITS.a had most of the bad taxa, and had the same single bad taxa found in ITSNEG.ITS.b, so it was only necessary to remove taxa from ITSNEG.ITS.a as I did below
cacao_ITS_neg = subset_samples(cacao_ITS, SampleID == "ITSNEG.ITS.a")
cacao_ITS_neg = filter_taxa(cacao_ITS_neg, function(x) sum(x) > 0,TRUE)
cacao_ITS_neg_tax <- as.data.frame(tax_table(cacao_ITS_neg))
cacao_ITS_neg_otu <- as.data.frame(otu_table(cacao_ITS_neg))
#View(cacao_ITS_neg_otu)
badTaxa <- row.names(cacao_ITS_neg_otu)
#View(badTaxa)
allTaxa = taxa_names(cacao_ITS)
myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
cacao_ITS_trim_0 = prune_taxa(myTaxa, cacao_ITS)
total_depth_trim_0 <- colSums((data.frame(otu_table(cacao_ITS_trim_0))))


cacao_ITS = subset_samples(cacao_ITS_trim_0, Amplicon== "ITS")
cacao_ITS = subset_samples(cacao_ITS, Description != "ITSMOCK")
# cacao_ITS = subset_samples(cacao_ITS, SampleID != "02.2.5.1.ITS.a")
cacao_ITS = subset_samples(cacao_ITS, days_ferment != "3" & days_ferment !="6")
#View(data.frame(sample_data(cacao_ITS)))

otu_ITS_df2 <- as.data.frame(otu_table(cacao_ITS))
#View(otu_ITS_df2)
# write.csv(otu_ITS_df2, "otu_ITS_df2.csv")

###Above, I'm trying to ensure that the correct metadata are associated with each sample, as I had to import the metadata separately, and I'm not sure how "smart" the phyloseq software is. When filtering by imported metadatata, it seems to be filtering out the correct samples from the otu table, as assessed by the otu_table sample id before and after filtering, and the abundance of otus within each sample.


rarecurve(t(otu_table(cacao_ITS)), step=50, cex=0.5)

cacao_ITS_depth <- colSums((data.frame(otu_table(cacao_ITS))))
#View(cacao_ITS_depth)

cacao_ITS <- prune_samples(sample_sums(cacao_ITS) >= 50, cacao_ITS)

rarecurve(t(otu_table(cacao_ITS)), step=50, cex=0.5)
cacao_ITS_depth <- colSums((data.frame(otu_table(cacao_ITS))))
#View(cacao_ITS_depth)

cacao_ITS_rare = rarefy_even_depth(cacao_ITS)
rarecurve(t(otu_table(cacao_ITS_rare)), step=50, cex=0.5)
total_depth_0_1 <- colSums((data.frame(otu_table(cacao_ITS_rare))))
#View(total_depth_0_1)

total_depth_cacao_ITS_rare <- colSums((data.frame(otu_table(cacao_ITS_rare))))
#View(total_depth_cacao_ITS_rare)

cacao_ITS_rareBray <- phyloseq::distance(cacao_ITS_rare, method = "bray")

cacao_ITS_rare_data = data.frame(sample_data(cacao_ITS_rare))
#View(cacao_ITS_rare_data)

########################################################
#Community composition
#Checking general model
set.seed(1943)
cacao_fun_adonis <- adonis(cacao_ITS_rareBray ~Location*days_ferment*tree_selection, data = cacao_ITS_rare_data, permutations=999, method ="bray")
cacao_fun_adonis
#plot(cacao_fun_adonis$aov.tab)

#blocked by phyllosphere/fermentation (days); get site
set.seed(1943)
cacao_fun_adonis_location <- adonis(cacao_ITS_rareBray ~Location, strata = cacao_ITS_rare_data$days_ferment, data = cacao_ITS_rare_data, permutations=999, method ="bray")
cacao_fun_adonis_location
#plot(cacao_fun_adonis_location$aov.tab)

#blocked by phyllosphere/fermentation (days); get selection
cacao_fun_adonis_location <- adonis(cacao_ITS_rareBray ~tree_selection, strata = cacao_ITS_rare_data$days_ferment, data = cacao_ITS_rare_data, permutations=999, method ="bray")
cacao_fun_adonis_selection
#plot(cacao_fun_adonis_location$aov.tab)

#blocked by site (cacao selection)
set.seed(1943)
cacao_fun_adonis_strata_select <- adonis(cacao_ITS_rareBray ~tree_selection+days_ferment*tree_selection, data = cacao_ITS_rare_data, strata = cacao_ITS_rare_data$Location, permutations=999, method ="bray")
cacao_fun_adonis_strata_select
#plot(cacao_fun_adonis_strata_select$aov.tab)

#blocked by site (phyllosphere vs fermentation)
set.seed(1943)
cacao_fun_adonis_strata_select <- adonis(cacao_ITS_rareBray ~days_ferment*tree_selection, data = cacao_ITS_rare_data, strata = cacao_ITS_rare_data$Location, permutations=999, method ="bray")
cacao_fun_adonis_strata_select
#plot(cacao_fun_adonis_strata_select$aov.tab)

#blocked by site (phyllosphere vs fermentation)
#This is the MS data comparing phyllosphere vs. fermentation
set.seed(1943)
cacao_fun_adonis_strata_ferm <- adonis(cacao_ITS_rareBray ~days_ferment, data = cacao_ITS_rare_data, strata = cacao_ITS_rare_data$Location, permutations=999, method ="bray")
cacao_fun_adonis_strata_ferm
#plot(cacao_fun_adonis_strata_ferm$aov.tab)

set.seed(1943)
cacao_fun_adonis_strata_ferm2 <- adonis(cacao_ITS_rareBray ~tree_selection, data = cacao_ITS_rare_data, strata = cacao_ITS_rare_data$days_ferment, permutations=999, method ="bray")
cacao_fun_adonis_strata_ferm2

###############################################
#Pair-wise adonis
#This is the data reported in the manuscript
set.seed(1943)
pairwise.adonis2(cacao_ITS_rareBray ~selection_days_ferm, data = cacao_ITS_rare_data, permutations=999, method ="bray")
pairwise.adonis2(cacao_ITS_rareBray ~location_days_ferm, data = cacao_ITS_rare_data, permutations=999, method ="bray")
###############################################

# richnessp_ITS <- plot_richness(cacao_ITS_rare, x="selection_days_ferm", shape = "group", color="Site", measures=c("Shannon"), color="Location")
# richnessp_ITS

cacao_ITS_sub_alpha_rare <- plot_richness(cacao_ITS_rare, x="selection_days_ferm", color="Site", measures = c("Shannon", "Observed")) 
cacao_ITS_sub_alpha_rare

cacao_ITS_sub_alphadt_rare <- data.table(cacao_ITS_sub_alpha_rare$data)
View(cacao_ITS_sub_alphadt_rare)
#write.csv(cacao_ITS_sub_alphadt_rare, "cacao_ITS_sub_alphadt_rare.csv")


rel_cacao_ITS_0 <- transform_sample_counts(cacao_ITS_rare, function(x) x/(sum(x)))
rel_cacao_ITS = filter_taxa(rel_cacao_ITS_0, function(x) sum(x) > 0, TRUE)

fung_cacao_relgen_all <- tax_glom(rel_cacao_ITS, taxrank="Genus")
fung_cacao_relgen_all2 <- transform_sample_counts(fung_cacao_relgen_all, function(x) x/(sum(x)))

fung_cacao_relgen_all3 <- filter_taxa(fung_cacao_relgen_all2, function(x) mean(x) > 1e-3, TRUE)
fung_cacao_relgen_all4 <- transform_sample_counts(fung_cacao_relgen_all3, function(x) x/(sum(x)))
plot_bar(fung_cacao_relgen_all4, fill = "Genus")

#write.csv(otu_table(fung_cacao_relgen_all4), "fung_cacao_relgen_all4.csv")
#write.csv(tax_table(fung_cacao_relgen_all4), "fung_cacao_taxa_relgen_all.csv")



# fung_cacao_relfam_gen_all <- tax_glom(rel_cacao_ITS, taxrank="Fam_Gen")
# fung_cacao_relfam_gen_all2 <- transform_sample_counts(fung_cacao_relfam_gen_all, function(x) x/(sum(x)))
# 
# fung_cacao_relfam_gen_all3 <- filter_taxa(fung_cacao_relfam_gen_all2, function(x) mean(x) > 1e-3, TRUE)
# fung_cacao_relfam_gen_all4 <- transform_sample_counts(fung_cacao_relfam_gen_all3, function(x) x/(sum(x)))
# plot_bar(fung_cacao_relfam_gen_all4, fill = "Fam_Gen")
# 
# write.csv(otu_table(fung_cacao_relfam_gen_all4), "fung_cacao_relfam_gen_all4.csv")
# write.csv(tax_table(fung_cacao_relfam_gen_all4), "fung_cacao_taxa_rel_fam_gen_all.csv")




#NMDS plot
rel_cacao_ITS_vegan <- t(otu_table(rel_cacao_ITS))
rel_cacao_ITS.dis <- vegdist(rel_cacao_ITS_vegan, method="bray")
rel_cacao_ITS.mdso <- isoMDS(rel_cacao_ITS.dis)
stressplot(rel_cacao_ITS.mdso, rel_cacao_ITS.dis)
rel_cacao_ITS.mds <- metaMDS(rel_cacao_ITS_vegan, trace = FALSE, distance = "bray")
rel_cacao_ITS.mds
rel_cacao_ITS.mds2 <- data.frame(MDS1=rel_cacao_ITS.mds$points[,1], MDS2=rel_cacao_ITS.mds$points[,2])
stressplot(rel_cacao_ITS.mds, rel_cacao_ITS.dis)
plot(rel_cacao_ITS.mds2)

#fit environmental data 
data.env <- data.frame(sample_data(rel_cacao_ITS))

# ef1 <- envfit(rel_eri.mds$points~, data.env, permu=999)
# ef1.df<- as.data.frame(ef1$vectors$arrows*sqrt(ef1$vectors$r))
# ef1.df$env <- rownames(ef1.df)
# ef1c.df <- as.data.frame(ef1$factors$centroids)
# ef1c.df$env2 <- rownames(ef1c.df)
# ef1

NMDSplot2 = ggplot(data = rel_cacao_ITS.mds2, aes(MDS1, MDS2)) + geom_point(aes(colour = data.env$Location, shape = data.env$selection_days_ferm), stroke=2, size=3) 
NMDSplot2


#PCoA Plot
fun_brayw_dmx = ordinate(rel_cacao_ITS, method="PCoA", distance="bray", weighted = TRUE)
fun_brayw = plot_ordination(rel_cacao_ITS, fun_brayw_dmx, type="samples", color="Location", shape="selection_days_ferm")
fun_brayw

# fun_uni_Wp = plot_ordination(rel_cacao_ITS, fun_uni_W, type = "samples", color = "Location", shape= "tree_selection")
# fun_uni_Wp
# fun_uni_UWp = plot_ordination(rel_cacao_ITS, fun_uni_UW, type = "samples", color = "Location", shape= "tree_selection")
# fun_uni_UWp

fun_brayw +
  theme_bw() +
  geom_point(size=3, show.legend=FALSE) +
  ggtitle("B) Fungi") +
  theme(text=element_text(size=12)) +
  scale_color_manual(values=cacao_colors) +
  scale_shape_manual(values=c(16, 17, 15, 8))

#Save the plot
ggsave("Figure1B-ordination.pdf", device="pdf", width=5, height=4)

svg("Figure1A.svg", width = , height = 4)
plot(fun_brayw)
dev.off()


#svg("fun_brayw.svg", width = , height = 4, pointsize=12)
#plot(fun_brayw)
#dev.off()





#Some exploratory analyses

cacao_ITS_dbrda <- dbrda(rel_cacao_ITS.dis ~ tree_selection + Location, data.env, scannf = TRUE, distance = "bray")
cacao_ITS_dbrdaplot <- plot(cacao_ITS_dbrda, type = "p", col = data.env$Col)
cacao_ITS_dbrdaplot
anova(cacao_ITS_dbrda, by="terms", permutations = 999)
cacao_ITS_dbrda

cacao_ITS_dbrdadf <- as.data.frame(cacao_ITS_dbrdaplot$sites)
cacao_ITS_dbrdavec <- as.data.frame(cacao_ITS_dbrdaplot$biplot)
cacao_ITS_dbrdavec$env <- rownames(cacao_ITS_dbrdavec) 


cacao_ITS_dbrda_p1 <- ggplot(data=cacao_ITS_dbrdadf, aes(dbRDA1, dbRDA2)) + geom_point(aes(colour = data.env$Location, shape = data.env$selection_days_ferm), stroke=1, size=3)
cacao_ITS_dbrda_p1

svg("cacao_ITS_dbrda_p1.svg", width = 8, height = 4, pointsize=12)
plot(cacao_ITS_dbrda_p1)
dev.off()

cacao_ITS_WF_rare = subset_samples(cacao_ITS_rare, Location == "Kualoa Ranch")

# cacao_ITS_WF_rare = rarefy_even_depth(cacao_ITS_WF)
rarecurve(t(otu_table(cacao_ITS_WF_rare)), step=50, cex=0.5)
total_depth_0_1_WF <- colSums((data.frame(otu_table(cacao_ITS_WF_rare))))
# View(total_depth_0_1_WF)

total_depth_cacao_ITS_WF_rare <- colSums((data.frame(otu_table(cacao_ITS_WF_rare))))
# View(total_depth_cacao_ITS_WF_rare)

cacao_ITS_WF_rareBray <- phyloseq::distance(cacao_ITS_WF_rare, method = "bray")

cacao_ITS_WF_rare_data = data.frame(sample_data(cacao_ITS_WF_rare))

# View(cacao_ITS_WF_rare_data)

cacao_fun_WF_adonis <- adonis(cacao_ITS_WF_rareBray ~days_ferment + tree_selection, data = cacao_ITS_WF_rare_data, permutations=999, method ="bray")
cacao_fun_WF_adonis
plot(cacao_fun_WF_adonis$aov.tab)

cacao_ITS_UG_rare = subset_samples(cacao_ITS_rare, Location ==
                                               "Urban Garden Center")

# cacao_ITS_UG_rare = rarefy_even_depth(cacao_ITS_UG)

rarecurve(t(otu_table(cacao_ITS_UG_rare)), step=50, cex=0.5)
total_depth_0\_1_UG \<-
  colSums((data.frame(otu_table(cacao_ITS_UG_rare)))) \#
View(total_depth_0\_1_UG)

total_depth_cacao_ITS_UG_rare \<-
  colSums((data.frame(otu_table(cacao_ITS_UG_rare)))) \#
View(total_depth_cacao_ITS_UG_rare)

cacao_ITS_UG_rareBray \<- phyloseq::distance(cacao_ITS_UG_rare, method =
                                               "bray")

cacao_ITS_UG_rare_data = data.frame(sample_data(cacao_ITS_UG_rare))

# View(cacao_ITS_UG_rare_data)

cacao_fun_UG_adonis \<- adonis(cacao_ITS_UG_rareBray \~days_ferment +
                                 tree_selection, data = cacao_ITS_UG_rare_data, permutations=999, method
                               ="bray") cacao_fun_UG_adonis plot(cacao_fun_UG_adonis\$aov.tab)
