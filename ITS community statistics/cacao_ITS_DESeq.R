#set working directory
setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1--pvx1vLSBijmxHMCZPrxKlyFBLp2jqF/Rick/cacao_combined/cacao_physeq_full_analysis/cacao_Fungal_ITS_ps")

#load libraries
library(phyloseq)
library(DESeq2)
library(tidyverse)

#Getting ready
theme_set(theme_bw())
cacao_colors = c("Ascomycota"="#FD8D3C", "Basidiomycota" = "#46AEA0")

#Import a phyloseq object (biom should have taxonomy attached), start with UNRAREFIED otu table
cacao_ITS = import_biom( "cacao_fungi.biom", refseqfilename = "cacao_fungi-only.fasta")

#Rename columns 
colnames(tax_table(cacao_ITS)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Save otu table (optional) -- this will save a simple OTU table without taxonomy
#cacao_ITS_otu <- data.frame(otu_table(cacao_ITS))
#write.csv(cacao_ITS_otu,"cacao_ITS_otu_initial.csv")

#Read in metadata; standard QIIME2 metadata table can be used here
cacao_ITS_metadata <- read.csv("cacao_metadata_combined.csv")

#Create another row for "SampleID", create a dataframe for the metadata
row.names(cacao_ITS_metadata) <- cacao_ITS_metadata$X.SampleID
cacao_ITS_metadata$SampleID <- cacao_ITS_metadata$X.SampleID
sample_data(cacao_ITS) <- cacao_ITS_metadata
# #View(data.frame(sample_data(cacao_ITS)))

#De-QIIME-ify the taxa table -- this will separate taxonomic ranks into each separate columns. This _0 table is necessary downstream.
tax_table_cacao_ITS_0 <- as.data.frame(tax_table(cacao_ITS))
#write.csv(tax_table_cacao_ITS_0, "tax_table_cacao_ITS_0.csv")

#Make a copy of the table
tax_table_cacao_ITS <- tax_table_cacao_ITS_0

tax_table_cacao_ITS$Kingdom<- gsub("k__", "", tax_table_cacao_ITS$Kingdom)
tax_table_cacao_ITS$Phylum <- gsub("p__", "", tax_table_cacao_ITS$Phylum)
tax_table_cacao_ITS$Class <- gsub("c__", "", tax_table_cacao_ITS$Class)
tax_table_cacao_ITS$Order <- gsub("o__", "", tax_table_cacao_ITS$Order)
tax_table_cacao_ITS$Family <- gsub("f__", "", tax_table_cacao_ITS$Family)
tax_table_cacao_ITS$Genus <- gsub("g__", "", tax_table_cacao_ITS$Genus)
tax_table_cacao_ITS$Species <- gsub("s__", "", tax_table_cacao_ITS$Species)
tax_table_cacao_ITS$Fam_Gen <- paste(tax_table_cacao_ITS$Family,"_",tax_table_cacao_ITS$Genus,sep="")
tax_table_cacao_ITS$Gen_Spec <- paste(tax_table_cacao_ITS$Genus,"_",tax_table_cacao_ITS$Species,sep="")
tax_table_cacao_ITS$Fam_Gen_Spec <- paste(tax_table_cacao_ITS$Family,"_",tax_table_cacao_ITS$Genus,"_",tax_table_cacao_ITS$Species,sep="")

#View(tax_table_cacao_ITS)

row.names(tax_table_cacao_ITS) <- row.names(tax_table_cacao_ITS_0)
tax_table(cacao_ITS) <- as.matrix(tax_table_cacao_ITS)
#View(data.frame(tax_table(cacao_ITS)))
# tax_table_cacao_ITS_1 <- as.data.frame(tax_table(cacao_ITS))
# write.csv(tax_table_cacao_ITS_1, "tax_table_cacao_ITS_1.csv")

# Remove OTUs from negative controls. After looking at the original otu table, ITSNEG.ITS.a had most of the bad taxa, and had the same single bad taxa found in ITSNEG.ITS.b, so it was only necessary to remove taxa from ITSNEG.ITS.a as I did below
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

#Subsetting dataset
cacao_ITS = subset_samples(cacao_ITS_trim_0, Amplicon== "ITS")
cacao_ITS = subset_samples(cacao_ITS, Description != "ITSMOCK")

# View(data.frame(otu_table(cacao_ITS)))
# rarecurve(t(otu_table(cacao_ITS)), step=50, cex=0.5)
# cacao_ITS_otu <- data.frame(otu_table(cacao_ITS))
# write.csv(cacao_ITS_otu,"cacao_ITS_otu.csv")

# Prune samples
cacao_ITS_0 <- prune_samples(sample_sums(cacao_ITS) >= 50, cacao_ITS)
cacao_ITS_0_otu <- data.frame(otu_table(cacao_ITS_0))
#write.csv(cacao_ITS_0_otu,"cacao_ITS_0_otu.csv")

# View(data.frame(sample_data(cacao_ITS_0)))
# View(data.frame(otu_table(cacao_ITS_0)))
# rarecurve(t(otu_table(cacao_ITS)), step=50, cex=0.5)
# rarecurve(t(otu_table(cacao_ITS_0)), step=50, cex=0.5)

##########################################
###DESeq Comparison of phyllosphere [Urban Garden Center (left) vs Kualoa Ranch (right)]
##########################################
### Note: DESeq will show positive fold change for OTUs/ASVs that are higher in the comparison group (not the level you set in relevel) and negative foldchange for those that are enriched in the baseline (relevel group)
# View(sample_data(cacao_ITS_0))

#Subsetting data: location, phyllosphere samples
cacao_ITS_D0_0 = subset_samples(cacao_ITS_0, Ferm_Group == "D0")

#Remove empty cells due to subsetting; empty cells can cause issues later
cacao_ITS_D0 = filter_taxa(cacao_ITS_D0_0, function(x) sum(x) > 0, TRUE)

#Convert phyloseq to DESeq Data Set object (dds)
cacao_ITS_D0dds <- phyloseq_to_deseq2(cacao_ITS_D0, ~Location)

#Determine which level should the dataset be set as the REFERENCE sample class
cacao_ITS_D0dds$Location <- relevel(cacao_ITS_D0dds$Location, "Urban Garden Center" )

#Perform the contrast using standard and recognized parameters and tests
cacao_ITS_D0dds = DESeq(cacao_ITS_D0dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
#In this case, the contrast will be in between "Location_Kualoa.Ranch_vs_Urban.Garden.Center"
resultsNames(cacao_ITS_D0dds)

#Performing the final calculations and extracting the points
rescacao_ITS_D0dds = results(cacao_ITS_D0dds, cooksCutoff = TRUE)

###Contrast reports results such that positive fold change means the first level is enriched with a specific taxa, and a negative fold change means the second level is enriched with a specific taxa. For instance, a positive log fold change in the results below indicates enrichment in "soil", and a negative fold change indicates enrichment in "mat".
#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#rescacao_16S_D0dds = lfcShrink(cacao_16S_D0dds, contrast = c("Location", "Kualoa Ranch", "Urban Garden Center"), res=rescacao_16S_D0dds, type="apeglm") #this did not work
#rescacao_16S_D0dds <- lfcShrink(dds = cacao_16S_D0dds, coef = 2, type = "apeglm")#this worked
#rescacao_ITS_D0dds = lfcShrink(cacao_ITS_D0dds,contrast = c("Location","Kualoa Ranch","Urban Garden Center"),res=rescacao_ITS_D0dds)

#choosing an alpha of 0.05 I feel is pretty conservative, especially because DESeq is already conservative, but I still typically go with it.
alpha = 0.05

#Extract information from the DESeq object. The objects designated as "sigtab" have a p-value < or equal to the alpha set above. The objects names "notsig" are the results that have p-values > the alpha. These latter results can provide insight into common OTUs/ASVs.
#finding the differential abundant for each ASVs -- between the treatments
sigtab_D0 = rescacao_ITS_D0dds[which(rescacao_ITS_D0dds$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sigtab_D0 = cbind(as(sigtab_D0, "data.frame"), as(tax_table(cacao_ITS_D0)[rownames(sigtab_D0), ], "matrix"))

#View(sigtab_D0)
#write.csv(sigtab_D0, "sigtab_D0.csv")

#Find ASVs that are not significant; ASVs that are common among the treatments
notsigtab_D0 = rescacao_ITS_D0dds[which(rescacao_ITS_D0dds$padj > alpha), ]

#Bind taxa names to tables of not significant taxa
notsigtab_D0 = cbind(as(notsigtab_D0, "data.frame"), as(tax_table(cacao_ITS_D0)[rownames(notsigtab_D0), ], "matrix"))
#View(notsigtab_D0)
#Write.csv(notsigtab_D0, "notsigtab_D0.csv")

#This will provide a list of the phyla, so you can make sure all of your phyla are in the colors in "scale_colour_manual" below
sigtab_D0sub_Phyla <- unique(sigtab_D0sub$Phylum)
#View(sigtab_D0sub_Phyla)

#Remove anything that do not have a family or phylum taxonomy
sigtab_D0sub <- subset(sigtab_D0, Family!="NA" & Phylum!= "N/A")

#Write this to a table and make manual edits on the genus column. I'm sure there is an easy way to do this in R but time constraints...
####write.csv(sigtab_D0sub, "sigtab_D0sub_manual_edits.csv")

#Read edited file back in
sigtab_D0sub <- read.csv("sigtab_D0sub_manual_edits.csv")
#View(sigtab_D0sub)

#Plot the logfold changes
sigtab_D0subp <- ggplot(sigtab_D0sub, aes(x=log2FoldChange, y=reorder(Genus,desc(Genus)), color=Phylum)) + 
  geom_point(size=2, stroke=0.6) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25)) + 
  theme(legend.position="none") + 
  ggtitle("Phyllosphere fungi: Urban Garden vs. Kualoa Farm") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black")

#Plot and beautify by faceting based on phyla using manual colors
sigtab_D0subp2 <- sigtab_D0subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=4, face="italic")) + 
  theme(axis.title.y=element_blank()) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(text = element_text(size = 8)) + 
  #scale_color_manual(values = c("Ascomycota"="gold1", "Basidiomycota"="orchid2"))
  #scale_color_manual(values = c("Ascomycota"="darkorange", "Basidiomycota"="#7a81ff"))
  scale_colour_manual(values = c("Ascomycota"="#FD8D3C", "Basidiomycota"="#46AEA0"))
sigtab_D0subp2

ggsave("Figure2B-DESeq-phyllosphere.pdf", device="pdf", width=3.5, height=5)

# svg("sigtab_D0subp2.svg", width = 7.5, height = 11, pointsize=10)
# plot(sigtab_D0subp2)
# dev.off()


###################################################################
##### Comparison of phyllosphere vs ferment at Urban Garden Center
###################################################################

#Subsetting data: location
cacao_ITS_UG = subset_samples(cacao_ITS_0, Location == "Urban Garden Center")

#Remove empty cells due to subsetting; empty cells can cause issues later
cacao_ITS_UG_D0D1 = filter_taxa(cacao_ITS_UG, function(x) sum(x) > 0, TRUE)
#View(sample_data(cacao_ITS_0))
#View(otu_table(cacao_ITS_UG_D0D1))
#View(sample_data(cacao_ITS_UG_D0D1))

#Convert phyloseq to DESeq Data Set object (dds)
cacao_ITS_UG_D0D1dds <- phyloseq_to_deseq2(cacao_ITS_UG_D0D1, ~Ferm_Group)

#Determine which level should the dataset be set as the REFERENCE sample class
cacao_ITS_UG_D0D1dds$Ferm_Group <- relevel(cacao_ITS_UG_D0D1dds$Ferm_Group, "D0")#D0=phyllosphere

#Perform the contrast using standard and recognized parameters and tests
cacao_ITS_UG_D0D1dds = DESeq(cacao_ITS_UG_D0D1dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
resultsNames(cacao_ITS_UG_D0D1dds)

#Performing the final calculations and extracting the points
rescacao_ITS_UG_D0D1dds = results(cacao_ITS_UG_D0D1dds, cooksCutoff = TRUE)

#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#rescacao_ITS_UG_D0D1dds <- lfcShrink(dds = cacao_ITS_UG_D0D1dds, coef = 2, type = "apeglm")

#Set an alpha
alpha = 0.05

#Extract information from the DESeq object
sigtab_UG_D0D1 = rescacao_ITS_UG_D0D1dds[which(rescacao_ITS_UG_D0D1dds$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sigtab_UG_D0D1 = cbind(as(sigtab_UG_D0D1, "data.frame"), as(tax_table(cacao_ITS_UG_D0D1)[rownames(sigtab_UG_D0D1), ], "matrix"))

#View(sigtab_UG_D0D1)
#write.csv(sigtab_UG_D0D1, "sigtab_UG_D0D1.csv")

# notsigtab_UG_D0D1 = rescacao_ITS_UG_D0D1dds[which(rescacao_ITS_UG_D0D1dds$padj > alpha), ]
# notsigtab_UG_D0D1 = cbind(as(notsigtab_UG_D0D1, "data.frame"), as(tax_table(cacao_ITS_UG_D0D1)[rownames(notsigtab_UG_D0D1), ], "matrix"))
# head(notsigtab_UG_D0D1)
# write.csv(notsigtab_UG_D0D1, "notsigtab_UG_D0D1.csv")

sigtab_UG_D0D1_Phyla <- unique(sigtab_UG_D0D1$Phylum)
sigtab_UG_D0D1_Phyla

#Remove anything that do not have a family or phylum taxonomy
sigtab_UG_D0D1sub <- subset(sigtab_UG_D0D1, Family!="N/A" & Phylum != "N/A")
#View(sigtab_UG_D0D1sub)

#Write this to a table and make manual edits on the genus column. I'm sure there is an easy way to do this in R but time constraints...
####write.csv(sigtab_UG_D0D1sub, "sigtab_UG_D0D1sub_manual_edits.csv")

#Read edited file back in
sigtab_UG_D0D1sub <- read.csv("sigtab_UG_D0D1sub_manual_edits.csv")

#Plot the logfold changes
sigtab_UG_D0D1subp <- ggplot(sigtab_UG_D0D1sub, aes(x=log2FoldChange, y=reorder(Species,desc(Species)), color=Phylum)) + 
  geom_point(size=2, stroke=0.3) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25)) + 
  theme(legend.position="none") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black")

#Plot and beautify by faceting based on phyla using manual colors
sigtab_UG_D0D1subp2 <- sigtab_UG_D0D1subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=4, face="italic")) + 
  theme(axis.title.y=element_blank()) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(text = element_text(size = 8)) + 
  scale_colour_manual(values = c("Ascomycota"="#FD8D3C", "Basidiomycota"="#46AEA0"))
sigtab_UG_D0D1subp2

ggsave("Figure3B-DESeq-phyllosphere-ferment.pdf", device="pdf", width=3, height=2.5)


###################################################################
##### Comparison of phyllosphere vs ferment at Kualoa Ranch
###################################################################
#Subsetting data: location
cacao_ITS_KR = subset_samples(cacao_ITS_0, Location == "Kualoa Ranch")

#Remove empty cells due to subsetting; empty cells can cause issues later
cacao_ITS_KR_D0D1 = filter_taxa(cacao_ITS_KR, function(x) sum(x) > 0, TRUE)
#View(as.data.frame(sample_data(cacao_ITS_0)))
#View(as.data.frame(otu_table(cacao_ITS_KR_D0D1)))
#View(data.frame(sample_data(cacao_ITS_KR_D0D1)))

#Convert phyloseq to DESeq Data Set object (dds)
cacao_ITS_KR_D0D1dds <- phyloseq_to_deseq2(cacao_ITS_KR_D0D1, ~Ferm_Group)

#Determine which level should the dataset be set as the REFERENCE sample class
cacao_ITS_KR_D0D1dds$Ferm_Group <- relevel( cacao_ITS_KR_D0D1dds$Ferm_Group, "D0" )

#Perform the contrast using standard and recognized parameters and tests
cacao_ITS_KR_D0D1dds = DESeq(cacao_ITS_KR_D0D1dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
resultsNames(cacao_ITS_KR_D0D1dds)

#Check to see if the comparision conditions are present
rescacao_ITS_KR_D0D1dds = results(cacao_ITS_KR_D0D1dds, cooksCutoff = TRUE)

#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#rescacao_ITS_KR_D0D1dds = lfcShrink(cacao_ITS_KR_D0D1dds,contrast = c("Ferm_Group","D1","D0"),res=rescacao_ITS_KR_D0D1dds)

#Set an alpha
alpha = 0.05

#Extract information from the DESeq object
sigtab_KR_D0D1 = rescacao_ITS_KR_D0D1dds[which(rescacao_ITS_KR_D0D1dds$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sigtab_KR_D0D1 = cbind(as(sigtab_KR_D0D1, "data.frame"), as(tax_table(cacao_ITS_KR_D0D1)[rownames(sigtab_KR_D0D1), ], "matrix"))
#View(sigtab_KR_D0D1)
#write.csv(sigtab_KR_D0D1, "sigtab_KR_D0D1.csv")

# notsigtab_KR_D0D1 = rescacao_ITS_KR_D0D1dds[which(rescacao_ITS_KR_D0D1dds$padj > alpha), ]
# notsigtab_KR_D0D1 = cbind(as(notsigtab_KR_D0D1, "data.frame"), as(tax_table(cacao_ITS_KR_D0D1)[rownames(notsigtab_KR_D0D1), ], "matrix"))
# head(notsigtab_KR_D0D1)
# write.csv(notsigtab_KR_D0D1, "notsigtab_KR_D0D1.csv")

sigtab_KR_D0D1_Phyla <- unique(sigtab_KR_D0D1$Phylum)
sigtab_KR_D0D1_Phyla

sigtab_KR_D0D1sub <- subset(sigtab_KR_D0D1, Family!="N/A" & Phylum!= "N/A")
#View(sigtab_KR_D0D1sub)

#Write this to a table and make manual edits on the genus column. I'm sure there is an easy way to do this in R but time constraints...
write.csv(sigtab_KR_D0D1sub, "sigtab_KR_D0D1sub_manual_edits.csv")

#Read edited file back in
sigtab_KR_D0D1sub <- read.csv("sigtab_KR_D0D1sub_manual_edits.csv")

#Plot the logfold changes
sigtab_KR_D0D1subp <- ggplot(sigtab_KR_D0D1sub, aes(x=log2FoldChange, y=reorder(Species,desc(Species)), color=Phylum)) + 
  geom_point(size=2, stroke=0.3) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25)) + 
  theme(legend.position="none") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black")

#Plot and beautify by faceting based on phyla using manual colors
sigtab_KR_D0D1subp2 <- sigtab_KR_D0D1subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=4, face="italic")) + 
  theme(axis.title.y=element_blank()) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(text = element_text(size = 8)) + 
  scale_colour_manual(values = c("Ascomycota"="#FD8D3C", "Basidiomycota"="#46AEA0"))
sigtab_KR_D0D1subp2

ggsave("Figure3A-DESeq-phyllosphere-ferment.pdf", device="pdf", width=3.1, height=2.5)


#==============================================================================================
#Other comparisons not used in the manuscript
#==============================================================================================
###DESeq Day 1 WF vs UG
cacao_ITS_D1_0 = subset_samples(cacao_ITS_0, Ferm_Group == "D1")

cacao_ITS_D1 = filter_taxa(cacao_ITS_D1_0, function(x) sum(x) > 0, TRUE)

# View(as.data.frame(sample_data(cacao_ITS_D1)))

# #View(as.data.frame(otu_table(cacao_ITS_D1)))

# #View(data.frame(sample_data(cacao_ITS_D1)))

cacao_ITS_D1dds <- phyloseq_to_deseq2(cacao_ITS_D1, ~Location)
cacao_ITS_D1dds$Location <- relevel(cacao_ITS_D1dds$Location, "Urban Garden Center")

cacao_ITS_D1dds = DESeq(cacao_ITS_D1dds, test="Wald", fitType="parametric")
rescacao_ITS_D1dds = results(cacao_ITS_D1dds, cooksCutoff = TRUE)


rescacao_ITS_D1dds = lfcShrink(cacao_ITS_D1dds,contrast = c("Location","Kualoa Ranch","Urban Garden Center"),res=rescacao_ITS_D1dds)

alpha = 0.05

sigtab_D1 = rescacao_ITS_D1dds[which(rescacao_ITS_D1dds$padj <= alpha), ]

sigtab_D1 = cbind(as(sigtab_D1, "data.frame"), as(tax_table(cacao_ITS_D1)[rownames(sigtab_D1), ], "matrix"))
head(sigtab_D1)
write.csv(sigtab_D1, "sigtab_D1.csv")

notsigtab_D1 = rescacao_ITS_D1dds[which(rescacao_ITS_D1dds$padj > alpha), ]
notsigtab_D1 = cbind(as(notsigtab_D1, "data.frame"), as(tax_table(cacao_ITS_D1)[rownames(notsigtab_D1), ], "matrix"))
head(notsigtab_D1)
write.csv(notsigtab_D1, "notsigtab_D1.csv")



sigtab_D1_RWL <- read.csv("sigtab_D1_RWL.csv")
sigtab_D1sub <- subset(sigtab_D1_RWL, Family!="N/A" & Phylum != "N/A")
#View(sigtab_D1sub)

sigtab_D1sub_Phyla <- unique(sigtab_D1sub$Phylum)
#View(sigtab_D1sub_Phyla)


sigtab_D1subp = ggplot(sigtab_D1sub, aes(x=log2FoldChange, y=reorder(TaxaPlot2,desc(TaxaPlot2)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs D1 UG vs WF") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_D1subp 

sigtab_D1subp2 <- sigtab_D1subp  + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_D1subp2

svg("sigtab_D1subp2.svg", width = 7.5, height = 7.5, pointsize=12)
plot(sigtab_D1subp2)
dev.off()

##########
# #### Urban Garden Day 1 of the fermentation
cacao_ITS_UG = subset_samples(cacao_ITS_0, Location == "Urban Garden Center")
cacao_ITS_UG_1_0 <- subset_samples(cacao_ITS_UG, days_ferment == "1")
cacao_ITS_UG_1 = filter_taxa(cacao_ITS_UG_1_0, function(x) sum(x) > 0, TRUE)


cacao_ITS_UG_1dds <- phyloseq_to_deseq2(cacao_ITS_UG_1, ~tree_selection)
cacao_ITS_UG_1dds$tree_selection <- relevel( cacao_ITS_UG_1dds$tree_selection, "3" )

cacao_ITS_UG_1dds = DESeq(cacao_ITS_UG_1dds, test="Wald", fitType="parametric")
rescacao_ITS_UG_1dds = results(cacao_ITS_UG_1dds, cooksCutoff = TRUE)


rescacao_ITS_UG_1dds = lfcShrink(cacao_ITS_UG_1dds,contrast = c("tree_selection","5","3"),res=rescacao_ITS_UG_1dds)

alpha = 0.05

sigtab_UG_1 = rescacao_ITS_UG_1dds[which(rescacao_ITS_UG_1dds$padj <= alpha), ]
sigtab_UG_1 = cbind(as(sigtab_UG_1, "data.frame"), as(tax_table(cacao_ITS_UG_1)[rownames(sigtab_UG_1), ], "matrix"))
head(sigtab_UG_1)
#write.csv(sigtab_UG_1, "sig_UG_1.csv")

# notsigtab_UG_1 = rescacao_ITS_UG_1dds[which(rescacao_ITS_UG_1dds$padj > alpha), ]
# notsigtab_UG_1 = cbind(as(notsigtab_UG_1, "data.frame"), as(tax_table(cacao_ITS_UG_1)[rownames(notsigtab_UG_1), ], "matrix"))
# head(notsigtab_UG_1)
# write.csv(notsigtab_UG_1, "notsig_UG_1.csv")



#########
#### WF Day 1 of the fermentation

cacao_ITS_KR = subset_samples(cacao_ITS_0, Location == "Kualoa Ranch")
cacao_ITS_KR_1_0 <- subset_samples(cacao_ITS_KR, days_ferment == "1")
cacao_ITS_KR_1 = filter_taxa(cacao_ITS_KR_1_0, function(x) sum(x) > 0, TRUE)


cacao_ITS_KR_1dds <- phyloseq_to_deseq2(cacao_ITS_KR_1, ~tree_selection)
cacao_ITS_KR_1dds$tree_selection <- relevel( cacao_ITS_KR_1dds$tree_selection, "3" )

cacao_ITS_KR_1dds = DESeq(cacao_ITS_KR_1dds, test="Wald", fitType="parametric")
rescacao_ITS_KR_1dds = results(cacao_ITS_KR_1dds, cooksCutoff = TRUE)


rescacao_ITS_KR_1dds = lfcShrink(cacao_ITS_KR_1dds,contrast = c("tree_selection","5","3"),res=rescacao_ITS_KR_1dds)

alpha = 0.05

### DA-ASVs for WF Day 1
sigtab_KR_1 = rescacao_ITS_KR_1dds[which(rescacao_ITS_KR_1dds$padj <= alpha), ]
sigtab_KR_1 = cbind(as(sigtab_KR_1, "data.frame"), as(tax_table(cacao_ITS_KR_1)[rownames(sigtab_KR_1), ], "matrix"))
head(sigtab_KR_1)
write.csv(sigtab_KR_1, "sig_KR_1.csv")

notsigtab_KR_1 = rescacao_ITS_KR_1dds[which(rescacao_ITS_KR_1dds$padj > alpha), ]
notsigtab_KR_1 = cbind(as(notsigtab_KR_1, "data.frame"), as(tax_table(cacao_ITS_KR_1)[rownames(notsigtab_KR_1), ], "matrix"))
head(notsigtab_KR_1)
write.csv(notsigtab_KR_1, "notsig_KR_1.csv")



############################
#Selection 3: UG vs WF Day 0

cacao_ITS_S3 = subset_samples(cacao_ITS_0, tree_selection == "3")
cacao_ITS_S3_0_0 <- subset_samples(cacao_ITS_S3, days_ferment == "0")
cacao_ITS_S3_0 = filter_taxa(cacao_ITS_S3_0_0, function(x) sum(x) > 0, TRUE)

# #View(data.frame(sample_data(cacao_ITS_S3_0)))

cacao_ITS_S3_0dds <- phyloseq_to_deseq2(cacao_ITS_S3_0, ~Location)
cacao_ITS_S3_0dds$tree_selection <- relevel( cacao_ITS_S3_0dds$Location, "Urban Garden Center" )

cacao_ITS_S3_0dds = DESeq(cacao_ITS_S3_0dds, test="Wald", fitType="parametric")
rescacao_ITS_S3_0dds = results(cacao_ITS_S3_0dds, cooksCutoff = TRUE)
rescacao_ITS_S3_0dds = lfcShrink(cacao_ITS_S3_0dds,contrast = c("Location","Kualoa Ranch","Urban Garden Center"),res=rescacao_ITS_S3_0dds)

alpha = 0.05
sigtab_S3_0 = rescacao_ITS_S3_0dds[which(rescacao_ITS_S3_0dds$padj <= alpha), ]
sigtab_S3_0 = cbind(as(sigtab_S3_0, "data.frame"), as(tax_table(cacao_ITS_S3_0)[rownames(sigtab_S3_0), ], "matrix"))
head(sigtab_S3_0)
write.csv(sigtab_S3_0, "sig_S3_0.csv")

notsigtab_S3_0 = rescacao_ITS_S3_0dds[which(rescacao_ITS_S3_0dds$padj > alpha), ]
notsigtab_S3_0 = cbind(as(notsigtab_S3_0, "data.frame"), as(tax_table(cacao_ITS_S3_0)[rownames(notsigtab_S3_0), ], "matrix"))
head(notsigtab_S3_0)
write.csv(notsigtab_S3_0, "notsig_S3_0.csv")

#Selection 3: UG vs WF Day 1

cacao_ITS_S3 = subset_samples(cacao_ITS_0, tree_selection == "3")
cacao_ITS_S3_1_0 <- subset_samples(cacao_ITS_S3, days_ferment == "1")
cacao_ITS_S3_1 = filter_taxa(cacao_ITS_S3_1_0, function(x) sum(x) > 0, TRUE)

# #View(data.frame(sample_data(cacao_ITS_S3_1)))

cacao_ITS_S3_1dds <- phyloseq_to_deseq2(cacao_ITS_S3_1, ~Location)
cacao_ITS_S3_1dds$tree_selection <- relevel( cacao_ITS_S3_1dds$Location, "Urban Garden Center" )

cacao_ITS_S3_1dds = DESeq(cacao_ITS_S3_1dds, test="Wald", fitType="parametric")
rescacao_ITS_S3_1dds = results(cacao_ITS_S3_1dds, cooksCutoff = TRUE)
rescacao_ITS_S3_1dds = lfcShrink(cacao_ITS_S3_1dds,contrast = c("Location","Kualoa Ranch","Urban Garden Center"),res=rescacao_ITS_S3_1dds)

alpha = 0.05
sigtab_S3_1 = rescacao_ITS_S3_1dds[which(rescacao_ITS_S3_1dds$padj <= alpha), ]
sigtab_S3_1 = cbind(as(sigtab_S3_1, "data.frame"), as(tax_table(cacao_ITS_S3_1)[rownames(sigtab_S3_1), ], "matrix"))
head(sigtab_S3_1)
write.csv(sigtab_S3_1, "sig_S3_1.csv")

notsigtab_S3_1 = rescacao_ITS_S3_1dds[which(rescacao_ITS_S3_1dds$padj > alpha), ]
notsigtab_S3_1 = cbind(as(notsigtab_S3_1, "data.frame"), as(tax_table(cacao_ITS_S3_1)[rownames(notsigtab_S3_1), ], "matrix"))
head(notsigtab_S3_1)
write.csv(notsigtab_S3_1, "notsig_S3_1.csv")






###############
#Urban Garden day 0
cacao_ITS_UG = subset_samples(cacao_ITS_0, Location == "Urban Garden Center")
cacao_ITS_UG_0_0 <- subset_samples(cacao_ITS_UG, days_ferment == "0")
cacao_ITS_UG_0 = filter_taxa(cacao_ITS_UG_0_0, function(x) sum(x) > 0, TRUE)


# #View(as.data.frame(otu_table(cacao_ITS_UG_0)))

# #View(data.frame(sample_data(cacao_ITS_UG_0)))

cacao_ITS_UG_0dds <- phyloseq_to_deseq2(cacao_ITS_UG_0, ~tree_selection)
cacao_ITS_UG_0dds$tree_selection <- relevel( cacao_ITS_UG_0dds$tree_selection, "3" )

cacao_ITS_UG_0dds = DESeq(cacao_ITS_UG_0dds, test="Wald", fitType="parametric")
rescacao_ITS_UG_0dds = results(cacao_ITS_UG_0dds, cooksCutoff = TRUE)


rescacao_ITS_UG_0dds = lfcShrink(cacao_ITS_UG_0dds,contrast = c("tree_selection","5","3"),res=rescacao_ITS_UG_0dds)

alpha = 0.05

sigtab_UG_0 = rescacao_ITS_UG_0dds[which(rescacao_ITS_UG_0dds$padj <= alpha), ]
sigtab_UG_0 = cbind(as(sigtab_UG_0, "data.frame"), as(tax_table(cacao_ITS_UG_0)[rownames(sigtab_UG_0), ], "matrix"))
head(sigtab_UG_0)
write.csv(sigtab_UG_0, "sig_UG_0.csv")

notsigtab_UG_0 = rescacao_ITS_UG_0dds[which(rescacao_ITS_UG_0dds$padj > alpha), ]
notsigtab_UG_0 = cbind(as(notsigtab_UG_0, "data.frame"), as(tax_table(cacao_ITS_UG_0)[rownames(notsigtab_UG_0), ], "matrix"))
head(notsigtab_UG_0)
write.csv(notsigtab_UG_0, "notsig_UG_0.csv")








###############
#Kualoa Ranch day 0

cacao_ITS_KR = subset_samples(cacao_ITS_0, Location == "Kualoa Ranch")
cacao_ITS_KR_0_0 <- subset_samples(cacao_ITS_KR, days_ferment == "0")
cacao_ITS_KR_0 = filter_taxa(cacao_ITS_KR_0_0, function(x) sum(x) > 0, TRUE)

#View(data.frame(sample_data(cacao_ITS_KR_0)))

cacao_ITS_KR_0dds <- phyloseq_to_deseq2(cacao_ITS_KR_0, ~tree_selection)
cacao_ITS_KR_0dds$tree_selection <- relevel( cacao_ITS_KR_0dds$tree_selection, "3" )

cacao_ITS_KR_0dds = DESeq(cacao_ITS_KR_0dds, test="Wald", fitType="parametric")
rescacao_ITS_KR_0dds = results(cacao_ITS_KR_0dds, cooksCutoff = TRUE)
rescacao_ITS_KR_0dds = lfcShrink(cacao_ITS_KR_0dds,contrast = c("tree_selection","5","3"),res=rescacao_ITS_KR_0dds)

alpha = 0.05
sigtab_KR_0 = rescacao_ITS_KR_0dds[which(rescacao_ITS_KR_0dds$padj <= alpha), ]
sigtab_KR_0 = cbind(as(sigtab_KR_0, "data.frame"), as(tax_table(cacao_ITS_KR_0)[rownames(sigtab_KR_0), ], "matrix"))
head(sigtab_KR_0)
write.csv(sigtab_KR_0, "sig_KR_0.csv")

notsigtab_KR_0 = rescacao_ITS_KR_0dds[which(rescacao_ITS_KR_0dds$padj > alpha), ]
notsigtab_KR_0 = cbind(as(notsigtab_KR_0, "data.frame"), as(tax_table(cacao_ITS_KR_0)[rownames(notsigtab_KR_0), ], "matrix"))
head(notsigtab_KR_0)
write.csv(notsigtab_KR_0, "notsig_KR_0.csv")



###############
###DA OTUs Urban Garden Day 0
sigtab_UG_0_2 <- sigtab_UG_0

sigtab_UG_0_2$abslog2FC <- abs(sigtab_UG_0_2$log2FoldChange)
sigtab_UG_0_2 <- subset(sigtab_UG_0_2, abslog2FC >= 0.945)
# #View(sigtab_UG_0_2)

sigtab_UG_0_3 <- sigtab_UG_0_2
sigtab_UG_0_3sub <- subset(sigtab_UG_0_3, Family!="N/A" & Phylum != "N/A")
#View(sigtab_UG_0_3sub)

sigtab_UG_0_3sub_Phyla <- unique(sigtab_UG_0_3sub$Phylum)
#View(sigtab_UG_0_3sub_Phyla)

#View(sigtab_UG_0_2)

sigtab_UG_0_3subp = ggplot(sigtab_UG_0_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Urban Garden Day 0: Selection 3 vs 5") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_UG_0_3subp

sigtab_UG_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())


###DA OTUs Urban Garden Day 1
sigtab_UG_1_2 <- sigtab_UG_1

sigtab_UG_1_2$abslog2FC <- abs(sigtab_UG_1_2$log2FoldChange)
sigtab_UG_1_2 <- subset(sigtab_UG_1_2, abslog2FC >= 0.945)
# #View(sigtab_UG_1_2)



sigtab_UG_1_3 <- sigtab_UG_1_2
sigtab_UG_1_3sub <- subset(sigtab_UG_1_3, Family!="N/A" & Phylum != "N/A")
#View(sigtab_UG_1_3sub)

sigtab_UG_1_3sub_Phyla <- unique(sigtab_UG_1_3sub$Phylum)
#View(sigtab_UG_1_3sub_Phyla)

#View(sigtab_UG_1_2)

sigtab_UG_1_3subp = ggplot(sigtab_UG_1_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Urban Garden Day 1: Selection 3 vs 5") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_UG_1_3subp

sigtab_UG_1_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())



### DA Kualoa Ranch (WF) Day 0

sigtab_KR_0_2 <- sigtab_KR_0

sigtab_KR_0_2$abslog2FC <- abs(sigtab_KR_0_2$log2FoldChange)
sigtab_KR_0_2 <- subset(sigtab_KR_0_2, abslog2FC >= 0.945)

sigtab_KR_0_3 <- sigtab_KR_0_2

sigtab_KR_0_3sub <- subset(sigtab_KR_0_3, Family!="N/A" & Phylum != "N/A")
#View(sigtab_KR_0_3sub)

sigtab_KR_0_3sub_Phyla <- unique(sigtab_KR_0_3sub$Phylum)
#View(sigtab_KR_0_3sub_Phyla)


sigtab_KR_0_3subp = ggplot(sigtab_KR_0_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Kualoa Ranch Day 0: Selection 3 vs 5") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_KR_0_3subp

sigtab_KR_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())


##DESEQ_Summary Selections compared at each Site at Day 0

sigtab_UG_0_3subp
sigtab_KR_0_3subp




#Selection 5: UG vs WF Day 0

cacao_ITS_S5 = subset_samples(cacao_ITS_0, tree_selection == "5")
cacao_ITS_S5_0_0 <- subset_samples(cacao_ITS_S5, days_ferment == "0")
cacao_ITS_S5_0 = filter_taxa(cacao_ITS_S5_0_0, function(x) sum(x) > 0, TRUE)

#View(data.frame(sample_data(cacao_ITS_S5_0)))

cacao_ITS_S5_0dds <- phyloseq_to_deseq2(cacao_ITS_S5_0, ~Location)
cacao_ITS_S5_0dds$tree_selection <- relevel( cacao_ITS_S5_0dds$Location, "Urban Garden Center" )

cacao_ITS_S5_0dds = DESeq(cacao_ITS_S5_0dds, test="Wald", fitType="parametric")
rescacao_ITS_S5_0dds = results(cacao_ITS_S5_0dds, cooksCutoff = TRUE)
rescacao_ITS_S5_0dds = lfcShrink(cacao_ITS_S5_0dds,contrast = c("Location","Kualoa Ranch","Urban Garden Center"),res=rescacao_ITS_S5_0dds)

alpha = 0.05
sigtab_S5_0 = rescacao_ITS_S5_0dds[which(rescacao_ITS_S5_0dds$padj <= alpha), ]
sigtab_S5_0 = cbind(as(sigtab_S5_0, "data.frame"), as(tax_table(cacao_ITS_S5_0)[rownames(sigtab_S5_0), ], "matrix"))
head(sigtab_S5_0)
write.csv(sigtab_S5_0, "sig_S5_0.csv")

notsigtab_S5_0 = rescacao_ITS_S5_0dds[which(rescacao_ITS_S5_0dds$padj > alpha), ]
notsigtab_S5_0 = cbind(as(notsigtab_S5_0, "data.frame"), as(tax_table(cacao_ITS_S5_0)[rownames(notsigtab_S5_0), ], "matrix"))
head(notsigtab_S5_0)
write.csv(notsigtab_S5_0, "notsig_S5_0.csv")


#Selection 5: UG vs WF Day 1

cacao_ITS_S5 = subset_samples(cacao_ITS_0, tree_selection == "5")
cacao_ITS_S5_1_0 <- subset_samples(cacao_ITS_S5, days_ferment == "1")
cacao_ITS_S5_1 = filter_taxa(cacao_ITS_S5_1_0, function(x) sum(x) > 0, TRUE)

# #View(data.frame(sample_data(cacao_ITS_S5_1)))

cacao_ITS_S5_1dds <- phyloseq_to_deseq2(cacao_ITS_S5_1, ~Location)
cacao_ITS_S5_1dds$tree_selection <- relevel( cacao_ITS_S5_1dds$Location, "Urban Garden Center" )

cacao_ITS_S5_1dds = DESeq(cacao_ITS_S5_1dds, test="Wald", fitType="parametric")
rescacao_ITS_S5_1dds = results(cacao_ITS_S5_1dds, cooksCutoff = TRUE)
rescacao_ITS_S5_1dds = lfcShrink(cacao_ITS_S5_1dds,contrast = c("Location","Kualoa Ranch","Urban Garden Center"),res=rescacao_ITS_S5_1dds)

alpha = 0.05
sigtab_S5_1 = rescacao_ITS_S5_1dds[which(rescacao_ITS_S5_1dds$padj <= alpha), ]
sigtab_S5_1 = cbind(as(sigtab_S5_1, "data.frame"), as(tax_table(cacao_ITS_S5_1)[rownames(sigtab_S5_1), ], "matrix"))
head(sigtab_S5_1)
write.csv(sigtab_S5_1, "sig_S5_1.csv")

notsigtab_S5_1 = rescacao_ITS_S5_1dds[which(rescacao_ITS_S5_1dds$padj > alpha), ]
notsigtab_S5_1 = cbind(as(notsigtab_S5_1, "data.frame"), as(tax_table(cacao_ITS_S5_1)[rownames(notsigtab_S5_1), ], "matrix"))
head(notsigtab_S5_1)
write.csv(notsigtab_S5_1, "notsig_S5_1.csv")


### DA Selection 3: UG vs WF Day 0

sigtab_S3_0_2 <- sigtab_S3_0

sigtab_S3_0_2$abslog2FC <- abs(sigtab_S3_0_2$log2FoldChange)
sigtab_S3_0_2 <- subset(sigtab_S3_0_2, abslog2FC >= 0.945)
#View(sigtab_S3_0_2)

sigtab_S3_0_3 <- sigtab_S3_0_2


sigtab_S3_0_3sub <- subset(sigtab_S3_0_3, Family!="N/A" & Phylum != "N/A")
#View(sigtab_S3_0_3sub)

sigtab_S3_0_3sub_Phyla <- unique(sigtab_S3_0_3sub$Phylum)
#View(sigtab_S3_0_3sub_Phyla)


sigtab_S3_0_3subp = ggplot(sigtab_S3_0_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Selection 3 Day 0: Urban Garden vs Kualoa Ranch") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_S3_0_3subp

sigtab_S3_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())


### DA Selection 3: UG vs WF Day 1

sigtab_S3_1_2 <- sigtab_S3_1

sigtab_S3_1_2$abslog2FC <- abs(sigtab_S3_1_2$log2FoldChange)
sigtab_S3_1_2 <- subset(sigtab_S3_1_2, abslog2FC >= 0.945)
#View(sigtab_S3_1_2)

sigtab_S3_1_3 <- sigtab_S3_1_2

sigtab_S3_1_3sub <- subset(sigtab_S3_1_3, Family!="N/A" & Phylum != "N/A")
#View(sigtab_S3_1_3sub)

sigtab_S3_1_3sub_Phyla <- unique(sigtab_S3_1_3sub$Phylum)
#View(sigtab_S3_1_3sub_Phyla)


sigtab_S3_1_3subp = ggplot(sigtab_S3_1_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Selection 3 Day 1: Urban Garden vs Kualoa Ranch") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_S3_1_3subp

sigtab_S3_1_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())


### DA Selection 5: UG vs WF Day 0

sigtab_S5_0_2 <- sigtab_S5_0

sigtab_S5_0_2$abslog2FC <- abs(sigtab_S5_0_2$log2FoldChange)
sigtab_S5_0_2 <- subset(sigtab_S5_0_2, abslog2FC >= 0.945)

sigtab_S5_0_3 <- sigtab_S5_0_2

sigtab_S5_0_3sub <- subset(sigtab_S5_0_3, Family!="N/A" & Phylum != "N/A")
#View(sigtab_S5_0_3sub)

sigtab_S5_0_3sub_Phyla <- unique(sigtab_S5_0_3sub$Phylum)
#View(sigtab_S5_0_3sub_Phyla)


sigtab_S5_0_3subp = ggplot(sigtab_S5_0_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Selection 5 Day 0: Urban Garden vs Kualoa Ranch") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_S5_0_3subp

sigtab_S5_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())



### DA Selection 5: UG vs WF Day 1

sigtab_S5_1_2 <- sigtab_S5_1

sigtab_S5_1_2$abslog2FC <- abs(sigtab_S5_1_2$log2FoldChange)
sigtab_S5_1_2 <- subset(sigtab_S5_1_2, abslog2FC >= 0.945)


sigtab_S5_1_3 <- sigtab_S5_1_2

sigtab_S5_1_3sub <- subset(sigtab_S5_1_3, Family!="N/A" & Phylum != "N/A")
#View(sigtab_S5_1_3sub)

sigtab_S5_1_3sub_Phyla <- unique(sigtab_S5_1_3sub$Phylum)
#View(sigtab_S5_1_3sub_Phyla)


sigtab_S5_1_3subp = ggplot(sigtab_S5_1_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Selection 5 Day 1: Urban Garden vs Kualoa Ranch") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_S5_1_3subp

sigtab_S5_1_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())


##Summary Plots

sigtab_UG_0_3plot <- sigtab_UG_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_UG_1_3plot <- sigtab_UG_1_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_KR_0_3plot <- sigtab_KR_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_S3_0_3plot <- sigtab_S3_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_S5_0_3plot <- sigtab_S5_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_S3_1_3plot <- sigtab_S3_1_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_S5_1_3plot <- sigtab_S5_1_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_UG_0_3plot
sigtab_UG_1_3plot
sigtab_KR_0_3plot
sigtab_S3_0_3plot
sigtab_S5_0_3plot
sigtab_S3_1_3plot
sigtab_S5_1_3plot


svg("sigtab_UG_0_3plot.svg", width = 7.5, height = 3, pointsize=12)
plot(sigtab_UG_0_3plot)
dev.off()

svg("sigtab_UG_1_3plot.svg", width = 7.5, height = 3, pointsize=12)
plot(sigtab_UG_1_3plot)
dev.off()


svg("sigtab_KR_0_3plot.svg", width = 7.5, height = 3, pointsize=12)
plot(sigtab_KR_0_3plot)
dev.off()

svg("sigtab_S3_0_3plot.svg", width = 7.5, height = 7.5, pointsize=12)
plot(sigtab_S3_0_3plot)
dev.off()

svg("sigtab_S5_0_3plot.svg", width = 7.5, height = 11, pointsize=12)
plot(sigtab_S5_0_3plot)
dev.off()

svg("sigtab_S3_1_3plot.svg", width = 7.5, height = 7.5, pointsize=12)
plot(sigtab_S3_1_3plot)
dev.off()

svg("sigtab_S5_1_3plot.svg", width = 7.5, height = 7.5, pointsize=12)
plot(sigtab_S5_1_3plot)
dev.off()