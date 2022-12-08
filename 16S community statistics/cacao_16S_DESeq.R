#set working directory
setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1--pvx1vLSBijmxHMCZPrxKlyFBLp2jqF/Rick/cacao_combined/cacao_physeq_full_analysis/cacao_16S_ps")

#load libraries
library(phyloseq)
library(DESeq2)
library(tidyverse)
theme_set(theme_bw())

#Import a phyloseq object (biom should have taxonomy attached), start with UNRAREFIED otu table
cacao_16S = import_biom("cacao_prok.biom", treefilename="cacao_rooted-tree-filtered2-prok.nwk", refseqfilename="cacao_prok-only.fasta")

#Rename columns 
colnames(tax_table(cacao_16S)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Save otu table (optional) -- this will save a simple OTU table without taxonomy
#write.csv(as.data.frame(otu_table(cacao_16S)),"otu_cacao_16S_initial.csv")

#Read in metadata; standard QIIME2 metadata table can be used here
cacao_16S_metadata <- read.csv("cacao_metadata_combined.csv")

#Create another row for "SampleID", create a dataframe for the metadata
row.names(cacao_16S_metadata) <- cacao_16S_metadata$X.SampleID
cacao_16S_metadata$SampleID <- cacao_16S_metadata$X.SampleID
sample_data(cacao_16S) <- cacao_16S_metadata
#View(data.frame(sample_data(cacao_16S)))

#De-QIIME-ify the taxa table -- this will separate taxonomic ranks into each separate columns. This _0 table is necessary downstream.
tax_table_cacao_16S_0 <- as.data.frame(tax_table(cacao_16S))
#write.csv(tax_table_cacao_16S_0, "tax_table_cacao_16S_0.csv")

#Make a copy of the table
tax_table_cacao_16S <- tax_table_cacao_16S_0

#Renaming the taxonomy if not standard (this may not be necessary depending on the taxonomic database used)
# tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Acidobacteriota", "Acidobacteria", x)}))
# tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Actinobacteriota", "Actinobacteria", x)}))
# tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Armatimonadota", "Armatimonadetes", x)}))
# tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Bacteroidota", "Bacteroidetes", x)}))
# tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Gemmatimonadota", "Gemmatimonadetes", x)}))
# tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Halobacterota", "Halobacteria", x)}))
# tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Planctomycetota", "Planctomycetes", x)}))
# tax_table_cacao_16S <- data.frame(lapply(tax_table_cacao_16S, function(x) {gsub("Verrucomicrobiota", "Verrucomicrobia", x)}))

#Remove taxonomic notations
tax_table_cacao_16S$Kingdom<- gsub("d__", "", tax_table_cacao_16S$Kingdom)
tax_table_cacao_16S$Phylum <- gsub("p__", "", tax_table_cacao_16S$Phylum)
tax_table_cacao_16S$Class <- gsub("c__", "", tax_table_cacao_16S$Class)
tax_table_cacao_16S$Order <- gsub("o__", "", tax_table_cacao_16S$Order)
tax_table_cacao_16S$Family <- gsub("f__", "", tax_table_cacao_16S$Family)
tax_table_cacao_16S$Genus <- gsub("g__", "", tax_table_cacao_16S$Genus)
tax_table_cacao_16S$Species <- gsub("s__", "", tax_table_cacao_16S$Species)
#tax_table_cacao_16S$Fam_Gen <- paste(tax_table_cacao_16S$Family,"_",tax_table_cacao_16S$Genus,sep="")

#View(tax_table_cacao_16S_0)
#View(tax_table_cacao_16S)

row.names(tax_table_cacao_16S) <- row.names(tax_table_cacao_16S_0)
tax_table(cacao_16S) <- as.matrix(tax_table_cacao_16S)
# ##View(data.frame(tax_table(cacao_16S)))
# tax_table_cacao_16S_1 <- as.data.frame(tax_table(cacao_16S))
# write.csv(tax_table_cacao_16S_1, "tax_table_cacao_16S_1.csv")

# can use for subsetting taxa: cacao_16S_sub0 = subset_taxa(cacao_16S, Kingdom=="Bacteria")
# Should we examine archaea?
cacao_16S = subset_samples(cacao_16S, Amplicon== "16S")
cacao_16S = subset_samples(cacao_16S, Description != "16SMOCK")
cacao_16S = subset_samples(cacao_16S, SampleID != "02.2.5.1.16S.a")

# View(data.frame(otu_table(cacao_16S)))
# rarecurve(t(otu_table(cacao_16S)), step=50, cex=0.5)

# Prune samples
cacao_16S_0 <- prune_samples(sample_sums(cacao_16S) >= 50, cacao_16S)
# #View(data.frame(sample_data(cacao_16S_0)))
# #View(data.frame(otu_table(cacao_16S_0)))
# rarecurve(t(otu_table(cacao_16S_0)), step=50, cex=0.5)
# cacao_16S_0_otu <- data.frame(otu_table(cacao_16S_0))
# write.csv(cacao_16S_0_otu,"cacao_16S_0_otu.csv")

##########################################
###DESeq Comparison of phyllosphere [Urban Garden Center (left) vs Kualoa Ranch (right)]
##########################################
#View(sample_data(cacao_16S_0))

#Subsetting data: location, phyllosphere samples
cacao_16S_D0_0 = subset_samples(cacao_16S_0, days_ferment == "0")

#Remove empty cells due to subsetting; empty cells can cause issues later
cacao_16S_D0 = filter_taxa(cacao_16S_D0_0, function(x) sum(x) > 0, TRUE)

#An error may occur later on when running DESeq because all OTUs have one 0. Replacing all 0's with 1's will solve this issue. If error does not occur, skip this step.
# Can replace zero with 1 in otu table using following code
# otu_table(cacao_16S_UG_0)[otu_table(cacao_16S_UG_0)==0] <- 1

#Convert phyloseq to DESeq Data Set object (dds)
cacao_16S_D0dds <- phyloseq_to_deseq2(cacao_16S_D0, ~Location)

#Determine which level should the dataset be set as the REFERENCE sample class
cacao_16S_D0dds$Location <- relevel(cacao_16S_D0dds$Location, "Urban Garden Center" )

#Perform the contrast using standard and recognized parameters and tests
cacao_16S_D0dds = DESeq(cacao_16S_D0dds, test="Wald", fitType="parametric")

#Check to see if the comparision conditions are present
#In this case, the contrast will be in between "Location_Kualoa.Ranch_vs_Urban.Garden.Center"
resultsNames(cacao_16S_D0dds)

#Performing the final calculations and extracting the points
rescacao_16S_D0dds = results(cacao_16S_D0dds, cooksCutoff = TRUE)

###Contrast reports results such that positive fold change means the first level is enriched with a specific taxa, and a negative fold change means the second level is enriched with a specific taxa. For instance, a positive log fold change in the results below indicates enrichment in "soil", and a negative fold change indicates enrichment in "mat".
#Reduce over estimation of fold changes in graphical format (makes the dots closer to better show on a figure)
#rescacao_16S_D0dds = lfcShrink(cacao_16S_D0dds, contrast = c("Location", "Kualoa Ranch", "Urban Garden Center"), res=rescacao_16S_D0dds, type="apeglm") #this did not work
#rescacao_16S_D0dds <- lfcShrink(dds = cacao_16S_D0dds, coef = 2, type = "apeglm")#this worked

#choosing an alpha of 0.05 I feel is pretty conservative, especially because DESeq is already conservative, but I still typically go with it.
alpha = 0.05

#Extract information from the DESeq object. The objects designated as "sigtab" have a p-value < or equal to the alpha set above. The objects names "notsig" are the results that have p-values > the alpha. These latter results can provide insight into common OTUs/ASVs.
#finding the differential abundant for each ASVs -- between the treatments
sigtab_D0 = rescacao_16S_D0dds[which(rescacao_16S_D0dds$padj <= alpha), ]

#Bind taxa names to tables of significant taxa
sigtab_D0 = cbind(as(sigtab_D0, "data.frame"), as(tax_table(cacao_16S_D0)[rownames(sigtab_D0), ], "matrix"))

#View(sigtab_D0)
#write.csv(sigtab_D0, "sig_D0.csv")

#Find ASVs that are not significant; ASVs that are common among the treatments
notsigtab_D0 = rescacao_16S_D0dds[which(rescacao_16S_D0dds$padj > alpha), ]

#Bind taxa names to tables of not significant taxa
notsigtab_D0 = cbind(as(notsigtab_D0, "data.frame"), as(tax_table(cacao_16S_D0)[rownames(notsigtab_D0), ], "matrix"))
#View(notsigtab_D0)
# write.csv(notsigtab_D0, "notsig_D0.csv")

#This will provide a list of the phyla, so you can make sure all of your phyla are in the colors in "scale_colour_manual" below
sigtab_D0_Phyla <- unique(sigtab_D0$Phylum)
sigtab_D0_Phyla

#Remove anything that do not have a family or phylum taxonomy
sigtab_D0sub <- subset(sigtab_D0, Family!="N/A" & Phylum != "N/A")

#Rename items in the Genus column for formatting. Each dataset will vary but this list below is a standard list; add to it as needed.
sigtab_D0sub$Genus <- gsub("uncultured", "unidentified bacterium", sigtab_D0sub$Genus)
sigtab_D0sub$Genus <- gsub("NA", "unidentified bacterium", sigtab_D0sub$Genus)
sigtab_D0sub$Genus <- gsub("Burkholderia-Caballeronia-Paraburkholderia", "Para/Burkholderia-Caballeronia", sigtab_D0sub$Genus)
sigtab_D0sub$Genus <- gsub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Allo/Neo/Para-Rhizobium", sigtab_D0sub$Genus)
sigtab_D0sub$Genus <- gsub("Methylobacterium-Methylorubrum", "Methylo-bacterium/-rubrum", sigtab_D0sub$Genus)
#View(sigtab_D0sub)

#Write this to a table and make manual edits on the genus column. I'm sure there is an easy way to do this in R but time constraints...
write.csv(sigtab_D0sub, "sigtab_D0sub_manual_edits.csv")

#Read edited file back in
sigtab_D0sub <- read.csv("sigtab_D0sub_manual_edits.csv")


#Plot the logfold changes
sigtab_D0subp = ggplot(sigtab_D0sub, aes(x=log2FoldChange, y=reorder(Genus,desc(Genus)), color=Phylum)) +
  geom_point(size=2, stroke = 0.6) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.25)) + 
  theme(legend.position="none") + 
  ggtitle("Phyllosphere bacteria: Urban Garden vs. Kualoa Farm") + 
  geom_vline(xintercept = 0, linetype = "solid", color = "black")

#Plot and beautify by faceting based on phyla using manual colors
sigtab_D0subp2 <- sigtab_D0subp + 
  facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + 
  scale_y_discrete(position = "right") + 
  theme(strip.text.y.left = element_text(angle = 0)) + 
  theme(axis.text.y=element_text(size=6, face="italic")) + 
  theme(axis.title.y=element_blank()) + 
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(text = element_text(size = 8)) + 
  scale_colour_manual(values = c("Acidobacteriota" = "red1", 
                                 "Actinobacteriota" = "darkorange", 
                                 "Armatimonadota" = "darkorange", 
                                 "Bacteroidota" = "gold1", 
                                 "Bdellovibrionota" = "#c5e908", 
                                 "Chloroflexi" = "#a4de02", 
                                 "Deinococcota" = "#86dc3d",
                                 "Firmicutes" = "seagreen", 
                                 "Myxococcota" = "darkturquoise", 
                                 "Planctomycetes" = "deepskyblue",
                                 "Proteobacteria" = "orchid2", 
                                 "Verrucomicrobiota" = "#7a81ff"))
sigtab_D0subp2

ggsave("Figure2A-DESeq-phyllosphere.pdf", device="pdf", width=4, height=5)

# svg("sigtab_D0subp2.svg", width = 6, height = 8, pointsize=12)
# plot(sigtab_D0subp2)
# dev.off()





















##########################################
###trying several other comparisons...did not use these in the manuscript
##########################################
#Urban Garden day 0

cacao_16S_UG = subset_samples(cacao_16S_0, Location == "Urban Garden Center")
cacao_16S_UG_0_0 <- subset_samples(cacao_16S_UG, days_ferment == "0")
cacao_16S_UG_0 = filter_taxa(cacao_16S_UG_0_0, function(x) sum(x) > 0, TRUE)

# Can replace zero with 1 in otu table using following code
# otu_table(cacao_16S_UG_0)[otu_table(cacao_16S_UG_0)==0] <- 1

# #View(as.data.frame(otu_table(cacao_16S_UG_0)))

cacao_16S_UG_0dds <- phyloseq_to_deseq2(cacao_16S_UG_0, ~tree_selection)
cacao_16S_UG_0dds$tree_selection <- relevel( cacao_16S_UG_0dds$tree_selection, "3" )

cacao_16S_UG_0dds = DESeq(cacao_16S_UG_0dds, test="Wald", fitType="parametric")
rescacao_16S_UG_0dds = results(cacao_16S_UG_0dds, cooksCutoff = TRUE)


rescacao_16S_UG_0dds = lfcShrink(cacao_16S_UG_0dds,contrast = c("tree_selection","5","3"),res=rescacao_16S_UG_0dds)

alpha = 0.05

sigtab_UG_0 = rescacao_16S_UG_0dds[which(rescacao_16S_UG_0dds$padj <= alpha), ]
sigtab_UG_0 = cbind(as(sigtab_UG_0, "data.frame"), as(tax_table(cacao_16S_UG_0)[rownames(sigtab_UG_0), ], "matrix"))
head(sigtab_UG_0)
write.csv(sigtab_UG_0, "sig_UG_0.csv")

notsigtab_UG_0 = rescacao_16S_UG_0dds[which(rescacao_16S_UG_0dds$padj > alpha), ]
notsigtab_UG_0 = cbind(as(notsigtab_UG_0, "data.frame"), as(tax_table(cacao_16S_UG_0)[rownames(notsigtab_UG_0), ], "matrix"))
head(notsigtab_UG_0)
write.csv(notsigtab_UG_0, "notsig_UG_0.csv")


#DESeq Waialua Farm day 0

cacao_16S_WF = subset_samples(cacao_16S_0, Location == "Waialua Farm")
cacao_16S_WF_0_0 <- subset_samples(cacao_16S_WF, days_ferment == "0")
cacao_16S_WF_0 = filter_taxa(cacao_16S_WF_0_0, function(x) sum(x) > 0, TRUE)

# #View(data.frame(sample_data(cacao_16S_WF_0)))

cacao_16S_WF_0dds <- phyloseq_to_deseq2(cacao_16S_WF_0, ~tree_selection)
cacao_16S_WF_0dds$tree_selection <- relevel( cacao_16S_WF_0dds$tree_selection, "3" )

cacao_16S_WF_0dds = DESeq(cacao_16S_WF_0dds, test="Wald", fitType="parametric")
rescacao_16S_WF_0dds = results(cacao_16S_WF_0dds, cooksCutoff = TRUE)
rescacao_16S_WF_0dds = lfcShrink(cacao_16S_WF_0dds,contrast = c("tree_selection","5","3"),res=rescacao_16S_WF_0dds)


alpha = 0.05
sigtab_WF_0 = rescacao_16S_WF_0dds[which(rescacao_16S_WF_0dds$padj <= alpha), ]
sigtab_WF_0 = cbind(as(sigtab_WF_0, "data.frame"), as(tax_table(cacao_16S_WF_0)[rownames(sigtab_WF_0), ], "matrix"))
head(sigtab_WF_0)
write.csv(sigtab_WF_0, "sig_WF_0.csv")

notsigtab_WF_0 = rescacao_16S_WF_0dds[which(rescacao_16S_WF_0dds$padj > alpha), ]
notsigtab_WF_0 = cbind(as(notsigtab_WF_0, "data.frame"), as(tax_table(cacao_16S_WF_0)[rownames(notsigtab_WF_0), ], "matrix"))
head(notsigtab_WF_0)
write.csv(notsigtab_WF_0, "notsig_WF_0.csv")

###DA OTUs Urban Garden Day 0 

sigtab_UG_0_2 <- sigtab_UG_0

sigtab_UG_0_2$abslog2FC <- abs(sigtab_UG_0_2$log2FoldChange)
sigtab_UG_0_2 <- subset(sigtab_UG_0_2, abslog2FC >= 0.945)


sigtab_UG_0_3 <- sigtab_UG_0_2
sigtab_UG_0_3sub <- subset(sigtab_UG_0_3, Family!="N/A" & Phylum != "N/A")
#View(sigtab_UG_0_3sub)

str(sigtab_UG_0_3sub)

sigtab_UG_0_3sub_Phyla <- unique(sigtab_UG_0_3sub$Phylum)
# #View(sigtab_UG_0_3sub_Phyla)

sigtab_UG_0_3subp = ggplot(sigtab_UG_0_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Urban Garden Day 0: Selection 3 vs 5") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_UG_0_3subp

sigtab_UG_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

#Note: color values have duplicate colors that are not duplicated in the graphs, but if this code is used elsewhere, one should double check: "Chloroflexi" = "green", "Deinococcota" = "green"

sigtab_UG_0_3subp2 <- sigtab_UG_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank()) + theme(panel.spacing = unit(0.1, "lines")) + theme(text = element_text(size = 8)) + scale_colour_manual(values = c("Acidobacteria" = "red1", "Actinobacteria" = "tomato", "Armatimonadetes" = "darkorange", "Bacteroidetes" = "gold1", "Chlamydiae" = "bisque3", "Chloroflexi" = "green", "Deinococcota" = "green", "Firmicutes" = "seagreen", "Gemmatimonadetes" = "darkturquoise", "Nitrospirae" = "blue", "OP3" = "darkblue", "Planctomycetes" = "deepskyblue","Proteobacteria" = "orchid2", "Thermi" = "hotpink", "Verrucomicrobia" = "purple1", "WS3" = "magenta1"))

sigtab_UG_0_3subp2

svg("sigtab_UG_0_3subp2.svg", width = 6, height = 4, pointsize=12)
plot(sigtab_UG_0_3subp2)
dev.off()


### DA Waialua Farm (WF)

sigtab_WF_0_2 <- sigtab_WF_0

sigtab_WF_0_2$abslog2FC <- abs(sigtab_WF_0_2$log2FoldChange)
sigtab_WF_0_2 <- subset(sigtab_WF_0_2, abslog2FC >= 0.945)

sigtab_WF_0_3 <- sigtab_WF_0_2

sigtab_WF_0_3sub <- subset(sigtab_WF_0_3, Family!="N/A" & Phylum != "N/A")
# #View(sigtab_WF_0_3sub)

sigtab_WF_0_3sub_Phyla <- unique(sigtab_WF_0_3sub$Phylum)
# #View(sigtab_WF_0_3sub_Phyla)


sigtab_WF_0_3subp = ggplot(sigtab_WF_0_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Waialua Farm Day 0: Selection 3 vs 5") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_WF_0_3subp

sigtab_WF_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_WF_0_3subp2 <- sigtab_WF_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank()) + theme(panel.spacing = unit(0.1, "lines")) + theme(text = element_text(size = 8)) + scale_colour_manual(values = c("Acidobacteria" = "red1", "Actinobacteria" = "tomato", "Armatimonadetes" = "darkorange", "Bacteroidetes" = "gold1", "Chlamydiae" = "bisque3", "Chloroflexi" = "green", "Deinococcota" = "green", "Firmicutes" = "seagreen", "Gemmatimonadetes" = "darkturquoise", "Nitrospirae" = "blue", "OP3" = "darkblue", "Planctomycetes" = "deepskyblue","Proteobacteria" = "orchid2", "Thermi" = "hotpink", "Verrucomicrobia" = "purple1", "WS3" = "magenta1"))

sigtab_WF_0_3subp2

svg("sigtab_WF_0_3subp2.svg", width = 6, height = 4, pointsize=12)
plot(sigtab_WF_0_3subp2)
dev.off()


##DESEQ_Summary Selections compared at each Site at Day 0

sigtab_UG_0_3subp2
sigtab_WF_0_3subp2


###DESeq Sites compared within each selection: examining the influence of site on each selection independentally 

# #View(data.frame(sample_data(cacao_16S_0)))

#DESeq Selection 3: UG vs WF Day 0

cacao_16S_S3 = subset_samples(cacao_16S_0, tree_selection == "3")
cacao_16S_S3_0_0 <- subset_samples(cacao_16S_S3, days_ferment == "0")
cacao_16S_S3_0 = filter_taxa(cacao_16S_S3_0_0, function(x) sum(x) > 0, TRUE)

# #View(data.frame(sample_data(cacao_16S_S3_0)))

cacao_16S_S3_0dds <- phyloseq_to_deseq2(cacao_16S_S3_0, ~Location)
cacao_16S_S3_0dds$tree_selection <- relevel( cacao_16S_S3_0dds$Location, "Urban Garden Center" )

cacao_16S_S3_0dds = DESeq(cacao_16S_S3_0dds, test="Wald", fitType="parametric")
rescacao_16S_S3_0dds = results(cacao_16S_S3_0dds, cooksCutoff = TRUE)
rescacao_16S_S3_0dds = lfcShrink(cacao_16S_S3_0dds,contrast = c("Location","Waialua Farm","Urban Garden Center"),res=rescacao_16S_S3_0dds)

alpha = 0.05
sigtab_S3_0 = rescacao_16S_S3_0dds[which(rescacao_16S_S3_0dds$padj <= alpha), ]
sigtab_S3_0 = cbind(as(sigtab_S3_0, "data.frame"), as(tax_table(cacao_16S_S3_0)[rownames(sigtab_S3_0), ], "matrix"))
head(sigtab_S3_0)
write.csv(sigtab_S3_0, "sig_S3_0.csv")

notsigtab_S3_0 = rescacao_16S_S3_0dds[which(rescacao_16S_S3_0dds$padj > alpha), ]
notsigtab_S3_0 = cbind(as(notsigtab_S3_0, "data.frame"), as(tax_table(cacao_16S_S3_0)[rownames(notsigtab_S3_0), ], "matrix"))
head(notsigtab_S3_0)
write.csv(notsigtab_S3_0, "notsig_S3_0.csv")

#DESeq Selection 5: UG vs WF Day 0 

cacao_16S_S5 = subset_samples(cacao_16S_0, tree_selection == "5")
cacao_16S_S5_0_0 <- subset_samples(cacao_16S_S5, days_ferment == "0")
cacao_16S_S5_0 = filter_taxa(cacao_16S_S5_0_0, function(x) sum(x) > 0, TRUE)

# #View(data.frame(sample_data(cacao_16S_S5_0)))

cacao_16S_S5_0dds <- phyloseq_to_deseq2(cacao_16S_S5_0, ~Location)
cacao_16S_S5_0dds$tree_selection <- relevel( cacao_16S_S5_0dds$Location, "Urban Garden Center" )

cacao_16S_S5_0dds = DESeq(cacao_16S_S5_0dds, test="Wald", fitType="parametric")
rescacao_16S_S5_0dds = results(cacao_16S_S5_0dds, cooksCutoff = TRUE)
rescacao_16S_S5_0dds = lfcShrink(cacao_16S_S5_0dds,contrast = c("Location","Waialua Farm","Urban Garden Center"),res=rescacao_16S_S5_0dds)

alpha = 0.05
sigtab_S5_0 = rescacao_16S_S5_0dds[which(rescacao_16S_S5_0dds$padj <= alpha), ]
sigtab_S5_0 = cbind(as(sigtab_S5_0, "data.frame"), as(tax_table(cacao_16S_S5_0)[rownames(sigtab_S5_0), ], "matrix"))
head(sigtab_S5_0)
# write.csv(sigtab_S5_0, "sig_S5_0.csv")

notsigtab_S5_0 = rescacao_16S_S5_0dds[which(rescacao_16S_S5_0dds$padj > alpha), ]
notsigtab_S5_0 = cbind(as(notsigtab_S5_0, "data.frame"), as(tax_table(cacao_16S_S5_0)[rownames(notsigtab_S5_0), ], "matrix"))
head(notsigtab_S5_0)
# write.csv(notsigtab_S5_0, "notsig_S5_0.csv")

### DA Selection 3: UG vs WF Day 0

sigtab_S3_0_2 <- sigtab_S3_0

sigtab_S3_0_2$abslog2FC <- abs(sigtab_S3_0_2$log2FoldChange)
sigtab_S3_0_2 <- subset(sigtab_S3_0_2, abslog2FC >= 0.945)
# #View(sigtab_S3_0_2)


sigtab_S3_0_3 <- sigtab_S3_0_2
sigtab_S3_0_3sub <- subset(sigtab_S3_0_3, Family!="N/A" & Phylum != "N/A")
# #View(sigtab_S3_0_3sub)

sigtab_S3_0_3sub_Phyla <- unique(sigtab_S3_0_3sub$Phylum)
# #View(sigtab_S3_0_3sub_Phyla)


sigtab_S3_0_3subp = ggplot(sigtab_S3_0_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Selection 3 Day 0: Urban Garden vs Waialua Farm") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_S3_0_3subp

sigtab_S3_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_S3_0_3subp2 <- sigtab_S3_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank()) + theme(panel.spacing = unit(0.1, "lines")) + theme(text = element_text(size = 8)) + scale_colour_manual(values = c("Acidobacteria" = "red1", "Actinobacteria" = "tomato", "Armatimonadetes" = "darkorange", "Bacteroidetes" = "gold1", "Chlamydiae" = "bisque3", "Chloroflexi" = "green", "Deinococcota" = "green", "Firmicutes" = "seagreen", "Gemmatimonadetes" = "darkturquoise", "Nitrospirae" = "blue", "OP3" = "darkblue", "Planctomycetes" = "deepskyblue","Proteobacteria" = "orchid2", "Thermi" = "hotpink", "Verrucomicrobia" = "purple1", "WS3" = "magenta1"))

sigtab_S3_0_3subp2

svg("sigtab_S3_0_3subp2.svg", width = 6, height = 6, pointsize=12)
plot(sigtab_S3_0_3subp2)
dev.off()


sigtab_S3_0_3sub2 <- subset(sigtab_S3_0_3, Phylum != "N/A")
# #View(sigtab_S3_0_3sub)

sigtab_S3_0_3sub2_Phyla <- unique(sigtab_S3_0_3sub2$Phylum)
# #View(sigtab_S3_0_3sub2_Phyla)


sigtab_S3_0_3sub2p = ggplot(sigtab_S3_0_3sub2, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Selection 3 Day 0: Urban Garden vs Waialua Farm") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_S3_0_3sub2p

sigtab_S3_0_3sub2p + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_S3_0_3sub2p2 <- sigtab_S3_0_3sub2p + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank()) + theme(panel.spacing = unit(0.1, "lines")) + theme(text = element_text(size = 8)) + scale_colour_manual(values = c("Acidobacteria" = "red1", "Actinobacteria" = "tomato", "Armatimonadetes" = "darkorange", "Bacteroidetes" = "gold1", "Chlamydiae" = "bisque3", "Chloroflexi" = "green", "Deinococcota" = "green", "Firmicutes" = "seagreen", "Gemmatimonadetes" = "darkturquoise", "Nitrospirae" = "blue", "OP3" = "darkblue", "Planctomycetes" = "deepskyblue","Proteobacteria" = "orchid2", "Thermi" = "hotpink", "Verrucomicrobia" = "purple1", "WS3" = "magenta1"))

sigtab_S3_0_3sub2p2

svg("sigtab_S3_0_3sub2p2.svg", width = 6, height = 6, pointsize=12)
plot(sigtab_S3_0_3sub2p2)
dev.off()





### DA Selection 5: UG vs WF Day 0

sigtab_S5_0_2 <- sigtab_S5_0

sigtab_S5_0_2$abslog2FC <- abs(sigtab_S5_0_2$log2FoldChange)
sigtab_S5_0_2 <- subset(sigtab_S5_0_2, abslog2FC >= 0.945)



sigtab_S5_0_3 <- sigtab_S5_0_2

sigtab_S5_0_3sub <- subset(sigtab_S5_0_3, Family!="N/A" & Phylum != "N/A")
# #View(sigtab_S5_0_3sub)

sigtab_S5_0_3sub_Phyla <- unique(sigtab_S5_0_3sub$Phylum)
# #View(sigtab_S5_0_3sub_Phyla)


sigtab_S5_0_3subp = ggplot(sigtab_S5_0_3sub, aes(x=log2FoldChange, y=reorder(Fam_Gen,desc(Fam_Gen)), color=Phylum))  + geom_point(size=2, stroke =0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.25)) + theme(legend.position="none") + ggtitle("DA OTUs Selection 5 Day 0: Urban Garden vs Waialua Farm") + geom_vline(xintercept = 0, linetype = "solid", color = "black")

sigtab_S5_0_3subp

sigtab_S5_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank())

sigtab_S5_0_3subp2 <- sigtab_S5_0_3subp + facet_grid(Phylum~., scales= "free_y", space="free_y", switch = "y") + theme(strip.text.y.left = element_text(angle = 0)) + scale_y_discrete(position = "right") + theme(axis.title.y=element_blank()) + theme(panel.spacing = unit(0.1, "lines")) + theme(text = element_text(size = 8)) + scale_colour_manual(values = c("Acidobacteria" = "red1", "Actinobacteria" = "tomato", "Armatimonadetes" = "darkorange", "Bacteroidetes" = "gold1", "Chlamydiae" = "bisque3", "Chloroflexi" = "green", "Deinococcota" = "green", "Firmicutes" = "seagreen", "Gemmatimonadetes" = "darkturquoise", "Nitrospirae" = "blue", "OP3" = "darkblue", "Planctomycetes" = "deepskyblue","Proteobacteria" = "orchid2", "Thermi" = "hotpink", "Verrucomicrobia" = "purple1", "WS3" = "magenta1"))

sigtab_S5_0_3subp2

svg("sigtab_S5_0_3subp2.svg", width = 6, height = 6, pointsize=12)
plot(sigtab_S5_0_3subp2)
dev.off()

###Day 0 Selection Summary

sigtab_S3_0_3subp
sigtab_S5_0_3subp

##Day 0 total summary

sigtab_UG_0_3subp2
sigtab_WF_0_3subp2
sigtab_S3_0_3subp2
sigtab_S5_0_3subp2


par(mfrow=c(2,1))
sigtab_UG_0_3subp2
sigtab_WF_0_3subp2

Site_DA <- ggarrange(sigtab_UG_0_3subp2, sigtab_WF_0_3subp2,
                     labels = c("A", "B"),
                     heights = c(7,7),
                     ncol = 1, nrow = 2,
                     align = "v")
Site_DA

