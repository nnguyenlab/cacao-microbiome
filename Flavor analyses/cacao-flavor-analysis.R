#set working directory
setwd("/Volumes/GoogleDrive/.shortcut-targets-by-id/1--pvx1vLSBijmxHMCZPrxKlyFBLp2jqF/Rick/cacao_combined/cacao_physeq_full_analysis/cacao_Fungal_ITS_ps")

#Load libraries
library(tidyverse)
library(fmsb)

#Read in data
flavor_panel <- read.csv("cacao_flavor_panel.csv") %>% 
  as_tibble ()

#Filter to just year
flavor2016 <- filter(flavor_panel, year == "2016")
flavor2015 <- filter(flavor_panel, year == "2015")

#----2016 dataset TESTING FOR NORMALITY ----------
# Using the Shapiro Test -- if p-value not significant = data is normal
shapiro.test(flavor2016$Cocoa)#NS
shapiro.test(flavor2016$Acidity)#NS
shapiro.test(flavor2016$Bitterness)#S
shapiro.test(flavor2016$Astringency)#S
shapiro.test(flavor2016$Sweet)#S
shapiro.test(flavor2016$Fresh_Fruit)#S
shapiro.test(flavor2016$Browned_Fruit)#NS
shapiro.test(flavor2016$Nutty)#S
shapiro.test(flavor2016$Floral)#S
shapiro.test(flavor2016$Woody)#S
shapiro.test(flavor2016$Spicy)#S
shapiro.test(flavor2016$Off_Flavors)#S
#Most are significant so cannot use ANOVA


##--------Using the Kruskal Wallis test--------
kruskal.test(Cocoa~Location, data=flavor2016)#NS
kruskal.test(Acidity~Location, data=flavor2016)#NS
kruskal.test(Bitterness~Location, data=flavor2016)#NS
kruskal.test(Astringency~Location, data=flavor2016)#NS
kruskal.test(Sweet~Location, data=flavor2016)#NS
kruskal.test(Fresh_Fruit~Location, data=flavor2016)#NS
kruskal.test(Browned_Fruit~Location, data=flavor2016)#NS
kruskal.test(Nutty~Location, data=flavor2016)#NS
kruskal.test(Floral~Location, data=flavor2016)#SIGNIFICANT
kruskal.test(Woody~Location, data=flavor2016)#NS
kruskal.test(Spicy~Location, data=flavor2016)#NS
kruskal.test(Off_Flavors~Location, data=flavor2016)#NS
#Only floral flavors were significant across the two sites

#--------testing for location x site interactions ----------
kruskal.test(Cocoa~LocationVar, data=flavor2016)#S
kruskal.test(Acidity~LocationVar, data=flavor2016)#NS
kruskal.test(Bitterness~LocationVar, data=flavor2016)#NS
kruskal.test(Astringency~LocationVar, data=flavor2016)#NS
kruskal.test(Sweet~LocationVar, data=flavor2016)#NS
kruskal.test(Fresh_Fruit~LocationVar, data=flavor2016)#NS
kruskal.test(Browned_Fruit~LocationVar, data=flavor2016)#NS
kruskal.test(Nutty~LocationVar, data=flavor2016)#NS
kruskal.test(Floral~LocationVar, data=flavor2016)#NS
kruskal.test(Woody~LocationVar, data=flavor2016)#NS
kruskal.test(Spicy~LocationVar, data=flavor2016)#NS
kruskal.test(Off_Flavors~LocationVar, data=flavor2016)#S

#pairwise test
pairwise.wilcox.test(flavor2016$Cocoa, flavor2016$LocationVar, p.adjust.method = "BH")#can't compute
pairwise.wilcox.test(flavor2016$Off_Flavors, flavor2016$LocationVar, p.adjust.method = "BH")#can't compute

#----2015 dataset - TESTING FOR NORMALITY ----------
# Using the Shapiro Test -- if p-value not significant = data is normal
shapiro.test(flavor2015$Cocoa)#S
shapiro.test(flavor2015$Acidity)#S
shapiro.test(flavor2015$Bitterness)#S
shapiro.test(flavor2015$Astringency)#S
shapiro.test(flavor2015$Sweet)#S
shapiro.test(flavor2015$Fresh_Fruit)#S
shapiro.test(flavor2015$Browned_Fruit)#S
shapiro.test(flavor2015$Nutty)#S
shapiro.test(flavor2015$Floral)#S
shapiro.test(flavor2015$Woody)#S
shapiro.test(flavor2015$Spicy)#S
shapiro.test(flavor2015$Off_Flavors)#S
#All are significant so cannot use ANOVA


##--------Using the Kruskal Wallis test--------
kruskal.test(Cocoa~Location, data=flavor2015)#NS
kruskal.test(Acidity~Location, data=flavor2015)#NS
kruskal.test(Bitterness~Location, data=flavor2015)#NS
kruskal.test(Astringency~Location, data=flavor2015)#NS
kruskal.test(Sweet~Location, data=flavor2015)#NS
kruskal.test(Fresh_Fruit~Location, data=flavor2015)#NS
kruskal.test(Browned_Fruit~Location, data=flavor2015)#NS
kruskal.test(Nutty~Location, data=flavor2015)#NS
kruskal.test(Floral~Location, data=flavor2015)#NS
kruskal.test(Woody~Location, data=flavor2015)#NS
kruskal.test(Spicy~Location, data=flavor2015)#NS
kruskal.test(Off_Flavors~Location, data=flavor2015)#NS
#No significant flavors across the two sites

#--------testing for location x site interactions ----------
kruskal.test(Cocoa~LocationVar, data=flavor2015)#NS
kruskal.test(Acidity~LocationVar, data=flavor2015)#NS
kruskal.test(Bitterness~LocationVar, data=flavor2015)#NS
kruskal.test(Astringency~LocationVar, data=flavor2015)#S
kruskal.test(Sweet~LocationVar, data=flavor2015)#NS
kruskal.test(Fresh_Fruit~LocationVar, data=flavor2015)#NS
kruskal.test(Browned_Fruit~LocationVar, data=flavor2015)#NS
kruskal.test(Nutty~LocationVar, data=flavor2015)#NS
kruskal.test(Floral~LocationVar, data=flavor2015)#NS
kruskal.test(Woody~LocationVar, data=flavor2015)#NS
kruskal.test(Spicy~LocationVar, data=flavor2015)#NS
kruskal.test(Off_Flavors~LocationVar, data=flavor2015)#S
#Only Astringency is significant

#--------Draw spider chart ----------
#Write to csv and make manual edits to combine categories
write.csv(flavor2016, "flavor_spider.csv")

#Read back in
flavor_spider <- read.csv("flavor_spider.csv", header=TRUE) %>% 
  as_tibble ()

#Subset data
flavor_spider_Kualoa <- filter(flavor_spider, Location == "Kualoa Ranch") %>% 
  select(-Location) #remove Location column

flavor_spider_Urban <- filter(flavor_spider, Location == "Urban Garden Center") %>% 
  select(-Location) #remove Location column

#Draw nets
radarchart(flavor_spider_Kualoa, axistype=1, 
           #custom polygon 
           pcol=rgb(0.27, 0.68, 0.63, 0.7), pfcol=rgb(0.27, 0.68, 0.63,0.3), plwd=2, 
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,4.5,1), cglwd=0.8,
            #custom labels
           vlcex=0.8)

flavor_spider_Urban_plot <-radarchart(flavor_spider_Urban, axistype=1, 
           #custom polygon 
           pcol=rgb(0.99, 0.55, 0.24, 0.7), pfcol=rgb(0.99, 0.55, 0.24,0.3), plwd=2, 
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,4.5,1), cglwd=0.8,
           #custom labels
           vlcex=0.8)
