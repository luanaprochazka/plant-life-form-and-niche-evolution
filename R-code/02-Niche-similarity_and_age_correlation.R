# -----------------------------------------------------------
# D and I metrics for niche similarity between species and
# ARC correlation
# (See Warren, 2008)
# -----------------------------------------------------------

# This is part of:

#Resource Availability and Disturbance Frequency Shape Plant Life Forms in Neotropical Habitats 

#Authors: Luana S. Prochazka, Suzana Alcantara, Juliana Gastaldello Rando, Thais Vasconcelos, Raquel C. Pizzardo and Anselmo Nogueira

#Journal: New Phytologist
#Article acceptance date: 30 January 2024


# uncomment to clean your R environment:

#rm(list = ls())

##packages: 

library(dplyr)
library(tidyverse)
library(phyloclim)
library(car)
library(geiger)
library(picante)
library(phytools)


## set working directory to the data folder in the data-paper:

phylo <- read.tree("/data-paper/data/tree_vasconcelos41spp.txt")

##Phylogeny:

plot(phylo)


#remove species that are not on the phylogeny:

spp <- c("Chamaecrista_basifolia", "Chamaecrista_debilis", "Chamaecrista_desvauxii_var_brevipes", "Chamaecrista_kunthiana", "Chamaecrista_ramosa_var_erythrocalyx", "Chamaecrista_serpens")

## set working directory with the maxent models results (.asc files): 

setwd("/maxent_models_results/")

lista <-list.files(pattern = ".asc")[-c(4, 13, 16, 24, 33, 40)]

###calc niche overlap: 

geo_overlap <- niche.overlap(lista)

setwd("/data-paper/Results/02-Niche_overlap_and_age_correlation/")

write.csv(geo_overlap, "geo_overlap.csv")


##Calc niche_verlap to each environmental layer: 

##PNO list: 

setwd("/data-paper/Results/01-pno")
pno.lista <- list.files(pattern = ".csv")

for(i in pno.lista) {

setwd("/data-paper/Results/01-pno")
pno <- read.csv(i)

#remove species that are not on the phylogeny:

pno2 <- pno[ , !(names(pno) %in% spp)]

##Calc Niche overlap: The upper triangle contains pairwise comparisons of niche overlap in terms of D, whereas the lower triangle contains values of I):

overlap <- niche.overlap(pno2)
 
setwd("/data-paper/Results/02-Niche_overlap_and_age_correlation")

write.csv(overlap, paste0("overlap_", i))

}

##############################################################
#Niche similarity and node age - Age range correlation - ARC 
#(See Warren, 2008)
##############################################################

#clean your R environment:

rm(list = ls())

#Packages:

library(dplyr)
library(tidyverse)
library(phyloclim)
library(car)
library(geiger)
library(picante)
library(phytools)

##Phylogeny:

phylo <- read.tree("/data-paper/data/tree_vasconcelos41spp.txt")


##Niche overlap files:

setwd("/data-paper/Results/02-Niche_overlap_and_age_correlation")
list <- list.files(pattern = ".csv")[-1]

for (i in list) {
  
  setwd("/data-paper/Results/02-Niche_overlap_and_age_correlation")
  
  overlap <- read.csv(i, header = T, row.names = 1)
  
  
  name <- gsub(".asc.csv", '', i)
  name <- gsub("overlap_pno_", '', name)
  
  overlap <-as.matrix(overlap)
  
  setwd("/data-paper/Results/02-Niche_overlap_and_age_correlation")
  
  ##Calc ARC: 
  
  corr.uper <- age.range.correlation(phylo, overlap, tri = "upper", n = 1000)
  
  saveRDS(corr.uper, file = paste0("age_correlation_upper_", name, ".rds"))
  
  corr.lower <- age.range.correlation(phylo, overlap, tri = "lower", n = 1000)
  
  saveRDS(corr.lower, file = paste0("age_correlation_lower_", name, ".rds"))
  
  
}
## Make ARC analyses to file geo_overlap: 


setwd("/data-paper/Results/02-Niche_overlap_and_age_correlation")

overlap <- read.csv("geo_overlap.csv",header = T, row.names = 1)


colnames(overlap) <- gsub(".asc", '', colnames(overlap))
colnames(overlap)

rownames(overlap) <- gsub(".asc", '', rownames(overlap))
rownames(overlap)


setdiff(phylo$tip.label, rownames(overlap))

rownames(overlap)[10] <- "Chamaecrista_cinerascens"
rownames(overlap)[26] <- "Chamaecrista_pascuorum"

setdiff(phylo$tip.label, colnames(overlap))

colnames(overlap)[10] <- "Chamaecrista_cinerascens"
colnames(overlap)[26] <- "Chamaecrista_pascuorum"


overlap <-as.matrix(overlap)

setwd("/data-paper/Results/02-Niche_overlap_and_age_correlation")

corr.uper <- age.range.correlation(phylo, overlap, tri = "upper", n = 1000)

saveRDS(corr.uper, file = "age_correlation_upper_geo.rds")

corr.uper$age.range.correlation

corr.lower <- age.range.correlation(phylo, overlap, tri = "lower", n = 1000)

saveRDS(corr.lower, file = "age_correlation_lower_geo.rds")

##end 