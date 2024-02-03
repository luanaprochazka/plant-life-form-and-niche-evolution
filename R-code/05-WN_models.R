# -------------------------------------------
# White noise models and AICc values
# -------------------------------------------

# This is part of:

#Resource Availability and Disturbance Frequency Shape Plant Life Forms in Neotropical Habitats 

#Authors: Luana S. Prochazka, Suzana Alcantara, Juliana Gastaldello Rando, Thais Vasconcelos, Raquel C. Pizzardo and Anselmo Nogueira

#Journal: New Phytologist
#Article acceptance date: 30 January 2024


#uncomment to clean your R environment:

#rm(list = ls())

##packages:  

library(car)
library(geiger)
library(picante)
library(phytools)
library(dplyr)
library(OUwie)
library(nloptr)
library(tibble)
library(data.table)
library(qpcR)
library(dplyr)
library(rlist)
library (ape)
library(OUwie)
library(ggplot2)
library(cowplot)
library(picante)
library(rlist)
library(qpcR)
library(tibble)
library(dplyr)
library(tidyr)

##set working directory to the data-paper data:

##Data:

##SIMMAP

SIMMAP<- readRDS("/data-paper/Results/03-Life_form_ancestral_reconstruction-SIMMAP/simmap_ER_model.rds")

bio17 <- read.csv("/data-paper/Results/04-OUwie_and_WN_models/bio17/bio17log_data.csv")

bio18 <- read.csv("/data-paper/Results/04-OUwie_and_WN_models/bio18/bio18log_data.csv")

fire <- read.csv("/data-paper/Results/04-OUwie_and_WN_models/fire/firelog_data.csv")

nitro <- read.csv("data-paper/Results/04-OUwie_and_WN_models/nitro/nitrolog_data.csv")


##Fit white noise in bio17:

bio17.WN <- bio17[,3]
names(bio17.WN)<- bio17$specie

bio17.SE <- bio17[,4]
names(bio17.SE)<- bio17$specie


fitWN <- list()
for(i in 1:length(SIMMAP)) {
  
  fitWN[[i]]<- fitContinuous(SIMMAP[[i]], bio17.WN, SE = bio17.SE, model="white")
  
}

saveRDS(fitWN, "/data-paper/Results/04-OUwie_and_WN_models/bio17/bio17_fitWN.rds")


##Fit white noise in bio18:

bio18.WN <- bio18[,3]
names(bio18.WN)<- bio18$specie

bio18.SE <- bio18[,4]
names(bio18.SE)<- bio18$specie

fitWN <- list()
for(i in 1:length(SIMMAP)) {
  
  fitWN[[i]]<- fitContinuous(SIMMAP[[i]], bio18.WN,SE = bio18.SE, model="white")
  
}

saveRDS(fitWN, "/data-paper/Results/04-OUwie_and_WN_models/bio18/bio18_fitWN.rds")

##Fit white noise in fire:

fire.WN <- fire[,3]
names(fire.WN)<- fire$specie

fire.SE <- fire[,4]
names(fire.SE)<- fire$specie

fitWN <- list()
for(i in 1:length(SIMMAP)) {
  
  fitWN[[i]]<- fitContinuous(SIMMAP[[i]],fire.WN, SE = fire.SE, model="white")
  
}

saveRDS(fitWN, "/data-paper/Results/04-OUwie_and_WN_models/fire/fire_fitWN.rds")

##Fit white noise in nitro:

nitro.WN <- nitro[,3]
names(nitro.WN)<- nitro$specie

nitro.SE <- nitro[,4]
names(nitro.SE)<- nitro$specie


fitWN <- list()
for(i in 1:length(SIMMAP)) {
  
  fitWN[[i]]<- fitContinuous(SIMMAP[[i]], nitro.WN,SE = nitro.SE, model="white")
  
}

saveRDS(fitWN, "/data-paper/Results/04-OUwie_and_WN_models/nitro/nitro_fitWN.rds")

##end