# ---------------------------------------------------------
# OUwie models with log and sd data - BIO17
# ---------------------------------------------------------

# This is part of:

#Resource Availability and Disturbance Frequency Shape Plant Life Forms in Neotropical Habitats 

#Authors: Luana S. Prochazka, Suzana Alcantara, Juliana Gastaldello Rando, Thais Vasconcelos, Raquel C. Pizzardo and Anselmo Nogueira

#Journal: New Phytologist
#Article acceptance date: 30 January 2024


#uncomment to clean your R environment:

#rm(list = ls())

##packages:

library(dplyr)
library(data.table)
library(phyloclim)
library(plotrix)
library(car)
library(geiger)
library(picante)
library(phytools)
library(OUwie)
library(nloptr)
library(tibble)
library(data.table)
library(qpcR)
library(rlist)


## bio17 pno and life form classification: 

bio17 <- read.csv("/data-paper/Results/01-pno/pno_bio17.csv")

life.form  <- read.csv("/data-paper/data/life_form.csv")

names <- colnames(bio17)

setdiff(names, life.form$specie)

##Remove species that are not in the phylogeny:

bio17 <-subset(bio17, select = -c(Chamaecrista_debilis, Chamaecrista_kunthiana, Chamaecrista_serpens, Chamaecrista_basifolia,Chamaecrista_desvauxii_var_brevipes, Chamaecrista_ramosa_var_erythrocalyx))


##PNO sample: 

ntax <- 41 ##species number
n <- 500   ## sampling amount

x <- matrix(nrow = n, ncol = ntax)

colnames(x) <- colnames(bio17)[-1]

##Make 500 samples:

for (i in 2:dim(bio17)[2]) { 
  x[, i - 1] <- sample(bio17[, 1], size = n, 
                       replace = TRUE, prob = bio17[, i])}

write.csv(x, "data-paper/Results/04-OUwie_and_WN_models/bio17/bio17_500samples_pno.csv", row.names = F)

 
##Log transformation, mean and sd calc: 

bio17.log <- log(x + 1)

bio17.mean <- apply(bio17.log, 2, mean)
bio17.sd <-  apply(bio17.log, 2, std.error)
data.bio17 <- cbind(bio17.mean, bio17.sd)
data.bio17 <- as.data.frame(data.bio17)
data.bio17 <- tibble::rownames_to_column(data.bio17, "specie")

##Join data:

data <- left_join(life.form, data.bio17)

write.csv(data, "data-paper/Results/04-OUwie_and_WN_models/bio17/bio17log_data.csv", row.names = F)

###SIMMAP data:

SIMMAP <- readRDS("data-paper/Results/03-Life_form_ancestral_reconstruction-SIMMAP/simmap_ER_model.rds")
rownames(data) <-data$specie


##set working directory to the data-paper folder:

setwd("/data-paper/Results/04-OUwie_and_WN_models/bio17/")

###Model fit:

fitBM1 <- list()
for(i in 1:length(SIMMAP)) {
  fitBM1[[i]] <-OUwie(SIMMAP[[i]], data, model="BM1", simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitBM1, "bio17_fitBM1.rds")

fitBMS <- list()
for(i in 1:length(SIMMAP)) {
  fitBMS[[i]] <-OUwie(SIMMAP[[i]], data, model="BMS",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitBMS,  "bio17_fitBMS.rds")

fitOU1 <- list()
for(i in 1:length(SIMMAP)) {
  fitOU1[[i]] <-OUwie(SIMMAP[[i]], data, model="OU1",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOU1, "bio17_fitOU1.rds")

fitOUM <- list()
for(i in 1:length(SIMMAP)) {
  fitOUM[[i]]<-OUwie(SIMMAP[[i]], data, model="OUM",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUM, "bio17_fitOUM.rds")

fitOUMV <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMV[[i]]<-OUwie(SIMMAP[[i]], data, model="OUMV",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMV, "bio17_fitOUMV.rds")

fitOUMA <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMA[[i]] <-OUwie(SIMMAP[[i]], data, model="OUMA",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMA, "bio17_fitOUMA.rds")

fitOUMVA <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMVA[[i]]<- OUwie(SIMMAP[[i]], data,model="OUMVA",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMVA, "bio17_fitOUMVA.rds")


# ---------------------------------------------------------
# OUwie models with log and sd data - BIO18
# ---------------------------------------------------------

rm(list = ls())

# bio18 pno and life form classification: 

bio18 <- read.csv("/data-paper/Results/01-pno/pno_bio18.asc.csv")

life.form  <- read.csv("/data-paper/data/life_form.csv")

names <- colnames(bio18)

setdiff(names, life.form$specie)

##Remove species that are not in the phylogeny:

bio18 <-subset(bio18, select = -c(Chamaecrista_debilis, Chamaecrista_kunthiana, Chamaecrista_serpens, Chamaecrista_basifolia,Chamaecrista_desvauxii_var_brevipes, Chamaecrista_ramosa_var_erythrocalyx))


##PNO sample: 

ntax <- 41 ##species number
n <- 500   ## sampling amount


x <- matrix(nrow = n, ncol = ntax)

colnames(x) <- colnames(bio18)[-1]

##Make 500 samples:

for (i in 2:dim(bio18)[2]) { 
  x[, i - 1] <- sample(bio18[, 1], size = n, 
                       replace = TRUE, prob = bio18[, i])}

write.csv(x, "data-paper/Results/04-OUwie_and_WN_models/bio18/bio18_500samples_pno.csv", row.names = F)

##Log transformation, mean and sd calc:  

bio18.log <- log(x + 1)

bio18.mean <- apply(bio18.log, 2, mean)
bio18.sd <-  apply(bio18.log, 2, std.error)
data.bio18 <- cbind(bio18.mean, bio18.sd)
data.bio18 <- as.data.frame(data.bio18)
data.bio18 <- tibble::rownames_to_column(data.bio18, "specie")

##Join data:

data <- left_join(life.form, data.bio18)

write.csv(data, "data-paper/Results/04-OUwie_and_WN_models/bio18/bio18log_data.csv", row.names = F)

###SIMMAP data:

SIMMAP <- readRDS("simmap_ER_model.rds")
rownames(data) <-data$specie

##set working directory to the data-paper folder:

setwd("/data-paper/Results/04-OUwie_and_WN_models/bio18/")

###Model fit:

fitBM1 <- list()
for(i in 1:length(SIMMAP)) {
  fitBM1[[i]] <-OUwie(SIMMAP[[i]], data, model="BM1", simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitBM1, "bio18_fitBM1.rds")

fitBMS <- list()
for(i in 1:length(SIMMAP)) {
  fitBMS[[i]] <-OUwie(SIMMAP[[i]], data, model="BMS",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitBMS,  "bio18_fitBMS.rds")

fitOU1 <- list()
for(i in 1:length(SIMMAP)) {
  fitOU1[[i]] <-OUwie(SIMMAP[[i]], data, model="OU1",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOU1, "bio18_fitOU1.rds")

fitOUM <- list()
for(i in 1:length(SIMMAP)) {
  fitOUM[[i]]<-OUwie(SIMMAP[[i]], data, model="OUM",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUM, "bio18_fitOUM.rds")

fitOUMV <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMV[[i]]<-OUwie(SIMMAP[[i]], data, model="OUMV",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMV, "bio18_fitOUMV.rds")

fitOUMA <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMA[[i]] <-OUwie(SIMMAP[[i]], data, model="OUMA",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMA, "bio18_fitOUMA.rds")

fitOUMVA <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMVA[[i]]<- OUwie(SIMMAP[[i]], data,model="OUMVA",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMVA, "bio18_fitOUMVA.rds")

# ---------------------------------------------------------
# OUwie models with log and sd data - FIRE
# ---------------------------------------------------------

rm(list = ls())

## Fire pno and life form classification: 

fire <- read.csv("/data-paper/Results/01-pno/pno_fire_frequency_0.01.csv")

life.form  <- read.csv("/data-paper/data/life_form.csv")

names <- colnames(fire)

setdiff(names, life.form$specie)

##Remove species that are not in the phylogeny: 

fire <-subset(fire, select = -c(Chamaecrista_debilis, Chamaecrista_kunthiana, Chamaecrista_serpens, Chamaecrista_basifolia,Chamaecrista_desvauxii_var_brevipes, Chamaecrista_ramosa_var_erythrocalyx))


##PNO sample: 

ntax <- 41 ##species number
n <- 500   ## sampling amount

x <- matrix(nrow = n, ncol = ntax)

colnames(x) <- colnames(fire)[-1]

##Make 500 samples:

for (i in 2:dim(fire)[2]) { 
  x[, i - 1] <- sample(fire[, 1], size = n, 
                       replace = TRUE, prob = fire[, i])}

write.csv(x, "data-paper/Results/04-OUwie_and_WN_models/fire/fire_500samples_pno.csv", row.names = F)

##Log transformation, mean and sd calc: 

fire.log <- log(x + 1)

fire.mean <- apply(fire.log, 2, mean)
fire.sd <-  apply(fire.log, 2, std.error)
data.fire <- cbind(fire.mean, fire.sd)
data.fire <- as.data.frame(data.fire)
data.fire <- tibble::rownames_to_column(data.fire, "specie")


##Join data:

data <- left_join(life.form, data.fire)

write.csv(data, "data-paper/Results/04-OUwie_and_WN_models/fire/firelog_data.csv", row.names = F)

###SIMMAP data:

SIMMAP <- readRDS("simmap_ER_model.rds")
rownames(data) <-data$specie

##set working directory to the data-paper folder:

setwd("/data-paper/Results/04-OUwie_and_WN_models/fire/")

###Model fit:

fitBM1 <- list()
for(i in 1:length(SIMMAP)) {
  fitBM1[[i]] <-OUwie(SIMMAP[[i]], data, model="BM1", simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitBM1, "fire_fitBM1.rds")

fitBMS <- list()
for(i in 1:length(SIMMAP)) {
  fitBMS[[i]] <-OUwie(SIMMAP[[i]], data, model="BMS",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitBMS,  "fire_fitBMS.rds")

fitOU1 <- list()
for(i in 1:length(SIMMAP)) {
  fitOU1[[i]] <-OUwie(SIMMAP[[i]], data, model="OU1",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOU1, "fire_fitOU1.rds")

fitOUM <- list()
for(i in 1:length(SIMMAP)) {
  fitOUM[[i]]<-OUwie(SIMMAP[[i]], data, model="OUM",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUM, "fire_fitOUM.rds")

fitOUMV <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMV[[i]]<-OUwie(SIMMAP[[i]], data, model="OUMV",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMV, "fire_fitOUMV.rds")

fitOUMA <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMA[[i]] <-OUwie(SIMMAP[[i]], data, model="OUMA",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMA, "fire_fitOUMA.rds")

fitOUMVA <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMVA[[i]]<- OUwie(SIMMAP[[i]], data,model="OUMVA",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMVA, "fire_fitOUMVA.rds")

# ---------------------------------------------------------
# OUwie models with log and sd data - Nitrogen
# ---------------------------------------------------------

rm(list = ls())

## Nitrogen pno and life form classification: 

nitro <- read.csv("/data-paper/Results/01-pno/pno_nitrogen_0_30cm.asc.csv")

life.form  <- read.csv("/data-paper/data/life_form.csv")

names <- colnames(nitro)

setdiff(names, life.form$specie)

##Remove species that are not in the phylogeny: 

nitro <-subset(nitro, select = -c(Chamaecrista_debilis, Chamaecrista_kunthiana, Chamaecrista_serpens, Chamaecrista_basifolia,Chamaecrista_desvauxii_var_brevipes, Chamaecrista_ramosa_var_erythrocalyx))


##PNO sample: 

ntax <- 41 ##species number
n <- 500   ## sampling amount

x <- matrix(nrow = n, ncol = ntax)

colnames(x) <- colnames(nitro)[-1]

##Make 500 samples:

for (i in 2:dim(nitro)[2]) { 
  x[, i - 1] <- sample(nitro[, 1], size = n, 
                       replace = TRUE, prob = nitro[, i])}

write.csv(x, "data-paper/Results/04-OUwie_and_WN_models/nitro/nitro_500samples_pno.csv", row.names = F)

##Log transformation, mean and sd calc: 

nitro.log <- log(x + 1)

nitro.mean <- apply(nitro.log, 2, mean)
nitro.sd <-  apply(nitro.log, 2, std.error)
data.nitro <- cbind(nitro.mean, nitro.sd)
data.nitro <- as.data.frame(data.nitro)
data.nitro <- tibble::rownames_to_column(data.nitro, "specie")

##Join data:

data <- left_join(life.form, data.nitro)

write.csv(data, "data-paper/Results/04-OUwie_and_WN_models/nitro/nitrolog_data.csv", row.names = F)

###SIMMAP data:

SIMMAP <- readRDS("data-paper/Results/03-Life_form_ancestral_reconstruction-SIMMAP/simmap_ER_model.rds")
rownames(data) <-data$specie


##set working directory to the data-paper folder:

setwd("/data-paper/Results/04-OUwie_and_WN_models/nitro/")

###Model fit:

fitBM1 <- list()
for(i in 1:length(SIMMAP)) {
  fitBM1[[i]] <-OUwie(SIMMAP[[i]], data, model="BM1", simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitBM1, "nitro_fitBM1.rds")

fitBMS <- list()
for(i in 1:length(SIMMAP)) {
  fitBMS[[i]] <-OUwie(SIMMAP[[i]], data, model="BMS",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitBMS,  "nitro_fitBMS.rds")

fitOU1 <- list()
for(i in 1:length(SIMMAP)) {
  fitOU1[[i]] <-OUwie(SIMMAP[[i]], data, model="OU1",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOU1, "nitro_fitOU1.rds")

fitOUM <- list()
for(i in 1:length(SIMMAP)) {
  fitOUM[[i]]<-OUwie(SIMMAP[[i]], data, model="OUM",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUM, "nitro_fitOUM.rds")

fitOUMV <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMV[[i]]<-OUwie(SIMMAP[[i]], data, model="OUMV",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMV, "nitro_fitOUMV.rds")

fitOUMA <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMA[[i]] <-OUwie(SIMMAP[[i]], data, model="OUMA",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMA, "nitro_fitOUMA.rds")

fitOUMVA <- list()
for(i in 1:length(SIMMAP)) {
  fitOUMVA[[i]]<- OUwie(SIMMAP[[i]], data,model="OUMVA",simmap.tree=TRUE, diagn = TRUE, mserr = "known", algorithm= "invert")}
saveRDS(fitOUMVA, "nitro_fitOUMVA.rds")

##end