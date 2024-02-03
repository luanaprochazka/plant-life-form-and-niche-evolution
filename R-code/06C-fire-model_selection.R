# -----------------------------------------------------------
# AICc values, model plausibility and
# Delta AICc calculation 
# -----------------------------------------------------------

# -----------------------------------------------------------
# AICc values and model plausibility - OUwie models - FIRE
# -----------------------------------------------------------

# This is part of:

#Resource Availability and Disturbance Frequency Shape Plant Life Forms in Neotropical Habitats 

#Authors: Luana S. Prochazka, Suzana Alcantara, Juliana Gastaldello Rando, Thais Vasconcelos, Raquel C. Pizzardo and Anselmo Nogueira

#Journal: New Phytologist
#Article acceptance date: 30 January 2024

library(dplyr)
library(tidyverse)

rm(list = ls())

path <- "C:/data-paper/"

setwd(paste0(path, "Results/04-OUwie_and_WN_models/"))

##Species data (in log): 

fire <- read.csv("fire/firelog_data.csv")

##Model Results: 

BM1   <- readRDS("fire/fire_fitBM1.rds")
BMS   <- readRDS("fire/fire_fitBMS.rds")
OU1   <- readRDS("fire/fire_fitOU1.rds")
OUM   <- readRDS("fire/fire_fitOUM.rds")
OUMA  <- readRDS("fire/fire_fitOUMA.rds")
OUMV  <- readRDS("fire/fire_fitOUMV.rds")
OUMVA <- readRDS("fire/fire_fitOUMVA.rds")

##model number: 

SIMMAP <- c(1:2000)


##Extract AICc value:

BM1_AICc <- as.vector(unlist(lapply(BM1, function (x) x[c('AICc')])))
BMS_AICc  <- as.vector(unlist(lapply(BMS, function (x) x[c('AICc')])))
OU1_AICc  <- as.vector(unlist(lapply(OU1, function (x) x[c('AICc')])))
OUM_AICc  <- as.vector(unlist(lapply(OUM, function (x) x[c('AICc')])))
OUMA_AICc  <- as.vector(unlist(lapply(OUMA, function (x) x[c('AICc')])))
OUMV_AICc  <- as.vector(unlist(lapply(OUMV, function (x) x[c('AICc')])))
OUMVA_AICc  <- as.vector(unlist(lapply(OUMVA, function (x) x[c('AICc')])))

fire_AICc <- cbind(SIMMAP, BM1_AICc, BMS_AICc, OU1_AICc, OUM_AICc, OUMA_AICc, OUMV_AICc, OUMVA_AICc)
fire_AICc <- as.data.frame(fire_AICc)

##Check the quality adjust for each model: 

## BM1 Model:

#BM models do not calc theta estimative. We use only eingvalues to decided about model plausivity:

##Verify eingValue (if eingValue is  >= 0 == TRUE, if not == FALSE): 
eigval.BM1 <- c()  
for (i in SIMMAP) {
  eigval1 <- all(BM1[[i]]$eigval >=0) 
  eigval.BM1 <- c(eigval.BM1, eigval1)}

summary(eigval.BM1)

BM1Results <- cbind(SIMMAP, eigval.BM1)
BM1Results <-as.data.frame(BM1Results)

##Wich models are plausible?

BM1.ok <- BM1Results[(BM1Results$eigval.BM1 == 1),]

BM1.ok

## Insert NA values in AICc table when the models are not plausible:

BM1Results2 <- subset(BM1Results, !(SIMMAP %in% BM1.ok$SIMMAP))
fire_AICc$BM1_AICc[fire_AICc$SIMMAP %in% BM1Results2$SIMMAP] <- "NA"

##Make the same procedure for the other models: 

##BMS Model:

eigval.BMS <- c()  
for (i in SIMMAP) {
  eigval1 <- all(BMS[[i]]$eigval >=0) 
  eigval.BMS <- c(eigval.BMS, eigval1)}

summary(eigval.BMS)

##BMs models do not calc theta estimative. We use only eingvalues to decided about model plausivity: 

BMSResults <- cbind(SIMMAP, eigval.BMS)
BMSResults <-as.data.frame(BMSResults)

##Wich models are plausible?

BMS.ok <- BMSResults[(BMSResults$eigval.BMS == 1),]

summary(BMS.ok)

## Insert NA values in AICc table when the models are not plausible:

BMSResults2 <- subset(BMSResults, !(SIMMAP %in% BMS.ok$SIMMAP))

fire_AICc$BMS_AICc[fire_AICc$SIMMAP %in% BMSResults2$SIMMAP] <- "NA"

## OU1 model:

eigval.OU1 <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OU1[[i]]$eigval >=0) 
  eigval.OU1 <- c(eigval.OU1, eigval1)}

summary(eigval.OU1)

##Does the estimated theta make biological sense?
##Verify theta value (if theta value is between fire minimum and fire maximum == TRUE,if not == FALSE):

THETA.OU1 <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OU1[[i]]$theta[,1] >= min(fire$fire.mean)) && 
            all(OU1[[i]]$theta[,1] <= max(fire$fire.mean))
  THETA.OU1 <- c(THETA.OU1, THETA1)}


OU1Results <- cbind(SIMMAP, eigval.OU1, THETA.OU1)
OU1Results <-as.data.frame(OU1Results)

OU1.ok <- OU1Results[(OU1Results$eigval.OU1 == 1) & (OU1Results$THETA.OU1 ==1), ]

summary(OU1.ok)

OU1Results2 <- subset(OU1Results, !(SIMMAP %in% OU1.ok$SIMMAP))
fire_AICc$OU1_AICc[fire_AICc$SIMMAP %in% OU1Results2$SIMMAP] <- "NA"


##OUM model:

eigval.OUM <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUM[[i]]$eigval >=0) 
  eigval.OUM <- c(eigval.OUM, eigval1)}
summary(eigval.OUM)

THETA.OUM <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUM[[i]]$theta[,1] >= min(fire$fire.mean)) &&
            all(OUM[[i]]$theta[,1] <= max(fire$fire.mean))
  THETA.OUM <- c(THETA.OUM, THETA1)}

OUMResults <- cbind(SIMMAP, eigval.OUM, THETA.OUM)
OUMResults <-as.data.frame(OUMResults)

OUM.ok <- OUMResults[(OUMResults$eigval.OUM == 1) & (OUMResults$THETA.OUM ==1), ]

OUMResults2 <- subset(OUMResults, !(SIMMAP %in% OUM.ok$SIMMAP))
fire_AICc$OUM_AICc[fire_AICc$SIMMAP %in% OUMResults2$SIMMAP] <- "NA"

##OUMA model:

eigval.OUMA <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUMA[[i]]$eigval >=0) 
  eigval.OUMA <- c(eigval.OUMA, eigval1)}
summary(eigval.OUMA)

THETA.OUMA <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUMA[[i]]$theta[,1] >= min(fire$fire.mean)) && 
            all(OUMA[[i]]$theta[,1] <= max(fire$fire.mean))
  THETA.OUMA <- c(THETA.OUMA, THETA1)}

OUMAResults <- cbind(SIMMAP, eigval.OUMA, THETA.OUMA)
OUMAResults <-as.data.frame(OUMAResults)

OUMA.ok <- OUMAResults[(OUMAResults$eigval.OUMA == 1) & (OUMAResults$THETA.OUMA ==1), ]
OUMAResults2 <- subset(OUMAResults, !(SIMMAP %in% OUMA.ok$SIMMAP))
fire_AICc$OUMA_AICc[fire_AICc$SIMMAP %in% OUMAResults2$SIMMAP] <- "NA"

##OUMV model:

eigval.OUMV <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUMV[[i]]$eigval >=0) 
  eigval.OUMV <- c(eigval.OUMV, eigval1)}
summary(eigval.OUMV)


THETA.OUMV <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUMV[[i]]$theta[,1] >= min(fire$fire.mean)) &&
            all(OUMV[[i]]$theta[,1] <= max(fire$fire.mean))
  THETA.OUMV <- c(THETA.OUMV, THETA1)}

OUMVResults <- cbind(SIMMAP, eigval.OUMV, THETA.OUMV)
OUMVResults <-as.data.frame(OUMVResults)

OUMV.ok <- OUMVResults[(OUMVResults$eigval.OUMV == 1) & (OUMVResults$THETA.OUMV ==1), ]

OUMVResults2 <- subset(OUMVResults, !(SIMMAP %in% OUMV.ok$SIMMAP))
fire_AICc$OUMV_AICc[fire_AICc$SIMMAP %in% OUMVResults2$SIMMAP] <- "NA"

##OUMVA model:

eigval.OUMVA <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUMVA[[i]]$eigval >=0) 
  eigval.OUMVA <- c(eigval.OUMVA, eigval1)}
summary(eigval.OUMVA)


THETA.OUMVA <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUMVA[[i]]$theta[,1] >= min(fire$fire.mean)) &&
            all(OUMVA[[i]]$theta[,1] <= max(fire$fire.mean))
  THETA.OUMVA <- c(THETA.OUMVA, THETA1)}

OUMVAResults <- cbind(SIMMAP, eigval.OUMVA, THETA.OUMVA)
OUMVAResults <-as.data.frame(OUMVAResults)

OUMVA.ok <- OUMVAResults[(OUMVAResults$eigval.OUMVA == 1) & (OUMVAResults$THETA.OUMVA ==1), ]

OUMVAResults2 <- subset(OUMVAResults, !(SIMMAP %in% OUMVA.ok$SIMMAP))
fire_AICc$OUMVA_AICc[fire_AICc$SIMMAP %in% OUMVAResults2$SIMMAP] <- "NA"

fire_AICc <- mutate_all(fire_AICc, function(x) as.numeric(as.character(x))) 

summary(fire_AICc)

##Save Aicc data:

path
write.csv(fire_AICc, paste0(path, "fire/fire_AICc.csv"), row.names = F)

# -------------------------------------------
# White noise models and AICc values
# -------------------------------------------

fire.WT <- readRDS(paste0(path, "Results/04-OUwie_and_WN_models/fire/fire_fitWN.rds"))

###AICc: 

fire_WT_AICc  <- as.vector(unlist(lapply(fire.WT, function (x) x[[4]][[7]])))

fire_AICc2 <- add_column(fire_AICc, WT_AICc = fire_WT_AICc, .after = 1)

write.csv(fire_AICc2, paste0(path, "fire/fire_AICc_all.csv"), row.names = F)

# ------------------------------------------
# Delta AICc and model selection - FIRE
# ------------------------------------------
rm(list = ls())

##set working directory to the  data-paper folder:


fire.aicc <- read.csv("data-paper/fire/fire_AICc_all.csv", row.names=1)

##Calc DELTA AICc: 

##fire: 

fire.D.AICc <- data.frame(matrix(ncol = 8, nrow = 2000))
x <- c("WN_DAICc", "BM1_DAICc", "BMS_DAICc", "OU1_DAICc", "OUM_DAICc", "OUMA_DAICc", "OUMV_DAICc", "OUMVA_DAICc")
colnames(fire.D.AICc) <- x

for (i in 1:2000) { ## For each row in table: 
  
  min.i <- min(fire.aicc[i,], na.rm = T) ## select the min value in the row
  
  fire.D.AICc[i,] <- fire.aicc[i,] - min.i ##subtract the min value of each element in the row
  
}

##Save Delta Aicc data:

write.csv(fire.D.AICc, "data-paper/fire/fire_Delta_AICc_all.csv", row.names = T)

###Calc the frequency at which each  model was best fitted:

#Compare with all WN model fits:

WN.DELTA <- drop_na(fire.D.AICc, WN_DAICc)

fittings <- WN.DELTA %>% 
  #select(1:7) %>% 
  is.na %>% 
  `!` %>% 
  colSums

model.names <- c("WN", "BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
best        <- c(0,0,0,0,0,0,0,0)  
tie         <- c(0,0,0,0,0,0,0,0) 


for (i in 1:2000) 
{
  hope <- which(WN.DELTA[i,] <=2) ##How many values in each row are less than 2?
  
  if (length(hope) == 1) ## If only 1 value is less than 2:  
    
  {
    best[hope[1]] = best[hope[1]] + 1 ##Add 1 at the position where the value is greater than 2 in the "best" vector.
    
  } else if (length(hope) > 1) ##If there is more than one value less than 2 
  {
    for (j in 1:length(hope))
    {
      tie[hope[j]] = tie[hope[j]] + 1 ##add 1 at each position in the "TIE" vector
    }
  }
}

best
tie

frequency <- (best*100)/ length(WN.DELTA$WN_DAICc) 

resultsWN <- rbind(model.names, fittings, best, tie, frequency)

##Make the same procedure for the other models: 

#Compare with all BM1 model fits:

BM1DELTA <- drop_na(fire.D.AICc, BM1_DAICc)

fittings <- BM1DELTA %>% 
  #select(1:7) %>% 
  is.na %>% 
  `!` %>% 
  colSums

model.names <- c("WN", "BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
best        <- c(0,0,0,0,0,0,0,0)  
tie         <- c(0,0,0,0,0,0,0,0) 

for (i in 1:2000)  
{
  hope <- which(BM1DELTA[i,] <=2) 
  
  if (length(hope) == 1) 
    
  {
    best[hope[1]] = best[hope[1]] + 1 
    
  } else if (length(hope) > 1) 
  {
    for (j in 1:length(hope))
    {
      tie[hope[j]] = tie[hope[j]] + 1 
    }
  }
}

best
tie


frequency <- (best*100)/ length(BM1DELTA$BM1_DAICc) 

resultsBM1 <- rbind(model.names, fittings, best, tie, frequency) 

#Compare with all BMS model fits:

BMSDELTA <- drop_na(fire.D.AICc, BMS_DAICc)


fittings <- BMSDELTA %>% 
  #select(1:7) %>% 
  is.na %>% 
  `!` %>% 
  colSums

model.names <- c("WN", "BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
best        <- c(0,0,0,0,0,0,0,0)  
tie         <- c(0,0,0,0,0,0,0,0) 


for (i in 1:2000)  
{
  hope <- which(BMSDELTA[i,] <=2) 
  
  if (length(hope) == 1)  
    
  {
    best[hope[1]] = best[hope[1]] + 1 
    
  } else if (length(hope) > 1)  
  {
    for (j in 1:length(hope))
    {
      tie[hope[j]] = tie[hope[j]] + 1 
    }
  }
}

best
tie



frequency <- (best*100)/ length(BMSDELTA$BMS_DAICc) 

resultsBMS <- rbind(model.names, fittings, best, tie, frequency) 

#Compare with all OU1 model fits:

OU1DELTA <- drop_na(fire.D.AICc, OU1_DAICc)


fittings <- OU1DELTA %>% 
  #select(1:7) %>% 
  is.na %>% 
  `!` %>% 
  colSums

model.names <- c("WN", "BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
best        <- c(0,0,0,0,0,0,0,0) 
tie         <- c(0,0,0,0,0,0,0,0) 


for (i in 1:2000)  
{
  hope <- which(OU1DELTA[i,] <=2) 
  
  if (length(hope) == 1) 
    
  {
    best[hope[1]] = best[hope[1]] + 1 
    
  } else if (length(hope) > 1) 
  {
    for (j in 1:length(hope))
    {
      tie[hope[j]] = tie[hope[j]] + 1 
    }
  }
}

best
tie



frequency <- (best*100)/ length(OU1DELTA$OU1_DAICc) 

resultsOU1 <- rbind(model.names, fittings, best, tie, frequency) 

#Compare with all OUM model fits:

OUMDELTA <- drop_na(fire.D.AICc, OUM_DAICc)

fittings <- OUMDELTA %>% 
  #select(1:7) %>% 
  is.na %>% 
  `!` %>% 
  colSums

model.names <- c("WN", "BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
best        <- c(0,0,0,0,0,0,0,0)  
tie         <- c(0,0,0,0,0,0,0,0) 


for (i in 1:2000) 
{
  hope <- which(OUMDELTA[i,] <=2) 
  
  if (length(hope) == 1)  
    
  {
    best[hope[1]] = best[hope[1]] + 1 
    
  } else if (length(hope) > 1)  
  {
    for (j in 1:length(hope))
    {
      tie[hope[j]] = tie[hope[j]] + 1 
    }
  }
}

best
tie

frequency <- (best*100)/ length(OUMDELTA$OUM_DAICc) 

resultsOUM <- rbind(model.names, fittings, best, tie, frequency) 

#Compare with all OUMA model fits:

OUMADELTA <- drop_na(fire.D.AICc, OUMA_DAICc)

fittings <- OUMADELTA %>% 
  #select(1:7) %>% 
  is.na %>% 
  `!` %>% 
  colSums

model.names <- c("WN", "BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
best        <- c(0,0,0,0,0,0,0,0) 
tie         <- c(0,0,0,0,0,0,0,0) 

for (i in 1:2000)  
{
  hope <- which(OUMADELTA[i,] <=2) 
  
  if (length(hope) == 1)  
    
  {
    best[hope[1]] = best[hope[1]] + 1 
    
  } else if (length(hope) > 1)  
  {
    for (j in 1:length(hope))
    {
      tie[hope[j]] = tie[hope[j]] + 1 
    }
  }
}

best
tie



frequency <- (best*100)/ length(OUMADELTA$OUMA_DAICc) 

resultsOUMA <- rbind(model.names, fittings, best, tie, frequency) 

#Compare with all OUMV model fits:

OUMVDELTA <- drop_na(fire.D.AICc, OUMV_DAICc)


fittings <- OUMVDELTA %>% 
  #select(1:7) %>% 
  is.na %>% 
  `!` %>% 
  colSums

model.names <- c("WN", "BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
best        <- c(0,0,0,0,0,0,0,0)  
tie         <- c(0,0,0,0,0,0,0,0) 


for (i in 1:2000)  
{
  hope <- which(OUMVDELTA[i,] <=2) 
  
  if (length(hope) == 1)  
    
  {
    best[hope[1]] = best[hope[1]] + 1 
    
  } else if (length(hope) > 1) 
  {
    for (j in 1:length(hope))
    {
      tie[hope[j]] = tie[hope[j]] + 1 
    }
  }
}

best
tie

frequency <- (best*100)/ length(OUMVDELTA$OUMV_DAICc) 

resultsOUMV <- rbind(model.names, fittings, best, tie, frequency) 

#Compare with all OUMVA model fits:

OUMVADELTA <- drop_na(fire.D.AICc, OUMVA_DAICc)

fittings <- OUMVADELTA %>% 
  #select(1:7) %>% 
  is.na %>% 
  `!` %>% 
  colSums

model.names <- c("WN", "BM1", "BMS", "OU1", "OUM", "OUMA", "OUMV", "OUMVA")
best        <- c(0,0,0,0,0,0,0,0)  
tie         <- c(0,0,0,0,0,0,0,0)


for (i in 1:2000) 
{
  hope <- which(OUMVADELTA[i,] <=2) 
  
  if (length(hope) == 1) 
    
  {
    best[hope[1]] = best[hope[1]] + 1 
    
  } else if (length(hope) > 1)  
  {
    for (j in 1:length(hope))
    {
      tie[hope[j]] = tie[hope[j]] + 1 
    }
  }
}

best
tie

frequency <- (best*100)/ length(OUMVADELTA$OUMVA_DAICc) 

resultsOUMVA <- rbind(model.names, fittings, best, tie, frequency) 

##Join all results:

resultadoVF <- rbind(resultsWN,resultsBM1, resultsBMS, resultsOU1, resultsOUM,resultsOUMA, resultsOUMV, resultsOUMVA)

write.csv(resultadoVF, "/data-paper/fire/fire_best_fit_all.csv", row.names = T)

# ------------------------------------------------------------
# Plot FIRE theta estimate values - Only the best models 
# ------------------------------------------------------------

rm(list = ls())
##set working directory to the  data-paper folder:

# The best model to fire is OUMV model: 

OUMV <- readRDS("data-paper/Results/04-OUwie_and_WN_models/fire/fire_fitOUMV.rds")

##Delta AICc data: 

fire.DAicc <- read.csv("data-paper/fire/fire_Delta_AICc_all.csv")

##Remove OUMV models with bad fit:

validos <- drop_na(fire.DAicc, OUMV_DAICc)

hope <- which(validos[1, 2:9] < 2)

best <- c()

for (i in 1:length(validos$OUMV_DAICc)) { 
  
  hope <- which(validos[i, 2:9] < 2) ##How many values in each row are less than 2?
  
  if (length(hope) == 1 && hope == 7) {### If only 1 value is less than 2 and it is in column 7 (OUMV model column)
    
    x <- validos[i, 1] ## X = SIMMAP ID
    
    best <- c(best, x)
    
  }  
}      

best

##Filter the best OUMV models:

OUMV.new <-list()

for (i in 1:2000){
  
  if (!(i %in% best)) next
  
  OUMV.new[[(length(OUMV.new) + 1)]] <- OUMV[[i]] 
  
}

##Theta value: 

x <- as.data.frame(OUMV.new[[1]]$tot.states)
OUMV.new[[1]]$theta

thetavf <- data.frame()

for (i in 1:length(OUMV.new)){
  
  states <-OUMV.new[[i]]$tot.states
  states <- as.data.frame(states)
  
  theta <- OUMV.new[[i]]$theta
  as.data.frame(theta)
  
  theta <- cbind(states, theta)
  
  thetavf <- rbind(thetavf, theta)
  
}

summary(thetavf)

#transform log data into original data (remove log):

thetavf2 <- exp(thetavf[,-1])
thetavf2 <- thetavf2 - 1

thetavf2 <- cbind(thetavf$states,thetavf2)
names(thetavf2)[1] <- "states"
names(thetavf2)[2] <- "theta"


summary(thetavf2)

##Transform theta values in 2 decimal digits: 
#install.packages("peRspective")
library("peRspective")
thetavf2$theta <- specify_decimal(thetavf2$theta, 2)
thetavf2$theta <- as.numeric(thetavf2$theta)

summary(thetavf2)

thero <- filter(thetavf2, states=="therophyte")
summary(thero$theta)

count(thero$theta == 0.00)
count(thero$theta == 0.01)
count(thero$theta == 0.02)
count(thero$theta == 0.03)
count(thero$theta == 0.04)
count(thero$theta == 0.05)
count(thero$theta == 0.06)

chama <- filter(thetavf2, states=="chamaephyte")
summary(chama$theta)

count(chama$theta == 0.00)
count(chama$theta == 0.01)
count(chama$theta == 0.02)
count(chama$theta == 0.03)
count(chama$theta == 0.04)
count(chama$theta == 0.05)
count(chama$theta == 0.06)

phane <- filter(thetavf2, states=="phanerophyte")
summary(phane$theta)

count(phane$theta == 0.00)
count(phane$theta == 0.01)
count(phane$theta == 0.02)
count(phane$theta == 0.03)
count(phane$theta == 0.04)
count(phane$theta == 0.05)
count(phane$theta == 0.06)

geo <- filter(thetavf2, states=="geophyte")
summary(geo$theta)

count(geo$theta == 0.00)
count(geo$theta == 0.01)
count(geo$theta == 0.02)
count(geo$theta == 0.03)
count(geo$theta == 0.04)
count(geo$theta == 0.05)
count(geo$theta == 0.06)



##Summary Theta Values by states: 
library("Rmisc")

fire.theta <- summarySE(thetavf2, measurevar= "theta", groupvars = "states", na.rm = T)
fire.theta

thetavf2 <- as.data.frame(thetavf2)
summary(thetavf2)


unique(thetavf2$theta)

#States order: 

fire.theta$states <- factor(fire.theta$states, levels=c("therophyte", "chamaephyte", "phanerophyte", "geophyte"))

thetavf2$states <- factor(thetavf2$state, levels=c("therophyte", "chamaephyte", "phanerophyte", "geophyte"))

##Plot fire theta: 
library(grid)

grob <- grobTree(textGrob("(E)", x=0.03,  y=0.96, hjust=0,
                          gp=gpar(col="black", fontsize=8, fontface="bold")))


ggplot(data = fire.theta, mapping = aes(x=states, y= theta)) + 
  geom_point(size = 1.5) +
  geom_jitter(data = thetavf2, mapping = aes(x= states, y= theta, color = states), size = 0.4, alpha = 0.4, height = 0.000, width = 0.1) +
  scale_color_manual(values=c("#7743DB","#FFC93C","#B1AEB9","#DC267F")) +
  geom_errorbar(aes(ymin=theta-sd, ymax=theta+sd), width=.150) +
  geom_point(size = 1.5) +
  annotation_custom(grob) +
  scale_y_continuous(name = "?? estimate in OUMV models \n Fire Frequency (average by year)",limits = c(0, 0.06), n.breaks = 6) +
  
  theme_bw() +
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text=element_text(size= 8.5 , color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title=element_text(size=12))

ggsave(width = 4, height = 4, dpi = 600, filename = "fire_OUMV_theta_BEST.png")
