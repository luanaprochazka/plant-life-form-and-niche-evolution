# -----------------------------------------------------------
# AICc values, model plausibility and
# Delta AICc calculation 
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

bio17 <- read.csv("bio17/bio17log_data.csv")

##Model Results: 

BM1   <- readRDS("bio17/bio17_fitBM1.rds")
BMS   <- readRDS("bio17/bio17_fitBMS.rds")
OU1   <- readRDS("bio17/bio17_fitOU1.rds")
OUM   <- readRDS("bio17/bio17_fitOUM.rds")
OUMA  <- readRDS("bio17/bio17_fitOUMA.rds")
OUMV  <- readRDS("bio17/bio17_fitOUMV.rds")
OUMVA <- readRDS("bio17/bio17_fitOUMVA.rds")

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

bio17_AICc <- cbind(SIMMAP, BM1_AICc, BMS_AICc, OU1_AICc, OUM_AICc, OUMA_AICc, OUMV_AICc, OUMVA_AICc)
bio17_AICc <- as.data.frame(bio17_AICc)

##Check the quality adjust for each model: 

## BM1 Model:

##BM models do not calc theta estimative. We use only eingvalues to decided about model plausivity:

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
bio17_AICc$BM1_AICc[bio17_AICc$SIMMAP %in% BM1Results2$SIMMAP] <- "NA"

##Make the same procedure for the other models: 

##BMS Model:

##BMs models do not calc theta estimative. We use only eingvalues to decided about model plausivity:

eigval.BMS <- c()  
for (i in SIMMAP) {
  eigval1 <- all(BMS[[i]]$eigval >=0) 
  eigval.BMS <- c(eigval.BMS, eigval1)}

summary(eigval.BMS)

BMSResults <- cbind(SIMMAP, eigval.BMS)
BMSResults <-as.data.frame(BMSResults)

##Wich models are plausible?

BMS.ok <- BMSResults[(BMSResults$eigval.BMS == 1),]

summary(BMS.ok)

## Insert NA values in AICc table when the models are not plausible:

BMSResults2 <- subset(BMSResults, !(SIMMAP %in% BMS.ok$SIMMAP))

bio17_AICc$BMS_AICc[bio17_AICc$SIMMAP %in% BMSResults2$SIMMAP] <- "NA"

## OU1 model:

eigval.OU1 <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OU1[[i]]$eigval >=0) 
  eigval.OU1 <- c(eigval.OU1, eigval1)}

summary(eigval.OU1)

##Does the estimated theta make biological sense?
##Verify theta value (if theta value is between bio17 minimum and bio17 maximum == TRUE,if not == FALSE):

THETA.OU1 <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OU1[[i]]$theta[,1] >= min(bio17$bio17.mean)) && 
            all(OU1[[i]]$theta[,1] <= max(bio17$bio17.mean))
  THETA.OU1 <- c(THETA.OU1, THETA1)}


OU1Results <- cbind(SIMMAP, eigval.OU1, THETA.OU1)
OU1Results <-as.data.frame(OU1Results)

OU1.ok <- OU1Results[(OU1Results$eigval.OU1 == 1) & (OU1Results$THETA.OU1 ==1), ]

summary(OU1.ok)

OU1Results2 <- subset(OU1Results, !(SIMMAP %in% OU1.ok$SIMMAP))
bio17_AICc$OU1_AICc[bio17_AICc$SIMMAP %in% OU1Results2$SIMMAP] <- "NA"


##OUM model:

eigval.OUM <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUM[[i]]$eigval >=0) 
  eigval.OUM <- c(eigval.OUM, eigval1)}
summary(eigval.OUM)

THETA.OUM <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUM[[i]]$theta[,1] >= min(bio17$bio17.mean)) &&
            all(OUM[[i]]$theta[,1] <= max(bio17$bio17.mean))
  THETA.OUM <- c(THETA.OUM, THETA1)}

OUMResults <- cbind(SIMMAP, eigval.OUM, THETA.OUM)
OUMResults <-as.data.frame(OUMResults)

OUM.ok <- OUMResults[(OUMResults$eigval.OUM == 1) & (OUMResults$THETA.OUM ==1), ]

OUMResults2 <- subset(OUMResults, !(SIMMAP %in% OUM.ok$SIMMAP))
bio17_AICc$OUM_AICc[bio17_AICc$SIMMAP %in% OUMResults2$SIMMAP] <- "NA"

##OUMA model:

eigval.OUMA <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUMA[[i]]$eigval >=0) 
  eigval.OUMA <- c(eigval.OUMA, eigval1)}
summary(eigval.OUMA)

THETA.OUMA <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUMA[[i]]$theta[,1] >= min(bio17$bio17.mean)) && 
            all(OUMA[[i]]$theta[,1] <= max(bio17$bio17.mean))
  THETA.OUMA <- c(THETA.OUMA, THETA1)}

OUMAResults <- cbind(SIMMAP, eigval.OUMA, THETA.OUMA)
OUMAResults <-as.data.frame(OUMAResults)

OUMA.ok <- OUMAResults[(OUMAResults$eigval.OUMA == 1) & (OUMAResults$THETA.OUMA ==1), ]
OUMAResults2 <- subset(OUMAResults, !(SIMMAP %in% OUMA.ok$SIMMAP))
bio17_AICc$OUMA_AICc[bio17_AICc$SIMMAP %in% OUMAResults2$SIMMAP] <- "NA"

##OUMV model:

eigval.OUMV <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUMV[[i]]$eigval >=0) 
  eigval.OUMV <- c(eigval.OUMV, eigval1)}
summary(eigval.OUMV)


THETA.OUMV <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUMV[[i]]$theta[,1] >= min(bio17$bio17.mean)) &&
            all(OUMV[[i]]$theta[,1] <= max(bio17$bio17.mean))
  THETA.OUMV <- c(THETA.OUMV, THETA1)}

OUMVResults <- cbind(SIMMAP, eigval.OUMV, THETA.OUMV)
OUMVResults <-as.data.frame(OUMVResults)

OUMV.ok <- OUMVResults[(OUMVResults$eigval.OUMV == 1) & (OUMVResults$THETA.OUMV ==1), ]

OUMVResults2 <- subset(OUMVResults, !(SIMMAP %in% OUMV.ok$SIMMAP))
bio17_AICc$OUMV_AICc[bio17_AICc$SIMMAP %in% OUMVResults2$SIMMAP] <- "NA"

##OUMVA model:

eigval.OUMVA <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUMVA[[i]]$eigval >=0) 
  eigval.OUMVA <- c(eigval.OUMVA, eigval1)}
summary(eigval.OUMVA)


THETA.OUMVA <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUMVA[[i]]$theta[,1] >= min(bio17$bio17.mean)) &&
            all(OUMVA[[i]]$theta[,1] <= max(bio17$bio17.mean))
  THETA.OUMVA <- c(THETA.OUMVA, THETA1)}

OUMVAResults <- cbind(SIMMAP, eigval.OUMVA, THETA.OUMVA)
OUMVAResults <-as.data.frame(OUMVAResults)

OUMVA.ok <- OUMVAResults[(OUMVAResults$eigval.OUMVA == 1) & (OUMVAResults$THETA.OUMVA ==1), ]

OUMVAResults2 <- subset(OUMVAResults, !(SIMMAP %in% OUMVA.ok$SIMMAP))
bio17_AICc$OUMVA_AICc[bio17_AICc$SIMMAP %in% OUMVAResults2$SIMMAP] <- "NA"

bio17_AICc <- mutate_all(bio17_AICc, function(x) as.numeric(as.character(x))) 

summary(bio17_AICc)

##Save Aicc data:

path
write.csv(bio17_AICc, paste0(path, "Results/04-OUwie_and_WN_models/bio17/bio17_AICc.csv"), row.names = F)

# -------------------------------------------
# White noise models and AICc values
# -------------------------------------------

bio17.WT <- readRDS(paste0(path, "Results/04-OUwie_and_WN_models/bio17/bio17_fitWN.rds"))

###AICc: 

bio17_WT_AICc  <- as.vector(unlist(lapply(bio17.WT, function (x) x[[4]][[7]])))

bio17_AICc2 <- add_column(bio17_AICc, WT_AICc = bio17_WT_AICc, .after = 1)

write.csv(bio17_AICc2, paste0(path, "Results/04-OUwie_and_WN_models/bio17/bio17_AICc_all.csv"), row.names = F)

# ------------------------------------------
# Delta AICc and model selection - bio17
# ------------------------------------------
rm(list = ls())

##set working directory to the  data-paper folder:

bio17.aicc <- read.csv(paste0(path, "/Results/04-OUwie_and_WN_models/bio17/bio17_AICc_all.csv", row.names=1))

##Calc DELTA AICc: 

##bio17: 

bio17.D.AICc <- data.frame(matrix(ncol = 8, nrow = 2000))
x <- c("WN_DAICc", "BM1_DAICc", "BMS_DAICc", "OU1_DAICc", "OUM_DAICc", "OUMA_DAICc", "OUMV_DAICc", "OUMVA_DAICc")
colnames(bio17.D.AICc) <- x

for (i in 1:2000) { ## For each row in table: 
  
  min.i <- min(bio17.aicc[i,], na.rm = T) ## select the min value in the row
  
  bio17.D.AICc[i,] <- bio17.aicc[i,] - min.i ##subtract the min value of each element in the row
  
}

##Save Delta Aicc data:

write.csv(bio17.D.AICc, "Results/04-OUwie_and_WN_models/bio17//bio17_Delta_AICc_all.csv", row.names = T)

###Calc the frequency at which each  model was best fitted:

#Compare with all WN model fits:

WN.DELTA <- drop_na(bio17.D.AICc, WN_DAICc)

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

BM1DELTA <- drop_na(bio17.D.AICc, BM1_DAICc)

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

BMSDELTA <- drop_na(bio17.D.AICc, BMS_DAICc)


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

OU1DELTA <- drop_na(bio17.D.AICc, OU1_DAICc)


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

OUMDELTA <- drop_na(bio17.D.AICc, OUM_DAICc)

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

OUMADELTA <- drop_na(bio17.D.AICc, OUMA_DAICc)

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

OUMVDELTA <- drop_na(bio17.D.AICc, OUMV_DAICc)


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

OUMVADELTA <- drop_na(bio17.D.AICc, OUMVA_DAICc)

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

write.csv(resultadoVF, paste0(path, "Results/04-OUwie_and_WN_models/bio17/bio17_best_fit_all.csv"), row.names = T)
