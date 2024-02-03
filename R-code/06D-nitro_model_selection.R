# -----------------------------------------------------------
# AICc values, model plausibility and
# Delta AICc calculation 
# -----------------------------------------------------------

# -----------------------------------------------------------
# AICc values and model plausibility - OUwie models - nitro
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

nitro <- read.csv("nitro/nitrolog_data.csv")

##Model Results:  

BM1   <- readRDS("nitro/nitro_fitBM1.rds")
BMS   <- readRDS("nitro/nitro_fitBMS.rds")
OU1   <- readRDS("nitro/nitro_fitOU1.rds")
OUM   <- readRDS("nitro/nitro_fitOUM.rds")
OUMA  <- readRDS("nitro/nitro_fitOUMA.rds")
OUMV  <- readRDS("nitro/nitro_fitOUMV.rds")
OUMVA <- readRDS("nitro/nitro_fitOUMVA.rds")

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

nitro_AICc <- cbind(SIMMAP, BM1_AICc, BMS_AICc, OU1_AICc, OUM_AICc, OUMA_AICc, OUMV_AICc, OUMVA_AICc)
nitro_AICc <- as.data.frame(nitro_AICc)

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
nitro_AICc$BM1_AICc[nitro_AICc$SIMMAP %in% BM1Results2$SIMMAP] <- "NA"

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

nitro_AICc$BMS_AICc[nitro_AICc$SIMMAP %in% BMSResults2$SIMMAP] <- "NA"

## OU1 model:

eigval.OU1 <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OU1[[i]]$eigval >=0) 
  eigval.OU1 <- c(eigval.OU1, eigval1)}

summary(eigval.OU1)

##Does the estimated theta make biological sense?
##Verify theta value (if theta value is between nitro minimum and nitro maximum == TRUE,if not == FALSE):

THETA.OU1 <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OU1[[i]]$theta[,1] >= min(nitro$nitro.mean)) && 
            all(OU1[[i]]$theta[,1] <= max(nitro$nitro.mean))
  THETA.OU1 <- c(THETA.OU1, THETA1)}


OU1Results <- cbind(SIMMAP, eigval.OU1, THETA.OU1)
OU1Results <-as.data.frame(OU1Results)

OU1.ok <- OU1Results[(OU1Results$eigval.OU1 == 1) & (OU1Results$THETA.OU1 ==1), ]

summary(OU1.ok)

OU1Results2 <- subset(OU1Results, !(SIMMAP %in% OU1.ok$SIMMAP))
nitro_AICc$OU1_AICc[nitro_AICc$SIMMAP %in% OU1Results2$SIMMAP] <- "NA"


##OUM model:

eigval.OUM <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUM[[i]]$eigval >=0) 
  eigval.OUM <- c(eigval.OUM, eigval1)}
summary(eigval.OUM)

THETA.OUM <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUM[[i]]$theta[,1] >= min(nitro$nitro.mean)) &&
            all(OUM[[i]]$theta[,1] <= max(nitro$nitro.mean))
  THETA.OUM <- c(THETA.OUM, THETA1)}

OUMResults <- cbind(SIMMAP, eigval.OUM, THETA.OUM)
OUMResults <-as.data.frame(OUMResults)

OUM.ok <- OUMResults[(OUMResults$eigval.OUM == 1) & (OUMResults$THETA.OUM ==1), ]

OUMResults2 <- subset(OUMResults, !(SIMMAP %in% OUM.ok$SIMMAP))
nitro_AICc$OUM_AICc[nitro_AICc$SIMMAP %in% OUMResults2$SIMMAP] <- "NA"

##OUMA model:

eigval.OUMA <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUMA[[i]]$eigval >=0) 
  eigval.OUMA <- c(eigval.OUMA, eigval1)}
summary(eigval.OUMA)

THETA.OUMA <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUMA[[i]]$theta[,1] >= min(nitro$nitro.mean)) && 
            all(OUMA[[i]]$theta[,1] <= max(nitro$nitro.mean))
  THETA.OUMA <- c(THETA.OUMA, THETA1)}

OUMAResults <- cbind(SIMMAP, eigval.OUMA, THETA.OUMA)
OUMAResults <-as.data.frame(OUMAResults)

OUMA.ok <- OUMAResults[(OUMAResults$eigval.OUMA == 1) & (OUMAResults$THETA.OUMA ==1), ]
OUMAResults2 <- subset(OUMAResults, !(SIMMAP %in% OUMA.ok$SIMMAP))
nitro_AICc$OUMA_AICc[nitro_AICc$SIMMAP %in% OUMAResults2$SIMMAP] <- "NA"

##OUMV model:

eigval.OUMV <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUMV[[i]]$eigval >=0) 
  eigval.OUMV <- c(eigval.OUMV, eigval1)}
summary(eigval.OUMV)


THETA.OUMV <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUMV[[i]]$theta[,1] >= min(nitro$nitro.mean)) &&
            all(OUMV[[i]]$theta[,1] <= max(nitro$nitro.mean))
  THETA.OUMV <- c(THETA.OUMV, THETA1)}

OUMVResults <- cbind(SIMMAP, eigval.OUMV, THETA.OUMV)
OUMVResults <-as.data.frame(OUMVResults)

OUMV.ok <- OUMVResults[(OUMVResults$eigval.OUMV == 1) & (OUMVResults$THETA.OUMV ==1), ]

OUMVResults2 <- subset(OUMVResults, !(SIMMAP %in% OUMV.ok$SIMMAP))
nitro_AICc$OUMV_AICc[nitro_AICc$SIMMAP %in% OUMVResults2$SIMMAP] <- "NA"

##OUMVA model:

eigval.OUMVA <- c()  
for (i in SIMMAP) {
  eigval1 <- all(OUMVA[[i]]$eigval >=0) 
  eigval.OUMVA <- c(eigval.OUMVA, eigval1)}
summary(eigval.OUMVA)


THETA.OUMVA <- c()  
for (i in SIMMAP) {
  THETA1 <- all(OUMVA[[i]]$theta[,1] >= min(nitro$nitro.mean)) &&
            all(OUMVA[[i]]$theta[,1] <= max(nitro$nitro.mean))
  THETA.OUMVA <- c(THETA.OUMVA, THETA1)}

OUMVAResults <- cbind(SIMMAP, eigval.OUMVA, THETA.OUMVA)
OUMVAResults <-as.data.frame(OUMVAResults)

OUMVA.ok <- OUMVAResults[(OUMVAResults$eigval.OUMVA == 1) & (OUMVAResults$THETA.OUMVA ==1), ]

OUMVAResults2 <- subset(OUMVAResults, !(SIMMAP %in% OUMVA.ok$SIMMAP))
nitro_AICc$OUMVA_AICc[nitro_AICc$SIMMAP %in% OUMVAResults2$SIMMAP] <- "NA"

nitro_AICc <- mutate_all(nitro_AICc, function(x) as.numeric(as.character(x))) 

summary(nitro_AICc)

##Save Aicc data:

path
write.csv(nitro_AICc, paste0(path, "nitro/nitro_AICc.csv"), row.names = F)

# -------------------------------------------
# White noise models and AICc values
# -------------------------------------------

nitro.WT <- readRDS(paste0(path, "Results/04-OUwie_and_WN_models/nitro/nitro_fitWN.rds"))

###AICc:

library(tidyverse) 

nitro_WT_AICc  <- as.vector(unlist(lapply(nitro.WT, function (x) x[[4]][[7]])))

nitro_AICc2 <- add_column(nitro_AICc, WT_AICc = nitro_WT_AICc, .after = 1)

write.csv(nitro_AICc2, paste0(path, "nitro/nitro_AICc_all.csv"), row.names = F)

# ------------------------------------------
# Delta AICc and model selection - nitro
# ------------------------------------------
rm(list = ls())

##set working directory to the  data-paper folder:


nitro.aicc <- read.csv("/data-paper/nitro/nitro_AICc_all.csv", row.names=1)

##Calc DELTA AICc: 

##nitro: 

nitro.D.AICc <- data.frame(matrix(ncol = 8, nrow = 2000))
x <- c("WN_DAICc", "BM1_DAICc", "BMS_DAICc", "OU1_DAICc", "OUM_DAICc", "OUMA_DAICc", "OUMV_DAICc", "OUMVA_DAICc")
colnames(nitro.D.AICc) <- x

for (i in 1:2000) { ## For each row in table: 
  
  min.i <- min(nitro.aicc[i,], na.rm = T) ## select the min value in the row
  
  nitro.D.AICc[i,] <- nitro.aicc[i,] - min.i ##subtract the min value of each element in the row
  
}

##Save Delta Aicc data:

write.csv(nitro.D.AICc, "/data-paper/nitro/nitro_Delta_AICc_all.csv", row.names = T)

###Calc the frequency at which each  model was best fitted:

#Compare with all WN model fits:

WN.DELTA <- drop_na(nitro.D.AICc, WN_DAICc)

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

BM1DELTA <- drop_na(nitro.D.AICc, BM1_DAICc)

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

BMSDELTA <- drop_na(nitro.D.AICc, BMS_DAICc)


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

OU1DELTA <- drop_na(nitro.D.AICc, OU1_DAICc)


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

OUMDELTA <- drop_na(nitro.D.AICc, OUM_DAICc)

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

OUMADELTA <- drop_na(nitro.D.AICc, OUMA_DAICc)

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

OUMVDELTA <- drop_na(nitro.D.AICc, OUMV_DAICc)


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

OUMVADELTA <- drop_na(nitro.D.AICc, OUMVA_DAICc)

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

write.csv(resultadoVF, "/data-paper/nitro/nitro_best_fit_all.csv", row.names = T)

# ------------------------------------------------------------
# Plot nitro theta estimate values - Only the best models 
# ------------------------------------------------------------

rm(list = ls())
##set working directory to the  data-paper folder:

# The best model to nitro is OUM model: 

OUM <- readRDS("/data-paper/Results/04-OUwie_and_WN_models/nitro/nitro_fitOUM.rds")

##Delta AICc data: 

nitro.DAicc <- read.csv("/data-paper/nitro/nitro_Delta_AICc_all.csv")

##Remove OUM models with bad fit:

validos <- drop_na(nitro.DAicc, OUM_DAICc)

hope <- which(validos[1, 2:9] < 2)

best <- c()

for (i in 1:length(validos$OUM_DAICc)) { 
  
  hope <- which(validos[i, 2:9] < 2) ##How many values in each row are less than 2?
  
  if (length(hope) == 1 && hope == 5) {### If only 1 value is less than 2 and it is in column 7 (OUM model column)
    
    x <- validos[i, 1] ## X = SIMMAP ID
    
    best <- c(best, x)
    
  }  
}      

best

##Filter the best OUM models:

OUM.new <-list()

for (i in 1:2000){
  
  if (!(i %in% best)) next
  
  OUM.new[[(length(OUM.new) + 1)]] <- OUM[[i]] 
  
}

##Theta value: 

x <- as.data.frame(OUM.new[[1]]$tot.states)
OUM.new[[1]]$theta

thetavf <- data.frame()

for (i in 1:length(OUM.new)){
  
  states <-OUM.new[[i]]$tot.states
  states <- as.data.frame(states)
  
  theta <- OUM.new[[i]]$theta
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


##Summary Theta Values by states: 
library("Rmisc")

nitro.theta <- summarySE(thetavf2, measurevar= "theta", groupvars = "states", na.rm = T)
nitro.theta

thetavf2 <- as.data.frame(thetavf2)
summary(thetavf2)


#States order: 

nitro.theta$states <- factor(nitro.theta$states, levels=c("therophyte", "chamaephyte", "phanerophyte", "geophyte"))

thetavf2$states <- factor(thetavf2$state, levels=c("therophyte", "chamaephyte", "phanerophyte", "geophyte"))

##Plot nitro theta: 
library(grid)

grob <- grobTree(textGrob("(C)", x=0.03,  y=0.96, hjust=0,
                          gp=gpar(col="black", fontsize=8, fontface="bold")))


ggplot(data = nitro.theta, mapping = aes(x=states, y= theta)) +
  geom_point(size = 1.5) +
  geom_jitter(data = thetavf2, mapping = aes(x= states, y= theta, color = states), size = 0.4, alpha = 0.4, height = 0.00, width = 0.1) +
  scale_color_manual(values=c("#7743DB","#FFC93C","#B1AEB9","#DC267F")) +
  geom_errorbar(aes(ymin=theta-sd, ymax=theta+sd), width=.150) +
  geom_point(size = 1.5) +
  annotation_custom(grob) +
  scale_y_continuous(name = "?? estimate in OUM models \n Soil nitrogen (cg/kg) depht 0-30cm", n.breaks = 10) + 
  
  theme_bw() +
  theme(axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text=element_text(size= 8.5 , color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.title=element_text(size=12))

ggsave(width = 4, height = 4, dpi = 600, filename = "nitro_OUM_theta_BEST.png")
