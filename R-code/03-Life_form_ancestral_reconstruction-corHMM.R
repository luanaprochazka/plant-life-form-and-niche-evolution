# ---------------------------------------------------
# CorHMM analysis and Plant Life form ancestral reconstruction
# ---------------------------------------------------

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
library(GMCM)
library(corHMM)

##set working directory to the data folder in the data-paper:

path <- "C:/data-paper"

##Data: 

##Phylogeny: 

phylo1 <- read.tree(paste0(path,"/data/tree_vasconcelos41spp.txt"))

phylo2 <- read.tree(paste0(path,"/data/tree200_vasconcelos41spp.txt"))

phylo1$tip.label

##life form classification: 

data <- read.csv(paste0(path,"/data/life_form.csv"), header = T)
data$life_form

#In corHMM discrete traits need be numeric: 

data$life_form <- gsub("therophyte", 1, data$life_form)
data$life_form <- gsub("chamaephyte", 2, data$life_form)

data$life_form <- gsub("phanerophyte", 3, data$life_form)
data$life_form <- gsub("geophyte", 4, data$life_form)

data$life_form <- as.numeric(data$life_form)

###########################################################
##Part I: Explore models with one phylogenetic tree

##Transition models in corHMM : 

##Models without hidden states: 

##model ER: 
fit.marginal.ER <- corHMM(phylo1, data,
                     rate.cat=1, 
                     model = "ER",
                     node.states = "marginal", 
                     ip = 1, ##DEFAULT
                     nstarts = 1000,
                     n.cores = 2,
                     get.tip.states = TRUE,
                     lewis.asc.bias = FALSE) ##DEFAULT

plotMKmodel(fit.marginal.ER)

##Model SYM: 
fit.marginal.SYM <- corHMM(phylo1, data,
                          rate.cat=1, 
                          model = "SYM",
                          node.states = "marginal", 
                          ip = 1, ##DEFAULT
                          nstarts = 1000,
                          n.cores = 2,
                          get.tip.states = TRUE,
                          lewis.asc.bias = FALSE) 

plotMKmodel(fit.marginal.SYM)

##Model ARD:
fit.marginal.ARD <- corHMM(phylo1, data,
                           rate.cat=1, 
                           model = "ARD",
                           node.states = "marginal", 
                           ip = 1, ##DEFAULT
                           nstarts = 1000,
                           n.cores = 2,
                           get.tip.states = TRUE,
                           lewis.asc.bias = FALSE)

plotMKmodel(fit.marginal.ARD)

fit.marginal.ARD$index.mat

##Model with 2 hidden states (ER/ER): 

fit.marginal.ER.2 <- corHMM(phylo1, data,
                          rate.cat=2, 
                          model = "ER",
                          node.states = "marginal", 
                          ip = 1, ##DEFAULT
                          nstarts = 1000,
                          n.cores = 2,
                          get.tip.states = TRUE,
                          lewis.asc.bias = FALSE)


plotMKmodel(fit.marginal.ER.2)

##Model with 2 hidden states (SYM/SYM):

fit.marginal.SYM.2 <- corHMM(phylo1, data,
                           rate.cat=2, 
                           model = "SYM",
                           node.states = "marginal", 
                           ip = 1, ##DEFAULT
                           nstarts = 1000,
                           n.cores = 2,
                           get.tip.states = TRUE,
                           lewis.asc.bias = FALSE)

plotMKmodel(fit.marginal.SYM.2)

##Model with 2 hidden states(ARD/ARD):

fit.marginal.ARD.2 <- corHMM(phylo1, data,
                           rate.cat=2, 
                           model = "ARD",
                           node.states = "marginal", 
                           ip = 1, ##DEFAULT
                           nstarts = 1000,
                           n.cores = 2,
                           get.tip.states = TRUE,
                           lewis.asc.bias = FALSE) 

plotMKmodel(fit.marginal.ARD.2)

## Model with 2 hidden states and with mix models (ER/ARD; ER/SYM; ARD/SYM):

##Make a custom rate matrices: 

##Matrix ER:
matrix.ER <- getStateMat4Dat(data, model = "ER", dual = FALSE, collapse = TRUE, indep = FALSE)
rate.mat.ER <- matrix.ER$rate.mat
rate.mat.ER

##This is the same than: 
fit.marginal.ER$index.mat

rate.mat.ER <- fit.marginal.ER$index.mat

##Matrix ARD: 
matrix.ARD <- getStateMat4Dat(data, model = "ARD", dual = FALSE, collapse = TRUE, indep = FALSE)
rate.mat.ARD <- matrix.ARD$rate.mat
rate.mat.ARD
##This is the same than:
fit.marginal.ARD$index.mat
 
rate.mat.ARD <- fit.marginal.ARD$index.mat

##Matrix SYM:
matrix.SYM <- getStateMat4Dat(data, model = "SYM", dual = FALSE, collapse = TRUE, indep = FALSE)
rate.mat.SYM <- matrix.SYM$rate.mat
rate.mat.SYM
##This is the same than:
fit.marginal.SYM$index.mat
 
rate.mat.SYM <- fit.marginal.SYM$index.mat


#transition rate from R1 to R2:
 
RateClassMat <- getRateCatMat(2) #
RateClassMat <- equateStateMatPars(RateClassMat, c(1,2))
RateClassMat

##corHMM matrix:   

#Model ER/ARD:
StateMats.ER.ARD <- list(rate.mat.ER, rate.mat.ARD)
StateMats.ER.ARD

FullMat.ER.ARD <- getFullMat(StateMats.ER.ARD, RateClassMat)
FullMat.ER.ARD

plotMKmodel(FullMat.ER.ARD, rate.cat = 2, display = "row", text.scale = 0.7)

#Model ER/SYM:
StateMats.ER.SYM <- list(rate.mat.ER, rate.mat.SYM)
StateMats.ER.SYM

FullMat.ER.SYM <- getFullMat(StateMats.ER.SYM, RateClassMat)
FullMat.ER.SYM

plotMKmodel(FullMat.ER.SYM, rate.cat = 2, display = "row", text.scale = 0.7)

#Model ARD/SYM:
StateMats.ARD.SYM <- list(rate.mat.ARD, rate.mat.SYM)
StateMats.ARD.SYM

FullMat.ARD.SYM <- getFullMat(StateMats.ARD.SYM, RateClassMat)
FullMat.ARD.SYM

plotMKmodel(StateMats.ARD.SYM, rate.cat = 2, display = "row", text.scale = 0.7)

##Make models:

#Model ER/ARD:   

fit.marginal.ER.ARD <- corHMM(phylo1, data,
                             rate.cat=2, 
                             rate.mat = FullMat.ER.ARD,
                             node.states = "marginal", 
                             ip = 1, ##DEFAULT
                             nstarts = 1000,
                             n.cores = 2,
                             get.tip.states = TRUE,
                             lewis.asc.bias = FALSE) 

plotMKmodel(fit.marginal.ER.ARD)

#Model ER/SYM:   

fit.marginal.ER.SYM <- corHMM(phylo1, data,
                              rate.cat=2, 
                              rate.mat = FullMat.ER.SYM,
                              node.states = "marginal", 
                              ip = 1, ##DEFAULT
                              nstarts = 1000,
                              n.cores = 2,
                              get.tip.states = TRUE,
                              lewis.asc.bias = FALSE) 

plotMKmodel(fit.marginal.ER.SYM)

#Model ARD/SYM:  

fit.marginal.ARD.SYM <- corHMM(phylo1, data,
                              rate.cat=2, 
                              rate.mat = FullMat.ARD.SYM,
                              node.states = "marginal", 
                              ip = 1, ##DEFAULT
                              nstarts = 1000,
                              n.cores = 2,
                              get.tip.states = TRUE,
                              lewis.asc.bias = FALSE) 

plotMKmodel(fit.marginal.ARD.SYM)


##Compare AICc among models: 

model.names <- c("ER", "SYM", "ARD", "ER.ER", "SYM.SYM", "ARD.ARD", "ER.ARD", "ER.SYM", "ARD.SYM")

AICc <- c(fit.marginal.ER$AICc, fit.marginal.SYM$AICc, fit.marginal.ARD$AICc, fit.marginal.ER.2$AICc, fit.marginal.SYM.2$AICc, fit.marginal.ARD.2$AICc, fit.marginal.ER.ARD$AICc, fit.marginal.ER.SYM$AICc, fit.marginal.ARD.SYM$AICc)

AICc <- rbind(model.names, AICc)
AICc

##The best model is the ER model without hidden states:
 
dev.off()
plotMKmodel(fit.marginal.ER, text.scale = 0.9, display = "square", color = 2)

####################################################
#Part 2: Make the same analysis to 200 phylogenetic tree:


AICc.vf <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(AICc.vf) <- c("ER", "SYM", "ARD", "ER.ER", "SYM.SYM", "ARD.ARD", "ER.ARD", "ER.SYM", "ARD.SYM")

trees <- c(1:200)

for (i in trees) {
  
  
  ##model ER: 
  fit.marginal.ER <- corHMM(phylo2[[i]], data,
                            rate.cat=1, 
                            model = "ER",
                            node.states = "marginal", 
                            ip = 1, ##DEFAULT
                            nstarts = 1000,
                            n.cores = 2,
                            get.tip.states = TRUE,
                            lewis.asc.bias = FALSE) ##DEFAULT
  
  #plotMKmodel(fit.marginal.ER)
  
  ##Model SYM: 
  fit.marginal.SYM <- corHMM(phylo2[[i]], data,
                             rate.cat=1, 
                             model = "SYM",
                             node.states = "marginal", 
                             ip = 1, ##DEFAULT
                             nstarts = 1000,
                             n.cores = 2,
                             get.tip.states = TRUE,
                             lewis.asc.bias = FALSE) 
  
  #plotMKmodel(fit.marginal.SYM)
  
  ##Model ARD:
  fit.marginal.ARD <- corHMM(phylo2[[i]], data,
                             rate.cat=1, 
                             model = "ARD",
                             node.states = "marginal", 
                             ip = 1, ##DEFAULT
                             nstarts = 1000,
                             n.cores = 2,
                             get.tip.states = TRUE,
                             lewis.asc.bias = FALSE)
  
  #plotMKmodel(fit.marginal.ARD)
  
  ##Model (ER/ER): 
  
  fit.marginal.ER.2 <- corHMM(phylo2[[i]], data,
                              rate.cat=2, 
                              model = "ER",
                              node.states = "marginal", 
                              ip = 1, ##DEFAULT
                              nstarts = 1000,
                              n.cores = 2,
                              get.tip.states = TRUE,
                              lewis.asc.bias = FALSE)
  
  
  #plotMKmodel(fit.marginal.ER.2)
  
  ##Model (SYM/SYM):
  
  fit.marginal.SYM.2 <- corHMM(phylo2[[i]], data,
                               rate.cat=2, 
                               model = "SYM",
                               node.states = "marginal", 
                               ip = 1, ##DEFAULT
                               nstarts = 1000,
                               n.cores = 2,
                               get.tip.states = TRUE,
                               lewis.asc.bias = FALSE)
  
  
  ##Model (ARD/ARD):
  
  fit.marginal.ARD.2 <- corHMM(phylo2[[i]], data,
                               rate.cat=2, 
                               model = "ARD",
                               node.states = "marginal", 
                               ip = 1, ##DEFAULT
                               nstarts = 1000,
                               n.cores = 2,
                               get.tip.states = TRUE,
                               lewis.asc.bias = FALSE) 
  
  #plotMKmodel(fit.marginal.ARD.2)
  
  
  ##Matrix ER:
  ##matrix.ER <- getStateMat4Dat(data, model = "ER", dual = FALSE, collapse = TRUE, indep = FALSE)
  #rate.mat.ER <- matrix.ER$rate.mat
  #rate.mat.ER
  ##Essa matrix é a mesma obtida em: 
  #fit.marginal.ER$index.mat
  ##Portanto vamos usar ela como rate mat: 
  rate.mat.ER <- fit.marginal.ER$index.mat
  
  ##Matrix ARD: 
  #matrix.ARD <- getStateMat4Dat(data, model = "ARD", dual = FALSE, collapse = TRUE, indep = FALSE)
  #rate.mat.ARD <- matrix.ARD$rate.mat
  #rate.mat.ARD
  ##Essa matrix é a mesma obtida em: 
  #fit.marginal.ARD$index.mat
  ##Portanto vamos usar ela como rate mat: 
  rate.mat.ARD <- fit.marginal.ARD$index.mat
  
  ##Matrix SYM:
  #matrix.SYM <- getStateMat4Dat(data, model = "SYM", dual = FALSE, collapse = TRUE, indep = FALSE)
  #rate.mat.SYM <- matrix.SYM$rate.mat
  #rate.mat.SYM
  ##Essa matrix é a mesma obtida em: 
  #fit.marginal.SYM$index.mat
  ##Portanto vamos usar ela como rate mat: 
  rate.mat.SYM <- fit.marginal.SYM$index.mat
  
##Matrix Rate: 
  
  RateClassMat <- getRateCatMat(2) #
  RateClassMat <- equateStateMatPars(RateClassMat, c(1,2))
  RateClassMat
  
  
  #Matrix ER/ARD:
  StateMats.ER.ARD <- list(rate.mat.ER, rate.mat.ARD)
  StateMats.ER.ARD
  
  FullMat.ER.ARD <- getFullMat(StateMats.ER.ARD, RateClassMat)
  FullMat.ER.ARD
  
  #plotMKmodel(FullMat.ER.ARD, rate.cat = 2, display = "row", text.scale = 0.7)
  
  #Matrix ER/SYM:
  StateMats.ER.SYM <- list(rate.mat.ER, rate.mat.SYM)
  StateMats.ER.SYM
  
  FullMat.ER.SYM <- getFullMat(StateMats.ER.SYM, RateClassMat)
  FullMat.ER.SYM
  
  #plotMKmodel(FullMat.ER.SYM, rate.cat = 2, display = "row", text.scale = 0.7)
  
  #Matrix ARD/SYM:
  StateMats.ARD.SYM <- list(rate.mat.ARD, rate.mat.SYM)
  StateMats.ARD.SYM
  
  FullMat.ARD.SYM <- getFullMat(StateMats.ARD.SYM, RateClassMat)
  FullMat.ARD.SYM
  
  #plotMKmodel(StateMats.ARD.SYM, rate.cat = 2, display = "row", text.scale = 0.7)
  
  
  #Matrix ER/ARD:   
  
  fit.marginal.ER.ARD <- corHMM(phylo2[[i]], data,
                                rate.cat=2, 
                                rate.mat = FullMat.ER.ARD,
                                node.states = "marginal", 
                                ip = 1, ##DEFAULT
                                nstarts = 1000,
                                n.cores = 2,
                                get.tip.states = TRUE,
                                lewis.asc.bias = FALSE) 
  
  #plotMKmodel(fit.marginal.ER.ARD)
  
  #Model ER/SYM:   
  
  fit.marginal.ER.SYM <- corHMM(phylo2[[i]], data,
                                rate.cat=2, 
                                rate.mat = FullMat.ER.SYM,
                                node.states = "marginal", 
                                ip = 1, ##DEFAULT
                                nstarts = 1000,
                                n.cores = 2,
                                get.tip.states = TRUE,
                                lewis.asc.bias = FALSE) 
  
  #plotMKmodel(fit.marginal.ER.SYM)
  
  #Model ARD/SYM:  
  
  fit.marginal.ARD.SYM <- corHMM(phylo2[[i]], data,
                                 rate.cat=2, 
                                 rate.mat = FullMat.ARD.SYM,
                                 node.states = "marginal", 
                                 ip = 1, ##DEFAULT
                                 nstarts = 1000,
                                 n.cores = 2,
                                 get.tip.states = TRUE,
                                 lewis.asc.bias = FALSE) 
  
  #plotMKmodel(fit.marginal.ARD.SYM)
  
  
  ##Compare AICc among models: 
  
  AICc <- c(fit.marginal.ER$AICc, fit.marginal.SYM$AICc, fit.marginal.ARD$AICc, fit.marginal.ER.2$AICc, fit.marginal.SYM.2$AICc, fit.marginal.ARD.2$AICc, fit.marginal.ER.ARD$AICc, fit.marginal.ER.SYM$AICc, fit.marginal.ARD.SYM$AICc)
  
AICc.vf <- rbind(AICc.vf, AICc)
  
}

write.csv(AICc.vf, paste0(path, "Results/03-Life_form_ancestral_reconstruction-SIMMAP/corHMM_AICc_models.csv"))

####################################################
##Part 3: Make SIMMAP with the best Model (ER model):

##In corHMM analysis the ER Model without hidden states had the Best Fit. We use this model to make Simmap: 

##life form classification: 

data <- read.csv(paste0(path, "/data-paper/data/life_form.csv", header = T))

life.form <- data[,-1]

names(life.form) <- data$specie

life.form

SIMMAP.ER <- make.simmap(phylo2, life.form, model="ER", nsim=10)

SIMMAP.ER.SUMMARY <- describe.simmap(SIMMAP.ER)

plot(SIMMAP.ER.SUMMARY)
nodelabels(cex = 0.5, font = 0.4)

##Save SIMMAP and summary Simmap in .rds file:

##set working directory to the data-paper folder:

saveRDS(SIMMAP.ER, "data-paper\Results\03-Life_form_ancestral_reconstruction-SIMMAP/simmap_ER_model.rds")

saveRDS(SIMMAP.ER.SUMMARY, "data-paper\Results\03-Life_form_ancestral_reconstruction-SIMMAP/simmap_ER_model_summary.rds")

##Data explore: 

##Transition number: 

count <- countSimmap(SIMMAP.ER, states=NULL, message=TRUE)

t <- as.data.frame(count$Tr)

write.csv(t, "transition_life_form.csv")

medias <- colMeans(t)
sort(medias, decreasing = FALSE)

medianas <- apply(t,2,median)
sort(medianas, decreasing = FALSE)

sd <- apply(t,2,sd)
sort(sd, decreasing = FALSE)

##posterior probabilities of each node being in each state:

ace <- SIMMAP.ER.SUMMARY$ace

write.csv(ace, "ace_lifeform.csv")

##End
