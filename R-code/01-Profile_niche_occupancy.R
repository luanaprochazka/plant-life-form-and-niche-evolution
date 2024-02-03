# ------------------------------------------
# Calc the Predict Niche Occupancy - PNO 
# (See Evans, 2009)
# ------------------------------------------

# This is part of:

#Resource Availability and Disturbance Frequency Shape Plant Life Forms in Neotropical Habitats 

#Authors: Luana S. Prochazka, Suzana Alcantara, Juliana Gastaldello Rando, Thais Vasconcelos, Raquel C. Pizzardo and Anselmo Nogueira

#Journal: New Phytologist
#Article acceptance date: 30 January 2024


#uncomment to clean your R environment:

#rm(list = ls())

##packages: 

library(raster)
library(phyloclim)
library(tibble)
library(tidyverse)
library(dplyr)
library(sp)
library(maptools)
library(ape)
library(dismo)


##See pno function code: 
 
pno

##Copy pno function code and bug correction - In line 43 "nrows|ncols" was raplace by "NROWS|NCOLS": 

pno2 <- function (path_bioclim, path_model, subset = NULL, bin_width = 1, 
                    bin_number = NULL) 
{
  hdr <- read.table(path_bioclim, row.names = 1, nrows = 6)
  bc.cells.dim <- hdr[grep("NROWS|NCOLS", rownames(hdr)), 
  ]
  nodata <- hdr[rownames(hdr) == "NODATA_value", ]
  bc <- scan(path_bioclim, what = "i", skip = 6, quiet = TRUE)
  bc <- as.numeric(bc)
  bc.na <- which(bc == nodata)
  bc[bc.na] <- NA
  bc_min <- min(bc, na.rm = TRUE)
  bc_max <- max(bc, na.rm = TRUE)
  if (is.null(bin_number)) {
    by <- bin_width
  }
  else {
    by <- ((bc_max - bc_min)/(bin_number - 1))
  }
  
  cats <- seq(from = bc_min, to = bc_max, by = by)
  
  cat("\nBioclimatic data ranges from", bc_min, "to", 
      bc_max, "and will be binned into", length(cats), 
      "categories with a bin width of", by, "\n")
  
  
  specs <- list.files(path = path_model, pattern = ".asc", 
                      full.names = TRUE)
  if (length(specs) == 0) 
    stop("\nNo models found. Please check your path_model ", 
         "argument!")
  clamping <- grep("clamping", specs)
  if (length(clamping) > 0) 
    specs <- specs[-clamping]
  if (!is.null(subset)) 
    specs <- specs[grep(paste(subset, collapse = "|"), 
                        specs)]
  label <- gsub(path_model, "", specs)
  label <- gsub("/|.asc", "", label)
  
  x <- matrix(nrow = length(cats), ncol = length(specs) + 1)
  colnames(x) <- c("variable", label)
  x[, 1] <- cats
  
  cum.prob.bin <- function(x, bc, enm, normalizer, bin) {
    id <- which(bc > x - bin/2 & bc <= x + bin/2)
    sum(enm[id], na.rm = TRUE)/normalizer
  }
  
  
  for (h in seq(along = specs)) {
    cat("\n\nReading ENM for", specs[h], "...")
    
    hdr <- read.table(specs[h], row.names = 1, nrows = 6)
    
    enm.cells.dim <- hdr[grep("nrows|ncols", rownames(hdr)), 
    ]
    
    nodata <- hdr[rownames(hdr) == "NODATA_value", 
    ]
    
    enm <- scan(specs[h], what = "i", skip = 6, quiet = TRUE)
    enm <- as.numeric(enm)
    enm.na <- which(enm == nodata)
    enm[enm.na] <- NA
    cat("finished.")
    
    if (!identical(bc.cells.dim, enm.cells.dim)) 
      stop("Resolution and/or extent of ENM and", 
           " bioclimatic layer do not match!")
    px <- abs(length(bc.na) - length(enm.na))
    if (px > 0) 
      cat("\n\nWARNING: Raster maps differ by", px, 
          "cells", "(perhaps due to differing coastlines.)", 
          "\nThese cells are treated as having zero", 
          "probability in the MAXENT distribution.")
    normalizer <- sum(enm, na.rm = TRUE)
    cat("\n\nBinning probabilities for", specs[h], 
        "...")
    prob <- sapply(cats, cum.prob.bin, bc = bc, enm = enm, 
                   normalizer = normalizer, bin = by)
    cat("finished.")
    remove(enm)
    x[, h + 1] <- prob
  }
  if (length(grep("Temperature|Diurnal", path_bioclim)) == 
      1) 
    x[, 1] <- x[, 1]/10
  if (length(grep("Temperature_Seasonality|Isothermality", 
                  path_bioclim)) == 1) 
    x[, 1] <- x[, 1]/100
  x
}


## set working directory with the environmental layers (.asc files):  

env <- "/environmental_layers/"

## set working directory with the maxent models results (.asc files): 

path_model <- "/maxent_models_results/"

##Make a list with each environmental layer at environmental_layers folder:

setwd(env)
list <- dir(pattern = ".asc")

##set working directory to save .csv pno files results:  

setwd("/data-paper/Results/01-pno")

##Calc pno for all environmental layers in env directory: 


for (i in list) {
  
  path_bioclim <- paste0(env, i)
  
  pno <- pno2(path_bioclim = path_bioclim,
              path_model = path_model, 
              subset = NULL, 
              bin_width = 1, bin_number = NULL)
  
  pno <- as.data.frame(pno)
  
  write.csv(pno, file = paste0("pno_", i, ".csv"), row.names = FALSE)
  
}

##Calc pno for fire layer with bin_width parameter = 0.01: 


setwd("/data-paper/Results/01-pno")

pno <- pno2(path_bioclim = path_bioclim,
            path_model = path_model, 
            subset = NULL, 
            bin_width = 0.01, bin_number = NULL)

pno <- as.data.frame(pno)

write.csv(pno, file = "pno_fire_frequency_0.01.csv", row.names = FALSE)

##end