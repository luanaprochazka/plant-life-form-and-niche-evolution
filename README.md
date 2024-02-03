# plant-life-form-and-niche-evolution:
This repository contain all our data and R scripts for the manuscript:

#Resource Availability and Disturbance Frequency Shape Plant Life Forms in Neotropical Habitats 

#Journal: New Phytologist

#Article acceptance date: 30 January 2024

##A brief summary of what the study is about:

In this study we investigate whether the occupation of different environmental niches occurred jointly with the evolution of plant life forms. For this we used trait-dependent evolution models, including pure Brownian motion (BM) and more complex Ornstein–Uhlenbeck (OU) models. Here we provide the R code that we use to run our analyzes and the main data generated.
For more details you can see the publication associated with this datasets. 

##Who is responsible for data and writing code:

The corresponding author was responsible for writing code.
For any issues, contact corresponding author: Luana S. Prochazka, prochazka.luana@gmail.com

FOLDERS AND FILE LIST: 

This dataset is organized in three main folders: 

*data folder: contain the primary files data:

	* life_form.csv: A .csv file with life form classification for each specie.
	* species_occurence.xlsx: A excel file with species occurrence data. The metadata for this file is available in the "metadata" tab within the file. 
	* tree_vasconcelos41spp.txt:  A phylogenetic tree file with 41 terminals. This is a prunned tree. The original tree is available in: Vasconcelos, Thais et al. (2020), Fast diversification through a mosaic of evolutionary histories characterizes the endemic flora of ancient Neotropical mountains, Dryad, Dataset, https://doi.org/10.5061/dryad.x69p8czf6
	* tree200_vasconcelos41spp.txt: 200 phylogenetic tree samples with 41 terminals. This is a prunned tree. The original tree is available in: Vasconcelos, Thais et al. (2020), Fast diversification through a mosaic of evolutionary histories characterizes the endemic flora of ancient Neotropical mountains, Dryad, Dataset, https://doi.org/10.5061/dryad.x69p8czf6

*R-code folder: contain 9 R-scripts that we use to run our analyzes. Details about each file are available within the README file available within the R-code folder:

	* 01-Profile_niche_occupancy.R
	* 02-Niche-similarity_and_age_correlation.R
	* 03-Life_form_ancestral_reconstruction-corHMM.R
	* 04-OUwie_models.R
	* 05-WN_models.R
	* 06A-Bio17-model_selection.R
	* 06B-Bio18-model_selection.R
	* 06C-fire-model_selection.R
	* 06D-nitro_model_selection.R

	* README.txt: Contain details about .R files and instructions to run R code. 

*Results folders: Contains 4 folders. The  intern folders contains  files that are generated by R scripts. Many of them are in .rds or .csv files and are the way to build the final results. You can find more information in README file at the R-code folder. 

	*01-pno folder: Contains four .csv files with a profile niche occupancy (pno files) for each environmental layer. This files are generated by 01-Profile_niche_occupancy.R script. In these files, the first column corresponds to environmental layer values. Each subsequent column corresponds to the probability of each species occupying the environmental value described in the first column. The sum of the probabilities of each species is 1. These files are: 
 
		*pno_bio17.csv. 
		*pno_bio18.csv. 
		*pno_fire_frequency_0.01.csv. 
		*pno_nitrogen_0_30cm.csv. 

*02-Niche_overlap_and_age_correlation folder: contain .csv files with niche overlap results (D and I metrics) in a species matrix. The .rds files contains the age range correlation result. All files in this folder are generated by 02-Niche-similarity_and_age_correlation.R script.These files are:

		*geo_overlap.csv                                  
		*overlap_pno_bio17.asc.csv                        
		*overlap_pno_bio18.asc.csv                        
		*overlap_pno_fire_frequency_0.01.csv              
		*overlap_pno_nitrogen_0_30cm.asc.csv
		*age_correlation_lower_bio17.rds                 
		*age_correlation_lower_bio18.rds                  
		*age_correlation_lower_fire_frequency_0.01.csv.rds
		*age_correlation_lower_geo.rds                    
		*age_correlation_lower_nitrogen_0_30cm.rds        
		*age_correlation_upper_bio17.rds                  
 		*age_correlation_upper_bio18.rds                  
		*age_correlation_upper_fire_frequency_0.01.csv.rds
		*age_correlation_upper_geo.rds                    
    *age_correlation_upper_nitrogen_0_30cm.rds        

*03-Life_form_ancestral_reconstruction-SIMMAP folder: Contain .csv file with AICc results to corHMM analyses and the .RDS files generated by 03-Life_form_ancestral_reconstruction-corHMM.R script. These files contain 2000 SIMMAP simulations and the summary results. These files are:

		*corHMM_AICc_models.csv: contain the AICC values to corHMM models.
		*simmap_ER_model.rds: Contain the 2000 SIMMAP simulation
		*simmap_ER_model_summary.rds: Contains the summary results of 2000 SIMMAP simulation.
		
*04-OUwie_and_WN_models folder: Contains folders for each environmental variable and inside them are the fit models results (.rds files) and the .csv files with best fit model selection (Aicc and deltaAicc tables). All the files in these folder are generated by scripts: 04, 05, 6A, 6B, 6C and 6D. For more information about the files please see the R code. These folders and files are:

		*bio17 folder: contain the results file about bio17 environmental axis
			*bio17_500samples_pno.csv
			*bio17_AICc.csv          
			*bio17_AICc_all.csv
      *bio17_best_fit_all.csv
   		*bio17_Delta_AICc_all.csv
 			*bio17_fitBM1.rds  
			*bio17_fitBMS.rds 
			*bio17_fitOU1.rds        
			*bio17_fitOUM.rds    
			*bio17_fitOUMA.rds    
			*bio17_fitOUMV.rds       
			*bio17_fitOUMVA.rds    
			*bio17_fitWN.rds     
			*bio17log_data.csv       

    *bio18 folder: contain the results file about bio18 environmental axis

			*bio18_500samples_pno.csv
			*bio18_AICc.csv          
			*bio18_AICc_all.csv
      *bio18_best_fit_all.csv
   		*bio18_Delta_AICc_all.csv
 			*bio18_fitBM1.rds  
			*bio18_fitBMS.rds 
			*bio18_fitOU1.rds        
			*bio18_fitOUM.rds    
			*bio18_fitOUMA.rds    
			*bio18_fitOUMV.rds       
			*bio18_fitOUMVA.rds    
			*bio18_fitWN.rds     
			*bio18log_data.csv   
		
		*fire folder: contain the results file about fire environmental axis

			*fire_500samples_pno.csv
			*fire_AICc.csv          
			*fire_AICc_all.csv
      *fire_best_fit_all.csv
   		*fire_Delta_AICc_all.csv
 			*fire_fitBM1.rds  
			*fire_fitBMS.rds 
			*fire_fitOU1.rds        
			*fire_fitOUM.rds    
			*fire_fitOUMA.rds    
			*fire_fitOUMV.rds       
			*fire_fitOUMVA.rds    
			*fire_fitWN.rds     
			*firelog_data.csv 

		*nitro folder: contain the results file about nitro environmental axis
 
			*nitro_500samples_pno.csv
			*nitro_AICc.csv          
			*nitro_AICc_all.csv
      *nitro_best_fit_all.csv
   		*nitro_Delta_AICc_all.csv
 			*nitro_fitBM1.rds  
			*nitro_fitBMS.rds 
			*nitro_fitOU1.rds        
			*nitro_fitOUM.rds    
			*nitro_fitOUMA.rds    
			*nitro_fitOUMV.rds       
			*nitro_fitOUMVA.rds    
			*nitro_fitWN.rds     
			*nitrolog_data.csv 
