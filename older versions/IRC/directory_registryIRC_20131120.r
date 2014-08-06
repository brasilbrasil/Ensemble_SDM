###################################
####Setting data directories#######
###################################

rootDir = "C:/USGS_Honolulu/" #for server rootDir = "Y:/"
dataDir = paste0(rootDir, "PICCC_data/") #where the original data exists for SDM (PA species and environmental data) - same for server
codeDir = paste0(rootDir, "PICCC_code/Ensemble_SDM/IRC/")  

analysisDir = paste0(rootDir, "PICCC_analysis/biomod2/") #same for server - assigns location for analysis directory (results)
resultsDir = paste0(analysisDir,"results/")

dir.create(resultsDir, showWarnings = F) #creates results directory if missing

#Location for species data
allSppNames = paste0(analysisDir, "spp_to_run_all.csv")

#Location for bioclim data
bioclimDataDir = paste0(dataDir,"bioclimData/") #points to location of bioclimate data
clim_data_2000 = paste0(bioclimDataDir,"all_grd/all_baseline/250m/")
clim_data_2100 = paste0(bioclimDataDir,"all_grd/all_future/500m/")
fitting_clim_data_dir = paste0(bioclimDataDir,"all_grd/all_baseline/100m") 
necessary_run_data = paste0(analysisDir,"necessary_run_data/") #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
