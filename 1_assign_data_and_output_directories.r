###################################
#FB specific directories
###################################

#set analysis output folder
codeDir = paste0(DR_code_S, "Ensemble_SDM/")  
resultsDir = paste0(analysisDir,"FB_analysis/model_results/biomod2/")
necessary_run_data = paste0(resultsDir,"necessary_run_data/") #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
dir.create(resultsDir, showWarnings=F)
#Location for species data
allSppNames = paste0(analysisDir, "spp_to_run_all.csv")
#Location for bioclim predictors
bioclimDataDir = paste0(bioClimDir,"allYrs_avg_rounded/bioclims_abs/") #points to location of bioclimate data
#bioclimDataDir = paste0(bioClimDir,"allYrs_avg_rounded/seasonal data/")
fitting_clim_data_dir = paste0(bioclimDataDir,"all_baseline/125m/") 
#clim_data_2000=paste0(bioclimDataDir,"all_baseline/500m/")
clim_data_2000=paste0(bioclimDataDir,"all_baseline/500m/")
#clim_data_2100=paste0(bioclimDataDir,"all_future/500m/")
clim_data_2100=paste0(bioclimDataDir,"all_future/500m/")
clim_data_2100rev="D:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg_revDelta_rounded/bioclims_abs/all_future/500m/"
clim_data_2000wettest="D:/GIS_Data/REnviroLayers/mixed_data_2000_250mwettest/"
clim_data_2000driest= "D:/GIS_Data/REnviroLayers/mixed_data_2000_250mdriest/"
clim_data_2100wettest="D:/GIS_Data/REnviroLayers/mixed_data_2100_250mwettest/"
clim_data_2100driest= "D:/GIS_Data/REnviroLayers/mixed_data_2100_250mdriest/"  
#original DD data
DR_3km_DD_data_S = paste0(dataDir,"climate_data/bioclim_recalc/2- DD data")