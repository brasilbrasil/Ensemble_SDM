###################################
####Setting data directories#######
###################################

rootDir = "D:" #for server rootDir = "Y:/"
dataDir = paste0(rootDir, "/PICCC_data") #where the original data exists for SDM (PA species and environmental data) - same for server
codeDir = paste0(rootDir, "/Dropbox/code/Ensemble_SDM")  
#codeDir = paste0(rootDir, "/Dropbox/code/Ensemble_SDM/IRC")  
analysisDir = paste0(rootDir, "/PICCC_analysis") #same for server - assigns location for analysis directory (results)
DR_code_S=paste0(rootDir, "/Dropbox/code/")
DR_PICCC_data_S=paste0(rootDir, "/PICCC_data")
resultsDir = paste0(analysisDir,"/FB_analysis/model_results/biomod2")

#checks if results directory exists and creates it if not
if (file.exists(resultsDir) == FALSE){
  dir.create(resultsDir)
}

#Location for species data
allSppNames = paste0(analysisDir, "/spp_to_run_all.csv")

#Location for bioclim data
bioclimDataDir = paste0(dataDir,"/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg_rounded/bioclims_abs/") #points to location of bioclimate data
#bioclimData2013Dir = paste0(dataDir,"/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg_rounded/bioclims_abs/")
bioclimData2013Dir = paste0(dataDir,"/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg_rounded/seasonal data/")

fitting_clim_data_dir = paste0(bioclimData2013Dir,"all_baseline/125m/") 
necessary_run_data = paste0(resultsDir,"/necessary_run_data/") #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)

clim_data_2000=paste0(bioclimData2013Dir,"all_baseline/500m/")
clim_data_2100=paste0(bioclimData2013Dir,"all_future/500m/")
clim_data_2000=paste0(bioclimData2013Dir,"all_baseline/500m/")
clim_data_2100=paste0(bioclimData2013Dir,"all_future/500m/")
clim_data_2000wettest="D:/GIS_Data/REnviroLayers/mixed_data_2000_250mwettest/"
clim_data_2000driest= "D:/GIS_Data/REnviroLayers/mixed_data_2000_250mdriest/"
clim_data_2100wettest="D:/GIS_Data/REnviroLayers/mixed_data_2100_250mwettest/"
clim_data_2100driest= "D:/GIS_Data/REnviroLayers/mixed_data_2100_250mdriest/"  

#original DD data
DR_3km_DD_data_S = paste0(dataDir,"/climate_data/bioclim_recalc/2- DD data")