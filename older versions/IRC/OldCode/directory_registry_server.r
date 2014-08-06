###################################
####Setting data directories#######
###################################

#PICCC data directory
rootDir = "Y:"
dataDir = paste0(rootDir, "/PICCC_data") #where the original data exists for SDM (PA species and environmental data) - same for server
codeDir = paste0(rootDir, "/PICCC_analysis/code/IRC/Ensemble_SDM")
analysisDir = paste0(rootDir, "/PICCC_analysis")
resultsDir = paste0(analysisDir, "/FB_analysis/model_results/biomod2" )

#Location for species data
allSppNames = paste0(resultsDir, "/spp_to_run_all.csv")

#Location for climate data
bioclimDataDir=paste0(dataDir,"/climate_data/full extent bioclim data") #points to location of bioclimate data
clim_data_2000=paste0(bioclimDataDir,"all_grd/all_baseline/250m/")
clim_data_2100=paste0(bioclimDataDir,"all_grd/all_future/500m/")

fitting_clim_data_dir=paste0(bioclimDataDir,"/all_baseline/100m") #where climate data i stored
necessary_run_data=paste0(resultsDir,'/necessary_run_data') #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)

