rm(list = ls()) #remove all past worksheet variables
source(paste0("Y:/PICCC_analysis/code/","directory_registry.r"))

###################################
####GENERAL MODEL CONFIGURATION####
###################################
local_config_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m/') #'C:/Users/lfortini/'
#spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
spp_nm=c("Oahu_Amakihi")#,"Akekee", "Anianiau", "Kauai_Amakihi", "Oahu_Amakihi","Hawaii_Akepa", "Hawaii_Elepaio", "Palila")
server=1
overwrite=1
models_to_run=c('GBM','MAXENT')
eval_stats=c("ROC")
plot_graphs=1
EM_fit=T
EM_ensemble=F
EM_project=F
memory.limit(size=24000000)

#source(paste0(DR_code_S,"Ensemble_SDM/0_SDM_run_config.r")) #this is where all configurations are at

if (server==1){
  working_dir=paste0(DR_FB_SDM_results_S,'test2_runs_500m_rounded/')
  fitting_clim_data_dir=paste0(DR_FB_clim_data,"all_grd/all_baseline/500m_test/") 
  necessary_run_data=paste0(DR_FB_SDM_results_S,'necessary_run_data/') #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/test/'
  necessary_run_data='C:/Users/lfortini/Data/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
  fitting_clim_data_dir="C:/Users/lfortini/Data/SDM_env_data/all_grd/all_baseline/100m/"
}

env_var_files=c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") 
crop_raster_dir=paste(working_dir, 'map_crop/',sep="")
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")

#################################
####CONFIG FOR SPECIFIC STEPS####
#################################
####fit config
remove_PA_abs=TRUE
NbRunEval=10

####ensemble config

####projection config
baseline_or_future=1 #0 for baseline, 1 for future
memory = T #keep.in.memory=memory
temp<-paste('Y:/temp/', baseline_or_future, '/', sep='') #dir for temp run data (to avoid memory errors)

if (server==1){
  clim_data_2000=paste0(DR_FB_clim_data,"all_grd/all_baseline/500m_test/")
  clim_data_2100=paste0(DR_FB_clim_data,"all_grd/all_future/500m/")
  clim_data_2000wettest="D:/GIS_Data/REnviroLayers/mixed_data_2000_250mwettest/"
  clim_data_2000driest= "D:/GIS_Data/REnviroLayers/mixed_data_2000_250mdriest/"
  clim_data_2100wettest="D:/GIS_Data/REnviroLayers/mixed_data_2100_250mwettest/"
  clim_data_2100driest= "D:/GIS_Data/REnviroLayers/mixed_data_2100_250mdriest/"  
}else{
  clim_data_2000="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/250m/"
  clim_data_2100="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_future/250m/"
}

##########################
####RUNNING SCRIPTS!!!####
##########################
maxentWDtmp = paste("maxentWDtmp_", baseline_or_future, sep = "")
dir.create(temp)

if (EM_fit){
  source(paste0(DR_code_S,"Ensemble_SDM/1_BM2_FB_SDM_fitting_w_cropping4.r")) #this is where all configurations are at
}
if (EM_ensemble){
  source(paste0(DR_code_S,"Ensemble_SDM/2_BM2_FB_SDM_EM_fitting2.r")) #this is where all configurations are at
}
if (EM_project){
  source(paste0(DR_code_S,"Ensemble_SDM/3_BM2_FB_SDM_EM_projection_with_crop4.r")) #this is where all configurations are at
}

