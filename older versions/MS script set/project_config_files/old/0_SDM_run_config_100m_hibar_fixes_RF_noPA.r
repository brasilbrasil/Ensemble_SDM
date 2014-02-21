rm(list = ls()) #remove all past worksheet variables
source(paste0("Y:/PICCC_analysis/code/","directory_registry.r"))

###################################
####GENERAL MODEL CONFIGURATION####
###################################
local_config_dir=DR_FB_SDM_results_S #'C:/Users/lfortini/'
#spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
spp_nm=c("Akekee", "Kauai_Amakihi", "Oahu_Amakihi","Hawaii_Akepa", "Palila")#, "Kauai_Amakihi", "Hawaii_Akepa", "Palila")
project_name='test_runs_100m_hibar_fixes_RF_noPA'
server=1
overwrite=0
models_to_run=c('GBM','MAXENT', 'RF')
eval_stats=c("ROC")
plot_graphs=1
EM_fit=T
EM_ensemble=T
EM_project=T
#memory.limit(size=24000000)
apply_biomod2_fixes=T #if running large models use this option


if (server==1){
  working_dir=paste0(DR_FB_SDM_results_S,project_name,'/')
  fitting_clim_data_dir=paste0(DR_FB_clim_data,"all_grd/all_baseline/100m/") 
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
NbRunEval=100

####ensemble config

####projection config
baseline_or_future=1 #1 for baseline, 4 for future
memory = T #keep.in.memory=memory
dir_for_temp_files<-paste('Y:/temp/', project_name,'/', baseline_or_future, '/', sep='') #dir for temp run data (to avoid memory errors)

if (server==1){
  clim_data_2000=paste0(DR_FB_clim_data,"all_grd/all_baseline/250m/")
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
if (apply_biomod2_fixes){
  maxentWDtmp = paste("maxentWDtmp_", baseline_or_future, sep = "")
  dir.create(dir_for_temp_files, showWarnings=F, recursive=T)
}#dir.create(paste('Y:/temp/', project_name,'/', sep=''), showWarnings=F)
dir.create(working_dir, showWarnings=F)

if (EM_fit){
  source(paste0(DR_code_S,"Ensemble_SDM/Em_code_before_PA_correction/1_BM2_FB_SDM_fitting_w_cropping4.r")) #this is where all configurations are at
}
if (EM_ensemble){
  source(paste0(DR_code_S,"Ensemble_SDM/Em_code_before_PA_correction/2_BM2_FB_SDM_EM_fitting2.r")) #this is where all configurations are at
}
if (EM_project){
  source(paste0(DR_code_S,"Ensemble_SDM/Em_code_before_PA_correction/3_BM2_FB_SDM_EM_projection_with_crop4.r")) #this is where all configurations are at
}

