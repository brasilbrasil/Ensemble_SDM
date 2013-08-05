local_config_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m/') #'C:/Users/lfortini/'
#spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Anianiau", "Kauai_Amakihi", "Hawaii_Elepaio", "Palila")
spp_nm=c("Akekee")#, "Anianiau", "Kauai_Amakihi", "Oahu_Amakihi","Hawaii_Akepa", "Hawaii_Elepaio", "Palila")
server=1
overwrite=1
models_to_run=c('GBM','MAXENT')
eval_stats=c("ROC")
plot_graphs=1

if (server==1){
  working_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m_rounded/')
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
