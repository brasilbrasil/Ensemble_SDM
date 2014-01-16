rm(list = ls()) #remove all past worksheet variables

###################################
####GENERAL MODEL CONFIGURATION####
###################################
#local_config_dir=resultsDir
#spp_nmS = c('Akekee', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi', 'Kauai_Amakihi', 'Hawaii_Amakihi', 'Apapane', 'Amakihi', 'Elepaio', 'Iiwi')
spp_nmS = c('Akekee', 'Hawaii_Akepa', 'Palila', 'Hawaii_Amakihi', 'Maui_Parrotbill')
spp_nm=spp_nmS[1]
sp_parallel_run=function(spp_nm){
  spp_nm=c(spp_nm)
  source(paste0("C:/Users/lfortini/","directory_registry.r"))
  project_name='finalmodel_P_PA_oldcode_less_PAs_all_models'
  server=1
  overwrite=0
  models_to_run=c('GLM','GBM','GAM','CTA','ANN', 'SRE','FDA','MARS','RF','MAXENT')
  eval_stats=c("ROC","TSS")
  plot_graphs=1
  EM_fit=T
  EM_ensemble=T
  EM_project=T
  create_response_curves=F
  apply_biomod2_fixes=T #if running large models use this option
  
  working_dir=paste0(resultsDir,project_name,'/')
  dir.create(working_dir,showWarnings=F)
  sink(file(paste0(working_dir,spp_nm[1],"_log.txt"), open="wt"))#######NEW
  
  env_var_files=c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif") 
  crop_raster_dir=paste(working_dir, 'map_crop/',sep="")
  csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")
  
  #################################
  ####CONFIG FOR SPECIFIC STEPS####
  #################################
  ####fit config (script#1)
  NbRunEval=3
  include_Abs=F #in test phase
  PAs_outside_CE=F #if T, will only consider PAs outside climate envelope of all points collected
  dens_PAs_outside_CE=1 #if 1 will create PA density that is equal to point density within surveyed areas
  PA.nb.rep=5
  PA.nb.absences = 10000 #only used if if PAs_outside_CE=F, this will be overridden! (n of PAs will be determined by P/A point density within CE 
  candidatePAperPA=100 #only used if if PAs_outside_CE=F, if value ==0, will use PA.nb.absences   
  PA.strategy = "random"
  equiv_100m=0.0009430131
  PA.dist.min = 5*equiv_100m #500 min distance from actual data points 
  do.full.models=T
  ####ensemble config (script#2)
  eval.metric.threshold = rep(0.5,length(eval_stats))
  
  ####projection config (script#3)
  baseline_or_future=1 #1 for baseline, 4 for future, 7 for hot
  memory = T #keep.in.memory=memory
  dir_for_temp_files<-paste(rootDir,'/temp/', project_name,'/', baseline_or_future, '/', sep='') #dir for temp run data (to avoid memory errors)
  
  if (server==1){
    clim_data_2004hottest=paste0(bioclimData2013Dir,"2004_hottest/500m/")
    clim_data_2104hottest=paste0(bioclimData2013Dir,"2104_hottest/500m/")
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
  library(stringr)
  if (apply_biomod2_fixes){
    maxentWDtmp = paste("maxentWDtmp_", baseline_or_future, sep = "")
    dir.create(dir_for_temp_files, showWarnings=F, recursive=T)
  }#dir.create(paste('Y:/temp/', project_name,'/', sep=''), showWarnings=F)
  dir.create(working_dir, showWarnings=F)
  
  setwd(working_dir)
  
  if (EM_fit){
    cat('starting fit', '\n')
    source(paste0(DR_code_S,"Ensemble_SDM/1_BM2_FB_SDM_fitting_w_cropping4.r"),local=T) #this is where all configurations are at
  }
  if (create_response_curves){
    cat('starting curves', '\n')
    source(paste0(DR_code_S,"Ensemble_SDM/2opt_BM2_FB_SDM_response_curves3.r"),local=T)
  }
  if (EM_ensemble){
    cat('starting ensemble', '\n')
    source(paste0(DR_code_S,"Ensemble_SDM/2_BM2_FB_SDM_EM_fitting2.r"),local=T) #this is where all configurations are at
  }
  if (EM_project){
    cat('starting projection', '\n')
    source(paste0(DR_code_S,"Ensemble_SDM/3_BM2_FB_SDM_EM_projection_with_crop4.r"),local=T) #this is where all configurations are at
  }
    
  sink.reset <- function(){
    for(i in seq_len(sink.number())){
      sink(NULL)
    }
  }
  sink.reset()
}


require(snowfall)
# Init Snowfall with settings from sfCluster
cpucores=as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
sfInit( parallel=T, cpus=cpucores) # 
sfExportAll() 
all_reps=sfLapply(spp_nmS,fun=sp_parallel_run)
sfRemoveAll( except=c( "rep_FX" ) )
sfStop()
