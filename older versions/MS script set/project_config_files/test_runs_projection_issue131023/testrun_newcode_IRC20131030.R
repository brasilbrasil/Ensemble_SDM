rm(list = ls()) #remove all past worksheet variables
set.seed(123)
source(paste0("C:/USGS_Honolulu/PICCC_code/Ensemble_SDM/IRC/directory_registryIRC_newcode.r"))
options(error=stop) #this keeps the code from running after errors 

###################################
####GENERAL MODEL CONFIGURATION####
###################################
#local_config_dir=resultsDir
#spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
spp_nm =  c('Maui_Parrotbill')  #c('Palila', 'Hawaii_Akepa', 'Maui_Parrotbill', 'Puaiohi')
project_name ='test20131030_newcode'
server = 1
overwrite = 0
paralelize = F
models_to_run = c('GBM','MAXENT')
eval_stats = c("ROC", "KAPPA", "TSS")
plot_graphs = 1
EM_fit = F
EM_ensemble = F
EM_project = F
create_response_curves = F
apply_biomod2_fixes = F #if running large models use this option

if (server == 1){
  working_dir = paste0(resultsDir, '/', project_name,'/')
  #fitting_clim_data_dir=paste0(DR_FB_clim_data_2013,"all_baseline/125m/") 
  #necessary_run_data=paste0(resultsDir,'necessary_run_data/') #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/test/'
  necessary_run_data='C:/Users/lfortini/Data/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
  fitting_clim_data_dir="C:/Users/lfortini/Data/SDM_env_data/all_grd/all_baseline/100m/"
  DR_code_S=paste0(rootDir, "/Dropbox/code/") #HAD TO ADD THIS TO READ THIS FROM OLD DIRECTORY FILE 
  
}

env_var_files = c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") 
crop_raster_dir = paste0(working_dir, 'map_crop/')
csv_dir=paste0(working_dir,"single_sp_CSVs/")

#################################
####CONFIG FOR SPECIFIC STEPS####
#################################
####fit config (script#1)
NbRunEval = 4
include_Abs = T #in test phase
PAs_outside_CE = T #if T, will only consider PAs outside climate envelope of all points collected
dens_PAs_outside_CE = 1 #if 1 will create PA density that is equal to point density within surveyed areas
PA.nb.rep = 5
PA.nb.absences = 10000 #only used if if PAs_outside_CE=F, this will be overridden! (n of PAs will be determined by P/A point density within CE 
candidatePAperPA = 50 #only used if if PAs_outside_CE=F, if value ==0, will use PA.nb.absences   
PA.strategy = "random"
equiv_100m = 0.0009430131
PA.dist.min = 5*equiv_100m #500 min distance from actual data points 
do.full.models = T

####ensemble config (script#2)
eval.metric.threshold = rep(0.5,length(eval_stats))

####projection config (script#3)
baseline_or_future = 1 #1 for baseline, 4 for future
memory = T #keep.in.memory=memory
dir_for_temp_files <- paste(rootDir,'/temp/', project_name,'/', baseline_or_future, '/', sep='') #dir for temp run data (to avoid memory errors)

if (server==1){
  #clim_data_2000=paste0(DR_FB_clim_data_2013,"all_baseline/500m/")
  #clim_data_2100=paste0(DR_FB_clim_data_2013,"all_future/500m/")
  #clim_data_2000wettest="D:/GIS_Data/REnviroLayers/mixed_data_2000_250mwettest/"
  #clim_data_2000driest= "D:/GIS_Data/REnviroLayers/mixed_data_2000_250mdriest/"
  #clim_data_2100wettest="D:/GIS_Data/REnviroLayers/mixed_data_2100_250mwettest/"
  #clim_data_2100driest= "D:/GIS_Data/REnviroLayers/mixed_data_2100_250mdriest/"  
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

#this code below will subset species into the right number of instances started with the bat file                        
Sys.sleep(6) #time for script process to show up on tasklist
n_instances =length(list.files(working_dir, pattern="^00instance"))
cpucores = as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
#n_instances=length(system('tasklist /FI "IMAGENAME eq Rscript.exe" ', intern = TRUE))-3
#rsession_instances=length(system('tasklist /FI "IMAGENAME eq rsession.exe" ', intern = TRUE))-3
#if (n_instances>0 & cpucores>1 & rsession_instances<1 & paralelize){
if (paralelize){
  if (cpucores>length(spp_nm)){cpucores=length(spp_nm)}
  jnkn = length(spp_nm)
  x=c(1:jnkn)
  chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
  groups=chunk(x,cpucores)
  jnk=groups[n_instances+1][[1]]
  spp_nm=spp_nm[jnk]
  spp_str=""
  for (sp_nm in spp_nm){
    spp_str=paste(spp_str,sp_nm,sep="__")
  }
  time=Sys.time()
  time=str_replace_all(time,":", ".")
  instance_file=paste0("00instance",spp_str,"_",time)
  file.create(paste0(working_dir,instance_file),showWarnings=F)  
}

if (EM_fit){
  source(paste0(codeDir,"/1_BM2_FB_SDM_fitting_w_cropping4_postmerge.r")) #this is where all configurations are at
}
if (create_response_curves){
  source(paste0(codeDir,"/2opt_BM2_FB_SDM_response_curves3.r"))
}
if (EM_ensemble){
  source(paste0(codeDir,"/2_BM2_FB_SDM_EM_fitting2_IRC.r")) #this is where all configurations are at
}
if (EM_project){
  source(paste0(codeDir,"/3_BM2_FB_SDM_EM_projection_with_crop4_IRC.r")) #this is where all configurations are at
}

if (paralelize){
  file.remove(paste0(working_dir,instance_file))
}
