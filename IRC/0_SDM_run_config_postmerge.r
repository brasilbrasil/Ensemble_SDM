rm(list = ls()) #remove all past worksheet variables
options(error = stop) #this keeps the code from running after errors 

###################################
####SET SOURCE LOCATION############
###################################

machine = 3 #use 1 for PICCC server, 2 for LF local drive, 3 for IRC local drive

#Assigns source file (directory_registry file) to appropriate hardware (according to machine number)
if (machine == 1){
  source("Y:/PICCC_analysis/code/","directory_registry.r")
} else {
  if (machine == 2) {
    source("C:/Users/lfortini/directory_registry.r") #not sure if this is the right location for Lucas's local file
  } else {
    if (machine == 3) {
      source("C:/USGS_Honolulu/PICCC_code/Ensemble_SDM/IRC/directory_registryIRC.r")
    } else {
      cat('\n','Error - invalid machine number')
    }
  }
}

###################################
####GENERAL MODEL CONFIGURATION####
###################################

#setting file locations 
project_name = "FB_test20130904" #assign project name to the current run
working_dir = paste(resultsDir, project_name, sep = "/") #assign working directory
crop_raster_dir = paste0(working_dir, "/map_crop") #assign directory for cropped raster files
csv_dir = paste0(working_dir,"/single_sp_CSVs") #assign directory for single species CSV's

#fitting_clim_data_dir=paste0(DR_FB_clim_data,"all_grd/all_baseline/100m/") #MOVED TO REGISTRY 
#necessary_run_data=paste0(DR_FB_SDM_results_S,'necessary_run_data/') #MOVED TO REGISTRY


if (file.exists(working_dir) == F){ 
  dir.create(working_dir, showWarnings = F) #creates working directory if missing
}

#choose species of interest - all (from CSV file) or subset listed
run_all_spp = 0 #if running all species enter "1" and if only subset enter "0"
spp_subset = c("Akekee","Kauai_Amakihi") # "Oahu_Amakihi","Hawaii_Akepa", "Palila") #if only subset, enter spp names here 

#spp_nm=(read.csv(paste0(resultsDir,'spp_to_run_all.csv'),header=F, stringsAsFactors=F)) #alternative way of running all species


#Biomod2 modelling options for species of interest
models_to_run = c('GBM','MAXENT') #choose biomod2 models to run - possibilities are: 'GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT' 
eval_stats = c('ROC') #choose evaluation methods - possibilties are: 'KAPPA','TSS','ROC'
env_var_files = c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") #choose bioclimatic variables of interest
plot_graphs = 1 #plot graphs of results (=1) or not (=0)
EM_fit = T #if you want to run the model fitting = T
EM_ensemble = T #if you want to run the ensemble modelling = T
EM_project = F #if you want to project the model results = T
apply_biomod2_fixes = F #if running large models use this option - solves memory problems
memory.limit(size = 4000) #increases memory allocation
overwrite=0

#Assigns the species names either according to csv file (all) or list
if (run_all_spp == 1){
  spp_nm = read.csv(allSppNames, header = F, stringsAsFactors = F)
} else {
  spp_nm = spp_subset
}

#################################
####CONFIG FOR SPECIFIC STEPS####
#################################
####fit config (script#1)
NbRunEval = 2 #number of evaluation runs for ensemble modeling
include_Abs = T #in test phase
PseudoAbs_outside_CE = T #if T, will only consider Pseudo Absences outside climate envelope of all points collected
dens_PAs_outside_CE=1 #if 1 will create PA density that is equal to point density within surveyed areas
PA.nb.rep = 2
PA.nb.absences = 1000 #asssign number of Pseudo absence points (if PseudoAbs_outside_CE = T, this will be overridden! (n of PAs will be determined by P/A point density within CE)) 
PA.strategy = "random" #strategy for selecting pseudo absences ('random', 'sre', 'disk' or 'user.defined')
equiv_100m = 0.0009430131
PA.dist.min = 5*equiv_100m #500m min distance from actual data points - only for 'disk' absences selection
do.full.models = T

####ensemble config (script#2)
eval.metric.threshold = rep(0.5,length(eval_stats)) #sets the minimum scores below which models will be excluding when building ensembles

####projection config (script#3)
baseline_or_future = 1 #1 for baseline, 4 for future
memory = T #keep.in.memory=memory
dir_for_temp_files <- paste(rootDir,'temp', project_name, baseline_or_future, sep = "/") #dir for temp run data (to avoid memory errors)

##########################
####RUNNING SCRIPTS!!!####
##########################
if (apply_biomod2_fixes){
maxentWDtmp = paste0("maxentWDtmp_", baseline_or_future)
dir.create(dir_for_temp_files, showWarnings=F, recursive=T)
}

###not in FWS code (multi instance automation)
#this code below will subset species into the right number of instances started with the bat file                        
Sys.sleep(6) #time for script process to show up on tasklist
n_instances=length(system('tasklist /FI "IMAGENAME eq Rscript.exe" ', intern = TRUE))-3
rsession_instances=length(system('tasklist /FI "IMAGENAME eq rsession.exe" ', intern = TRUE))-3
cpucores=as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
if (n_instances>0 & cpucores>1 & rsession_instances<1){
  n_instances=1
  jnkn=length(spp_nm)
  x=c(1:jnkn)
  chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
  groups=chunk(x,cpucores)
  jnk=groups[n_instances][[1]]
  spp_nm=spp_nm[jnk]
}


####Runs script for model fitting, creating ensemble models, and projecting models according to settings above.
if (EM_fit){
  source(paste0(codeDir,"/IRC/1_BM2_FB_SDM_fitting_w_cropping4_postmerge.r")) #this is where all configurations are at
}
if (EM_ensemble){
  source(paste0(codeDir,"/IRC/2_BM2_FB_SDM_EM_fitting2_IRC.r")) #this is where all configurations are at
}
if (EM_project){
  source(paste0(codeDir,"/IRC/3_BM2_FB_SDM_EM_projection_with_crop4_IRC.r")) #this is where all configurations are at
}

##add bit here at the end to save config file copy (renamed) in folder
