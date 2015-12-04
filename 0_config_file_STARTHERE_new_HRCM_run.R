rm(list = ls()) #remove all past worksheet variables
#options(error=stop) #this keeps the code from running after errors 

library(stringr)

###################################
####SET SOURCE LOCATION############
###################################
#Assigns source file (directory_registry file) to appropriate hardware (according to machine number)
rootDir = "D:/"
#location of repository:
codeDir = paste0(rootDir, "Dropbox/code/", "Ensemble_SDM/")  
#were to place project directory where all output will be saved to:
resultsDir = paste0(rootDir, "PICCC_analysis/","FB_analysis/model_results/biomod2/") 
#location of example_data/necessary_run_data/
necessary_run_data = paste0(resultsDir,"necessary_run_data/") #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)

###Location for clim predictors
#bioclimDataDir is location of example_data/clim_data/
#bioclimDataDir = paste0(resultsDir,"example_data/clim_data/") 
bioclimDataDir ="D:/PICCC_data/climate_data/HRCM_151110/complete_rasters/allYrs_avg/bioclims_abs/"
fitting_clim_data_dir = paste0(bioclimDataDir,"all_baseline/125m/") 
clim_data_2000=paste0(bioclimDataDir,"all_baseline/500m/")
clim_data_2100=paste0(bioclimDataDir,"all_future/500m/")
#clim_data_2100rev="D:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg_revDelta_rounded/bioclims_abs/all_future/500m/"

###################################
####GENERAL MODEL CONFIGURATION####
###################################
project_name = "FB_new_HRCM_projection_few_runs" #assign project name to the current run
#choose species of interest. Due to data sensitivity, pseudo data is provided for Akekee and Akikiki as example 
spp_nm = c('Akiapolauu', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Iiwi', 'Akekee', 'Akikiki', 'Puaiohi', 'Akohekohe', 'Maui_Alauahio', 'Maui_Parrotbill')#, , 'Akikiki', 'Anianiau',  'Oahu_Amakihi', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Oahu_Elepaio', 'Puaiohi', 'Kauai_Amakihi',  ) #'Amakihi', 'Elepaio', 

#Biomod2 modelling options
models_to_run = c("GBM","MAXENT") #choose biomod2 models to run - possibilities are: 'GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT' 
eval_stats = c("ROC") #choose evaluation methods - possibilties are: 'KAPPA','TSS','ROC'
env_var_files = c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif") #choose bioclimatic variables of interest - if using new clim data use ".tif" instead
plot_graphs = F #plot graphs of results (T) or not (F)
apply_biomod2_fixes = T #if running large models use this option - solves memory problems
overwrite = F #T if want to overwrite past results
paralelize = T #turn on multi instance execution (1 species per thread)
cpucores=10 #set to a very high number to use total number of threads available in computer

# main scripts (options below)
EM_fit = T #if you want to run the model fitting = T
EM_ensemble = T  #if you want to run the ensemble modelling = T
EM_project = T #if you want to project the model results = T

#optional scripts (options below)
merge_all_var_importance_and_model_eval = T
model_eval_graph = T
var_importance_graph = T
create_response_curves = T
#scripts below need two projections for comparison (typically baseline and future)
create_analog_climate_map = T
raster_output_creation = T #additional/ overridding configurations within files
distribution_shift_calculations = T 
spp_ensemble_maps = T

#################################
####CONFIG FOR SPECIFIC STEPS####
#################################
####fit config (script#1)
NbRunEval = 6 #number of evaluation runs for ensemble modeling
include_Abs = F #in test phase
PseudoAbs_outside_CE = F #if T, will only consider Pseudo Absences outside climate envelope of all points collected
dens_PAs_outside_CE = 1 #if 1 will create PA density that is equal to point density within surveyed areas
PA.nb.rep = 6
PA.nb.absences = 1000 #asssign number of Pseudo absence points (if PseudoAbs_outside_CE = T, this will be overridden! (n of PAs will be determined by P/A point density within CE)) 
candidatePAperPA = 200 #only used if if PAs_outside_CE = F, if value == 0, will use PA.nb.absences   
PA.strategy = "random" #strategy for selecting pseudo absences ('random', 'sre', 'disk' or 'user.defined')
equiv_100m = 0.0009430131
PA.dist.min = 5*equiv_100m #500m min distance from actual data points - only for 'disk' absences selection
do.full.models = T

####ensemble config (script#2)
eval.metric.threshold = rep(0.5,length(eval_stats)) #sets the minimum scores below which models will be excluding when building ensembles

####projection config (script#3)
proj_by_island=F #this is an option to project our models by island, to avoid memory issues.
baseline_or_future = 4 #1 for baseline, 4 for future
clampingMask = F #if T clamping mask will be saved
memory = T #keep.in.memory = memory; if T and clamping Mask = T, clamping mask will be saved to hard drive 

#create_analog_climate_map config
toCompareWithCurrentClimate=4 #4 for future

#raster_output_creation config, shift calculation parameters, multi spp maps
spp_ensemble_type="wmean" #for raster creation/ shift calc/ multi spp maps 
spp_ensemble_eval_stats=c('ROC') #for raster creation/ shift calc/ multi spp maps

comp_projects=c('baseline', 'future') #for raster creation/ shift calc
plot_spp_ensemble_CV=T #for raster creation
masked_spp_ensemble_map=FALSE #for raster creation
model_resolution=0.5 #inkm #for shift calc
exclude_areas_beyond_primary_habitat=T #for shift calc
habitat_overlay=T #for multi species maps
BPS=F #for multi species maps #overlay with prehistorical habitat distribution (landfire BPS)
veg_overlay_eval_stat="ROC"

##########################
####RUNNING SCRIPTS!!!####
##########################
working_dir = paste0(resultsDir, project_name, "/") #assign working directory
crop_raster_dir = paste0(working_dir, "map_crop/") #assign directory for cropped raster files
csv_dir = paste0(working_dir,"single_sp_CSVs/") #assign directory for single species CSV's
dir.create(working_dir, showWarnings = F) #creates working directory if missing
dir.create(paste0(working_dir,"output_rasters/"), showWarnings = F, recursive = T) #creates working directory if missing
setwd(working_dir)

##Plotting options depending on if server or not
useRasterDef = TRUE
interpolateDef = FALSE

dir_for_temp_files <- paste0(rootDir,'temp/', project_name, "/", baseline_or_future, "/") #dir for temp run data (to avoid memory errors)

if (apply_biomod2_fixes){
  maxentWDtmp = paste0("maxentWDtmp_", baseline_or_future)
  dir.create(dir_for_temp_files, showWarnings = F, recursive = T)
}

#assigning which projected climate data set to use depending on scenario
if (baseline_or_future == 1){
  clim_data_dir = clim_data_2000 
  proj_nm = 'baseline'}
if (baseline_or_future == 2){
  clim_data_dir = clim_data_2000wettest
  proj_nm ='baseline_wettest'}
if (baseline_or_future == 3){
  clim_data_dir = clim_data_2000driest
  proj_nm = 'baseline_driest'}
if (baseline_or_future == 4){
  clim_data_dir = clim_data_2100 
  proj_nm = 'future'}
if (baseline_or_future == 5){
  clim_data_dir = clim_data_2100wettest
  proj_nm = 'future_wettest'}
if (baseline_or_future == 6){
  clim_data_dir = clim_data_2100driest
  proj_nm = 'future_driest'}
if (baseline_or_future == 7){
  clim_data_dir = clim_data_2100rev
  proj_nm = 'future_rev'}

####Runs script for model fitting, creating ensemble models, and projecting models according to settings above.
##start the clock
ptmOverallStart <- proc.time()

if (EM_fit){ #runs fitting code
  source(paste0(codeDir,"1_BM2_FB_SDM_fitting.r")) 
}
if (EM_ensemble){ #runs ensemble code
  source(paste0(codeDir,"2_BM2_FB_SDM_EM_creation.r")) 
}
if (EM_project){ #runs projection code
  source(paste0(codeDir,"3_BM2_FB_SDM_EM_projection_byIsland.r"))   
}

#auxiliary scripts
if (merge_all_var_importance_and_model_eval){source(paste0(codeDir,"1opt_merge_all_var_importance_and_model_eval.R"))}
if (model_eval_graph){source(paste0(codeDir,"1opt_model_eval_graph.R"))}
if (var_importance_graph){source(paste0(codeDir,"1opt_var_importance_graph.R"))}
if (create_response_curves){source(paste0(codeDir,"2opt_BM2_FB_SDM_response_curves.r"))}
if (create_analog_climate_map){source(paste0(codeDir,"4opt_create_analog_climate_map.R"))}
if (raster_output_creation){source(paste0(codeDir,"4_SDM_raster_output_creation_newBiomVer.R"))}
if (distribution_shift_calculations){source(paste0(codeDir,"6_distribution_shift_calculations_newbiomod2.r"))}
if (spp_ensemble_maps){source(paste0(codeDir,"7_spp_ensemble_maps.r"))}


##stop the clock
ptmOverallElaps = proc.time() - ptmOverallStart #calculates time it took to run all code
jnk = as.numeric(ptmOverallElaps[3])/length(spp_nm) #assigns temporary variable to the numeric value of the time elapsed per species
jnk = jnk/60 #converts elapsed time into minutes
cat('\n','It took ', jnk, "minutes (on average) to model each species with",
    length(models_to_run), "model types") #sign-posting

