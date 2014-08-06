rm(list = ls()) #remove all past worksheet variables
#options(error=stop) #this keeps the code from running after errors 

library(stringr)

###################################
####SET SOURCE LOCATION############
###################################
machine = 3 #use 1 for LF laptop, 2 for LF hard drive, 3 for IRC local drive
#Assigns source file (directory_registry file) to appropriate hardware (according to machine number)
if (machine == 1){
  source(paste0("C:/Users/lfortini/","directory_registry.r"))
  server = F
} else {
  if (machine == 2) {
    source(paste0("C:/Users/lfortini/","directory_registry.r")) #actually same as above
    server = F
  } else {
    if (machine == 3) {
      source(paste0("C:/USGS_Honolulu/PICCC_code/Ensemble_SDM/IRC/directory_registryIRC_20131120.r"))
      server = F
    } else {
      cat('\n','Error - invalid machine number')
    }
  }
}

###################################
####GENERAL MODEL CONFIGURATION####
###################################
#setting file locations 
project_name = "FB_test20131216" #assign project name to the current run

#choose species of interest - all (from CSV file) or subset listed
run_all_spp = F #if running all species enter "T" and if only subset enter "F"
spp_subset = c("Kauai_Amakihi","Akekee") # "Oahu_Amakihi","Hawaii_Akepa", "Palila") #if only subset, enter spp names here 

#Biomod2 modelling options for species of interest
models_to_run = c("GBM","MAXENT") #choose biomod2 models to run - possibilities are: 'GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT' 
eval_stats = c("ROC", "KAPPA") #choose evaluation methods - possibilties are: 'KAPPA','TSS','ROC'
env_var_files = c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") #choose bioclimatic variables of interest - if using new clim data use ".tif" instead
plot_graphs = T #plot graphs of results (T) or not (F)
EM_fit = F #if you want to run the model fitting = T
EM_ensemble = F  #if you want to run the ensemble modelling = T
EM_project = T #if you want to project the model results = T
create_response_curves = F
apply_biomod2_fixes = F #if running large models use this option - solves memory problems
overwriteData = F #T if want to overwrite and F if not
paralelize = F #turn on multi instance auto start

#################################
####CONFIG FOR SPECIFIC STEPS####
#################################
####fit config (script#1)
NbRunEval = 2 #number of evaluation runs for ensemble modeling
include_Abs = T #in test phase
PseudoAbs_outside_CE = T #if T, will only consider Pseudo Absences outside climate envelope of all points collected
dens_PAs_outside_CE = 1 #if 1 will create PA density that is equal to point density within surveyed areas
PA.nb.rep = 2
PA.nb.absences = 1000 #asssign number of Pseudo absence points (if PseudoAbs_outside_CE = T, this will be overridden! (n of PAs will be determined by P/A point density within CE)) 
candidatePAperPA = 50 #only used if if PAs_outside_CE = F, if value == 0, will use PA.nb.absences   
PA.strategy = "random" #strategy for selecting pseudo absences ('random', 'sre', 'disk' or 'user.defined')
equiv_100m = 0.0009430131
PA.dist.min = 5*equiv_100m #500m min distance from actual data points - only for 'disk' absences selection
do.full.models = T

####ensemble config (script#2)
eval.metric.threshold = rep(0.5,length(eval_stats)) #sets the minimum scores below which models will be excluding when building ensembles

####projection config (script#3)
baseline_or_future = 4 #1 for baseline, 4 for future
clampingMask = F #if T clamping mask will be saved
memory = T #keep.in.memory = memory; if T and clamping Mask = T, clamping mask will be saved to hard drive 


##########################
####RUNNING SCRIPTS!!!####
##########################

working_dir = paste0(resultsDir, project_name, "/") #assign working directory
crop_raster_dir = paste0(working_dir, "map_crop/") #assign directory for cropped raster files
csv_dir = paste0(working_dir,"single_sp_CSVs/") #assign directory for single species CSV's
dir.create(working_dir, showWarnings = F) #creates working directory if missing

#Assigns the species names either according to csv file (all) or list
if (run_all_spp){
  spp_nm = read.csv(allSppNames, header = F, stringsAsFactors = F)
} else {
  spp_nm = spp_subset
}

##Plotting options depending on if server or not
if (server == TRUE){
  useRasterDef = FALSE
  interpolateDef = TRUE
} else {
  useRasterDef = TRUE
  interpolateDef = FALSE
}

dir_for_temp_files <- paste0(rootDir,'/temp/', project_name, "/", baseline_or_future, "/") #dir for temp run data (to avoid memory errors)

if (apply_biomod2_fixes){
  maxentWDtmp = paste0("maxentWDtmp_", baseline_or_future)
  dir.create(dir_for_temp_files, showWarnings = F, recursive = T)
}

#this code below will subset species into the right number of instances started with the bat file                        
n_instances = length(list.files(working_dir, pattern = "^00instance"))
cpucores = as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
if (paralelize){
  if (cpucores > length(spp_nm)){cpucores = length(spp_nm)}
  jnkn = length(spp_nm)
  x = c(1:jnkn)
  chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
  groups = chunk(x,cpucores)
  jnk = groups[n_instances+1][[1]]
  spp_nm = spp_nm[jnk]
  spp_str = ""
  for (sp_nm in spp_nm){
    spp_str = paste(spp_str, sp_nm, sep = "__")
  }
  time = Sys.time()
  time = str_replace_all(time,":", ".")
  instance_file = paste0("00instance", spp_str, "_", time)
  file.create(paste0(working_dir, "/", instance_file), showWarnings=F)  
}


####Runs script for model fitting, creating ensemble models, and projecting models according to settings above.
##start the clock
ptmOverallStart <- proc.time()

if (EM_fit){ #runs fitting code
  source(paste0(codeDir,"1_BM2_FB_SDM_fitting_w_cropping4_newPackageFunctions.r")) 
}
if (create_response_curves){
  source(paste0(codeDir,"2opt_BM2_FB_SDM_response_curves3.r"))
}
if (EM_ensemble){ #runs ensemble code
  source(paste0(codeDir,"2_BM2_FB_SDM_EM_fitting2_IRC_newPackageFunctions.r")) 
}
if (EM_project){ #runs projection code
  source(paste0(codeDir,"3_BM2_FB_SDM_EM_projection_with_crop4_IRC.r")) 
}

##stop the clock
ptmOverallElaps = proc.time() - ptmOverallStart #calculates time it took to run all code
jnk = as.numeric(ptmOverallElaps[3])/length(spp_nm) #assigns temporary variable to the numeric value of the time elapsed per species
jnk = jnk/60 #converts elapsed time into minutes
cat('\n','It took ', jnk, "minutes (on average) to model each species with",
    length(models_to_run), "model types") #sign-posting

if (paralelize){
  file.remove(paste0(working_dir, instance_file))
}
