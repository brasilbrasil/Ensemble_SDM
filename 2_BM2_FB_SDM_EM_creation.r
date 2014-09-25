###USER CONFIGURATION
#see 0_sdm_config file.r

####START UNDERHOOD
setwd(working_dir) #sets working directory
require(snowfall)

#memory.limit(size = 4095) #increases memory limit size
sp_nm = spp_nm[1] #resets so the first species to run is the first one listed in config file or csv
sp_parallel_run=function(sp_nm){  ##Loading package libraries
  library(biomod2)
  library(stringr)
  sp_nm = as.character(sp_nm) #defines the species name as a character string - not needed if it is already a text name
  sp_dir = str_replace_all(sp_nm,"_", ".") #replaces "_" with "." in sp_nm
  sink(file(paste0(working_dir,sp_dir,"/",sp_dir,Sys.Date(),"_EM_creation_log.txt"), open="wt"))#######NEW
  cat('\n',sp_nm,'ensemble creation...')
  workspace_name = paste0(sp_nm,"_FB_modelfitting.RData") #set name of file to load workspace data from model run
  workspace_name_out = paste0(sp_nm,"_FB_EM_fit.RData") #set name of file to load workspace data from model run
  if (file.exists(workspace_name_out) == F | overwriteData == T){ #run only if file does not already exist and overwrite is turned off (in config file)
    # Start the clock!
    ptmModule2Start <- proc.time()
    
    load(workspace_name) #loads the model results from the species of interest
    sp_nm = str_replace_all(sp_nm,"_", ".") #replaces any "_" with "." in the species name
    spp_info = read.csv(paste0(csv_dir, "FB_spp_data.csv")) #creates data frame from species info csv file
    
    ###################################################
    ### code chunk number 8: modeling_summary
    ###################################################
    myBiomodModelOut #returns summary of biomod2 model fit run for species
    myBiomodModelOut[models.computed]
    #showClass("myBiomodModelOut")
    
    ###################################################
    ### code chunk number 9: modeling_model_evaluation
    ###################################################
    # get all models evaluation
    myBiomodModelEval <- get_evaluations(myBiomodModelOut) #creates an array with model evaluation results for each species model
    
    ###################################################
    ### code chunk number 10: modeling_variable_importance
    ###################################################
    # print variable importances
    get_variables_importance(myBiomodModelOut) #returns an array with model variable importances (i.e bio1, bio7, etc)
    
    ###################################################
    ###new code- remove models with bad cutoffs
    myBiomodModelEval_sum_cutoffs=apply(myBiomodModelEval, c(2,3,4,5), sum)
    myBiomodModelEval_sum_cutoffs=c(myBiomodModelEval_sum_cutoffs[2, , , ])
    good_cutoffs = which(!is.na(myBiomodModelEval_sum_cutoffs))
    models.computed=myBiomodModelOut@models.computed
    remaining_models=models.computed[good_cutoffs]
    
    ###################################################
    ### code chunk number 11: ensemble_modeling
    ###################################################
    #combines models and make ensemble predictions built with BIOMOD_Modeling in module 1
    myBiomodEM <- BIOMOD_EnsembleModeling( 
      modeling.output = myBiomodModelOut, #the output from BIOMO_Modeling in previous step
      chosen.models = remaining_models, #these are not model types (e.g., GBM), but model runs (e.g., PA1_RF)
      em.by='all', #Available values are 'PA_dataset+repet' (default), 'PA_dataset+algo', 'PA_dataset', 'algo' and 'all'
      eval.metric = eval_stats, #c('TSS', 'ROC', 'KAPPA'); 'all', #c('TSS', 'ROC'),
      eval.metric.quality.threshold = eval.metric.threshold,
      prob.mean = T,
      prob.cv = T,
      prob.ci = T,
      prob.ci.alpha = 0.05, #signficance level for estimating confidence interval
      prob.median = T,
      committee.averaging = T,
      prob.mean.weight = T,
      prob.mean.weight.decay = 'proportional' ) #defines relative importance of weights; default is 'proportional'
    cat('\n',sp_nm,'ensemble done...')

    ###################################################
    ### code chunk number 12: ensemble_modeling_outputs
    ###################################################
    # print summary
    myBiomodEM #returns summary of ensemble modeling
    
    # get evaluation scores
    get_evaluations(myBiomodEM) #returns evaluation stats (testing.data, cutoff, sensitivity, and specificity) for mean, cv, etc.    
    save("myBiomodEM", "myBiomodModelOut","remaining_models", file = workspace_name_out)   #save workspace
    
    #Stop the clock
    ptmModule2Elaps = proc.time() - ptmModule2Start #calculates time it took to run all code
    jnk = as.numeric(ptmModule2Elaps[3]) #assigns temporary variable to the numeric value of the time elapsed
    jnk = jnk/60 #converts elapsed time into minutes
    cat('\n','It took ', jnk, "minutes to ensemble model", sp_nm) #sign-posting
    
    }else{
    cat('\n',sp_nm,'ensemble previously done...') #if file already exists for this run
  }
  sink(NULL)
}

if (is.null(cpucores)){
  cpucores=as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))  
}else{
  cpucores=min(cpucores, as.integer(Sys.getenv('NUMBER_OF_PROCESSORS')))
}
sfInit( parallel=T, cpus=cpucores) # 
sfExportAll() 
system.time((sfLapply(spp_nm,fun=sp_parallel_run)))
#system.time(sfClusterApplyLB(iter_strings,fun=sp_parallel_run)) #why not alway us LB? Reisensburg2009_TutParallelComputing_Knaus_Porzelius.pdf
sfRemoveAll()
sfStop()
