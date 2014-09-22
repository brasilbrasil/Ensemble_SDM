###USER CONFIGURATION
#see 0_sdm_config file.r

####START UNDERHOOD
setwd(working_dir) #sets working directory

##Loading package libraries
library(biomod2)
library(stringr)

#memory.limit(size = 4095) #increases memory limit size
sp_nm = spp_nm[1] #resets so the first species to run is the first one listed in config file or csv
for (sp_nm in spp_nm){
  sp_nm = as.character(sp_nm) #defines the species name as a character string - not needed if it is already a text name
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
  
    jnk = myBiomodModelEval[,2,,,]
    NAs = which(is.na(jnk))
    all_models = list()
    if (length(dimnames(jnk)) == 4) {
      for (d in dimnames(jnk)[[4]]) {
        for (c in dimnames(jnk)[[3]]) {
          for (b in dimnames(jnk)[[2]]) {
            for (a in dimnames(jnk)[[1]]) {
              jnk_str = paste(sp_nm,d,c,b,a, sep = "_")
              jnk_str2 = paste(sp_nm,d,c,b, sep = "_")
              all_models[length(all_models)+1] = jnk_str2
            }
          }
        }
      }
    } else {
      #turn the "NbRunEval" variable into the Run names (e.g. a NbRunEval = 2 would run as "RUN1", "RUN2", "Full")
      RUNnames = c()
      for (RUNnum in 1: NbRunEval) {
        RUNstr = paste0("RUN", RUNnum)
        RUNnames = c(RUNnames, RUNstr)
      }
      RUNnames = c(RUNnames, "Full")
      
      #generate the PA# variable names by adding "PA" before each repetition
      PAnames = c()
      for (PAnum in 1:PA.nb.rep) {
        PAstr = paste0("PA", PAnum)
        PAnames = c(PAnames, PAstr)  
      }
      
      for (d in PAnames) {
        for (c in RUNnames){
          for (b in models_to_run) {
            for (a in eval_stats) {
              jnk_str2 = paste(sp_nm,d,c,b, sep = "_")
              all_models[length(all_models)+1] = jnk_str2
            }
          }
        }
      }
    }
    
    bad_models_short = all_models[NAs]
    bad_models_short = unique(bad_models_short)
    jnk_good =! (all_models %in% bad_models_short)
    remaining_models = all_models[jnk_good]
    remaining_models = unlist(unique(remaining_models))    

    
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
}