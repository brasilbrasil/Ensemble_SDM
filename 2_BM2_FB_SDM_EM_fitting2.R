###USER CONFIGURATION
#see 0_sdm_config file.r

####START UNDERHOOD
setwd(working_dir)

library(biomod2)
library(stringr)
#sp_nm="Akepa" #debug

memory.limit(size=4095)
sp_nm=spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'ensemble creation...')
  sp_nm0=sp_nm
  workspace_name=paste(sp_nm0,"_FB_modelfitting.RData", sep = "") #set name of file to load workspace data from model run
  workspace_name_out=paste(sp_nm0,"_FB_EM_fit.RData", sep = "") #set name of file to load workspace data from model run
  if (file.exists(workspace_name_out)==F | overwrite==1){
    load(workspace_name)
    sp_nm=str_replace_all(sp_nm,"_", ".")
    
    spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))
    ###################################################
    ### code chunk number 8: modeling_summary
    ###################################################
    myBiomodModelOut
    
    
    ###################################################
    ### code chunk number 9: modeling_model_evaluation
    ###################################################
    # get all models evaluation
    myBiomodModelEval <- getModelsEvaluations(myBiomodModelOut)
    
    ###################################################
    ### code chunk number 10: modeling_variable_importance
    ###################################################
    # print variable importances
    getModelsVarImport(myBiomodModelOut)
    
    ###################################################
    ### code chunk number 11: ensemble_modeling
    ###################################################
    myBiomodEM <- BIOMOD_EnsembleModeling(
      modeling.output = myBiomodModelOut,
      chosen.models = 'all', #these are not model types (e.g., GBM), but model runs (e.g., PA1_RF)
      em.by='all',
      eval.metric = eval_stats, #c('TSS', 'ROC', 'KAPPA'); 'all', #c('TSS', 'ROC'),
      eval.metric.quality.threshold = eval.metric.threshold,
      prob.mean = T,
      prob.cv = T,
      prob.ci = T,
      prob.ci.alpha = 0.05,
      prob.median = T,
      committee.averaging = T,
      prob.mean.weight = T,
      prob.mean.weight.decay = 'proportional' )
    cat('\n',sp_nm,'ensemble done...')

    ###################################################
    ### code chunk number 12: ensemble_modeling_outputs
    ###################################################
    # print summary
    myBiomodEM
    # get evaluation scores
    getEMeval(myBiomodEM)    
    save("myBiomodEM", "myBiomodModelOut", file=workspace_name_out)   #save workspace
    
    }else{
    cat('\n',sp_nm,'ensemble previously done...')
  }
}
