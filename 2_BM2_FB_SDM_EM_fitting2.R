###USER CONFIGURATION
#see 0_sdm_config file.r

####START UNDERHOOD
setwd(working_dir)

library(biomod2)
library(stringr)
#sp_nm="Akepa" #debug

#memory.limit(size=24000000)

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
    ###new code- remove models with bad cutoffs
    jnk=myBiomodModelEval[,2,,,]
    NAs=which(is.na(jnk))
    all_models=list()
    for (d in dimnames(jnk)[[4]]){
      for (c in dimnames(jnk)[[3]]){
        for (b in dimnames(jnk)[[2]]){
          for (a in dimnames(jnk)[[1]]){
            jnk_str=paste(sp_nm,d,c,b,a,sep="_")
            jnk_str2=paste(sp_nm,d,c,b,sep="_")
            all_models[length(all_models)+1]=jnk_str2
          }}}}
    bad_models_short=all_models[NAs]
    bad_models_short=unique(bad_models_short)
    jnk_good=!(all_models %in% bad_models_short)
    remaining_models=all_models[jnk_good]
    remaining_models=unlist(unique(remaining_models))
    
    
    ###################################################
    ### code chunk number 11: ensemble_modeling
    ###################################################
    myBiomodEM <- BIOMOD_EnsembleModeling(
      modeling.output = myBiomodModelOut,
      chosen.models = remaining_models, #these are not model types (e.g., GBM), but model runs (e.g., PA1_RF)
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
    save("myBiomodEM", "myBiomodModelOut","remaining_models", file=workspace_name_out)   #save workspace
    
    }else{
    cat('\n',sp_nm,'ensemble previously done...')
  }
}
