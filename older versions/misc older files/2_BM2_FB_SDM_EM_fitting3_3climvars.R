rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
server=1
plot_graphs=1
local_config_dir='C:/Users/lfortini/'
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Akekee")#, "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi")#, "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi")   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
if(server==1){
  working_dir='Y:/FB SDM/biomod2/'
}else{
  working_dir='C:/Users/lfortini/Forest bird SDM/biomod2/'

}
#env_var_files=c("slope.grd", "trasp.grd", "tri.grd", "bio15.grd", "bio16.grd", "bio2.grd", "bio3.grd", "bio4.grd", "bio5.grd", "bio6.grd") ###DEBUG!!
#env_var_files=c("Tmax.grd", "Tmin.grd", "ppt.grd", "slope.grd", "trasp.grd", "tri.grd") ###DEBUG!!
env_var_files=c("Tmin.tif", "ppt.tif", "bio4.tif", "bio16.tif")

csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")
overwrite=1

####START UNDERHOOD
setwd(working_dir)

library(biomod2)
library(stringr)
#sp_nm="Akepa" #debug
var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}
#memory.limit(size=24000000)

#sp_nm=spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  workspace_name=paste(sp_nm,"_FB_run.RData", sep = "") #set name of file to load workspace data from model run
  workspace_name_out=paste(sp_nm,"_FB_EM_fit.RData", sep = "") #set name of file to load workspace data from model run
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
      eval.metric = c('ROC'), #c('TSS', 'ROC', 'KAPPA'); 'all', #c('TSS', 'ROC'),
      eval.metric.quality.threshold = c(0.5),
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
    
    #load(workspace_name_out)
    save.image("temp_workspace2.RData")   #to save workspace
    rm(list=c("spp_info","sp_nm","local_config_dir", "spp_nm", "models_to_run", "working_dir", 
              "clim_data_dir0", "env_var_files", "csv_dir", "spp_info", "var_name",
              "eval_stats0", "spp_nm0", "clim_surface_to_use", "proj_nm0", "overwrite", 
              "plot_graphs", "local_config_dir","spp_nm", "clim_data_2000", 
              "clim_data_2100", "working_dir", "env_var_files", "csv_dir", "eval_stats"))      
    save.image(workspace_name_out)   #to save workspace
    load("temp_workspace2.RData")
  }else{
    cat('\n',sp_nm,'ensemble previously done...')
  }
}
