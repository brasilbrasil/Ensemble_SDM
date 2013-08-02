rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
server=1
plot_graphs=1
source(paste0("Y:/PICCC_analysis/code/","directory_registry.r"))
local_config_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m/') #'C:/Users/lfortini/'
#spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
spp_nm=c("Akekee")#, "Anianiau", "Kauai_Amakihi", "Oahu_Amakihi","Hawaii_Akepa", "Hawaii_Elepaio", "Palila")
models_to_run=c('GBM','MAXENT')
eval_stats=c("ROC")

if(server==1){
  working_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m_rounded/')
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/test/'
}

csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")
overwrite=0

####START UNDERHOOD
setwd(working_dir)

library(biomod2)
library(stringr)
#sp_nm="Akepa" #debug

memory.limit(size=4095)
sp_nm=spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  workspace_name=paste(sp_nm0,"_FB_run.RData", sep = "") #set name of file to load workspace data from model run
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
    ### response curves
    ###################################################
    
    myGLMs <- BIOMOD_LoadModels(myBiomodModelOut, models="GBM")
    
    # 4.2 plot 2D response plots
    myRespPlot2D <- response.plot2(models  = myGLMs,
                                   Data = getModelsInputData(myBiomodModelOut,'expl.var'), 
                                   show.variables= getModelsInputData(myBiomodModelOut,'expl.var.names'),
                                   do.bivariate = FALSE,
                                   fixed.var.metric = 'mean',
                                   save.file="no", 
                                   name="response_curve", 
                                   ImageSize=480, 
                                   plot=TRUE)
    
    # 4.2 plot 3D response plots
    ## here only for a lone model (i.e "MyocastorCoypus_PA1_RUN1_GLM")
#     myRespPlot3D <- response.plot2(models  = myGLMs[1],
#                                    Data = getModelsInputData(myBiomodModelOut,'expl.var'), 
#                                    show.variables= getModelsInputData(myBiomodModelOut,'expl.var.names'),
#                                    do.bivariate = TRUE,
#                                    fixed.var.metric = 'mean',
#                                    save.file="no", 
#                                    name="response_curve", 
#                                    ImageSize=480, 
#                                    plot=TRUE)
    
    ### all the values used to produce this plot are stored into the returned object
    ### you can redo plots by yourself and customised them
    dim(myRespPlot2D)
    dimnames(myRespPlot2D)
    
    dev.off()
    for (asd in 1:dim(myRespPlot2D)[3]){
      ymax_lim=max(myRespPlot2D[,2,asd,])
      xmax_lim=max(myRespPlot2D[,1,asd,])
      xmin_lim=max(myRespPlot2D[,2,asd,])
      
    }
    
    var=myRespPlot2D[,1,1,1]
    pred=myRespPlot2D[,2,1,1]
    plot(var,pred, type="l")
    
    #dim(myRespPlot3D)
    #dimnames(myRespPlot3D)
    
    
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
      eval.metric = eval_stats, #c('TSS', 'ROC', 'KAPPA'); 'all', #c('TSS', 'ROC'),
      eval.metric.quality.threshold = rep(0.15,length(eval_stats)),
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
              "clim_data_dir0", "csv_dir", "spp_info", "var_name",
              "eval_stats0", "spp_nm0", "clim_surface_to_use", "proj_nm0", "overwrite", 
              "plot_graphs", "local_config_dir","spp_nm", "clim_data_2000", 
              "clim_data_2100", "working_dir", "env_var_files", "csv_dir", "eval_stats"))      
    save.image(workspace_name_out)   #to save workspace
    load("temp_workspace2.RData")
  }else{
    cat('\n',sp_nm,'ensemble previously done...')
  }
}
