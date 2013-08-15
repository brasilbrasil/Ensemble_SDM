rm(list = ls()) #remove all past worksheet variables
source(paste0("Y:/PICCC_analysis/code/","directory_registry.r"))

###USER CONFIGURATION
plot_graphs=1
#local_config_dir='C:/Users/lfortini/'
#spp_nm=(read.csv(paste(local_config_dir,'spp_to_run.csv', sep = ""),header=F, stringsAsFactors=F))
spp_nm=c("Akekee", "Palila", "Hawaii_Akepa")#"Kauai_Amakihi", "Anianiau", "Apapane", "Iiwi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi")   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
clim_data_2000=paste0(DR_FB_clim_data,"all_grd/all_baseline/250m/")
clim_data_2100=paste0(DR_FB_clim_data,"all_grd/all_future/500m/")
clim_surface_to_use=clim_data_2000 
proj_nm0='baseline' 
models_to_run=c('GBM','RF','MAXENT')
overwrite=0 #if 1, will overwrite past results

project_name='test_runs_old_code_new_package'
working_dir=paste0(DR_FB_SDM_results_S,project_name,'/')
env_var_files=c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") 
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")
crop_raster_dir=paste(working_dir, 'map_crop/',sep="")
eval_stats=c("ROC") 

####START UNDERHOOD
setwd(working_dir)
spp_nm0=spp_nm
eval_stats0=eval_stats
library(biomod2)
library(stringr)
#sp_nm="Akepa" #debug
var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}
memory.limit(size=4095)
sp_nm=spp_nm[1]
spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))

for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  workspace_name=paste(sp_nm,"_FB_EM_fit.RData", sep = "") #set name of file to load workspace data from model run
  #if (file.exists(workspace_name)==F){ 
  #  workspace_name=paste(sp_nm,"_FB_run.RData", sep = "")}
  load(workspace_name)
  
  #model run specific variables that must not be saved to workspace
  spp_nm=spp_nm0
  eval_stats=eval_stats0  
  clim_data_dir0=clim_surface_to_use
  proj_nm=proj_nm0 
  
  sp_nm=str_replace_all(sp_nm,"_", ".")
  workspace_name_out=paste(sp_nm,"_FB_EM_proj_", proj_nm, ".RData", sep = "")
  #jnk=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.TSS", sep = "")
  #jnk2=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.TSS", sep = "")
  
  if (file.exists(workspace_name_out)==F | overwrite==1){
    #raster_based_env_grid:
    sp_index=which(spp_info[,"Species"]==sp_nm0)
    raster_res= spp_info[sp_index,"rasterdir"]
    clim_data_dir=clim_data_dir0 
    jnk0=length(env_var_files)
    crop_raster=raster(paste(crop_raster_dir,raster_res,".grd",sep=""))
    predictors = raster( paste(clim_data_dir, env_var_files[1], sep=""))
    predictors=crop(predictors,  crop_raster)
    for (jj in 2:jnk0){
      temp=raster(paste(clim_data_dir, env_var_files[jj], sep=""))
      temp=crop(temp,  crop_raster)
      predictors = addLayer(predictors, temp)
    }
    names(predictors)<- var_name
    rm("crop_raster" ,"temp") 
    predictors
    cat('\n',sp_nm,'projection raster stack created...')
    
    workspace_name_out0=paste(sp_nm,"_FB_all_model_proj_", proj_nm, ".RData", sep = "")
    if (file.exists(workspace_name_out0)==F | overwrite==1){  
      myBiomomodProj_baseline <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = predictors, #error: additional stack fx
        proj.name = proj_nm,
        selected.models = 'all',
        binary.meth = eval_stats,
        compress = 'xz',
        clamping.mask = F)
      
      cat('\n',sp_nm,'projection complete...')
      save.image("temp_workspace3.RData")   #to save workspace
      rm(list=c("spp_info", "eval_stats0", "spp_nm0", "clim_surface_to_use", "proj_nm0", "overwrite", 
                "plot_graphs", "local_config_dir","spp_nm", "clim_data_2000", 
                "clim_data_2100", "working_dir", "env_var_files", "csv_dir", "eval_stats"))      
      save.image(workspace_name_out0)   #to save workspace
      load("temp_workspace3.RData")
    }else{
      load(workspace_name_out0)
      cat('\n',sp_nm,'projection of individual models loaded from past run...')
    }
    ###################################################
    ### code chunk number 17: projection_current_plot
    ###################################################
    # make some plots, sub-selected by str.grep argument
    #plot(myBiomomodProj_baseline, str.grep = 'GBM')
    if (plot_graphs==1){
      for (model in models_to_run){
        jpeg_name=paste(sp_nm0,"_", model, "_", proj_nm, "runs.jpg", sep = "")
        jpeg(jpeg_name,
             width = 10, height = 10, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        par(mfrow=c(1,2))
        try(plot(myBiomomodProj_baseline, str.grep = model), TRUE)
        dev.off()
      }
    }
    
    #print sample binary model for species
    if (plot_graphs==1){
      for (model in models_to_run){
        jpeg_name=paste(sp_nm0,"_", model, "_BIN_model", proj_nm, "runs.jpg", sep = "")
        if (file.exists(jpeg_name)==F | overwrite==1){
          sp_bin_file=paste(proj_nm, "_", sp_nm, "_bin_ROC_RasterStack" , sep = "")
          sp_bin_file=paste(sp_nm,"/proj_", proj_nm, "/", sp_bin_file , sep = "")
          if (file.exists(sp_bin_file)){
            jnk=load(sp_bin_file) #current_BI_Akepa_bin_ROC_RasterStack
            sp_bin_stack=get(jnk)
            sample=c()
            try(sample <- raster(sp_bin_stack, layer=paste(sp_nm, "_AllData_Full_",model,".bin" , sep = "")), TRUE)
          }
          sp_bin_file=paste(proj_nm, "_", sp_nm, "_AllData_Full_", model,"_bin_ROC_RasterLayer" , sep = "")
          sp_bin_file=paste(sp_nm,"/proj_", proj_nm, "/", sp_bin_file , sep = "")          
          if (file.exists(sp_bin_file)){
            jnk=load(sp_bin_file) #current_BI_Akepa_bin_ROC_RasterStack
            sample=get(jnk)
          }
          
          sp_bin_file=paste(proj_nm, "_", sp_nm, "_ROCbin.grd", sep = "")
          sp_bin_file=paste(sp_nm,"/proj_", proj_nm, "/proj_", sp_bin_file , sep = "") #current_BI_Akepa_bin_ROC_RasterStack
          if (file.exists(sp_bin_file)){
            sample=raster(sp_bin_file)
            #plot(jnk)  
          }
          
          jpeg(jpeg_name,
               width = 10, height = 10, units = "in",
               pointsize = 12, quality = 90, bg = "white", res = 300)
          #par(mfrow=c(1,2))
          try(plot(sample), TRUE)
          dev.off()
        }
      }
    }
    
    cat('\n',sp_nm,'projection graphs done...')
    
    
    ###################################################
    ### code chunk number 18: EnsembleForecasting_future
    ###################################################
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      projection.output = myBiomomodProj_baseline,
      EM.output = myBiomodEM, binary.meth=eval_stats) ###DEBUG### , binary.meth=c('TSS', 'KAPPA', 'ROC')
    cat('\n',sp_nm,'ensemble projection done...')
    
    ###################################################
    ### code chunk number 19: EnsembleForecasting_loading_res
    ###################################################
    #eval_stats=c("TSS") ###DEBUG###
    if (plot_graphs==1){
      for (eval_stat in eval_stats){
        try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat, sep = "")), TRUE)    
        try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.",eval_stat, sep = "")), TRUE)
        jpeg_name=paste(sp_nm0,"_", eval_stat,"_ensemble_", proj_nm, "runs.jpg", sep = "")
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        par(mfrow=c(1,2))
        try(plot(get(paste(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat, sep = ""))), TRUE)
        try(plot(get(paste(sp_nm,"_AllData_AllRun_EM.",eval_stat, sep = ""))), TRUE)
        sp_bin_file=paste(proj_nm, "_", sp_nm, "_TotalConsensus_EMby", eval_stat,".grd", sep = "")
        sp_bin_file=paste(sp_nm,"/proj_", proj_nm, "/proj_", sp_bin_file , sep = "") #current_BI_Akepa_bin_ROC_RasterStack          
        try(plot(raster(sp_bin_file)), TRUE)
        dev.off()
        eval_stat0=eval_stat
        for (eval_stat in eval_stats){
          try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat0,".bin.",eval_stat, sep = "")), TRUE)    
          try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.",eval_stat0,".bin.",eval_stat, sep = "")), TRUE)
          jpeg_name=paste(sp_nm0,"_", eval_stat0,"_ensemble_", proj_nm, "_bin_",eval_stat,"runs.jpg", sep = "")
          jpeg(jpeg_name,
               width = 10, height = 8, units = "in",
               pointsize = 12, quality = 90, bg = "white", res = 300)
          #par(mfrow=c(1,2))
          try(plot(get(paste(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat0,".bin.",eval_stat, sep = ""))), TRUE)
          try(plot(get(paste(sp_nm,"_AllData_AllRun_EM.",eval_stat0,".bin.",eval_stat, sep = ""))), TRUE)
          sp_bin_file=paste(proj_nm, "_", sp_nm, "_TotalConsensus_EMby", eval_stat0,"_",eval_stat, "bin.grd", sep = "")
          sp_bin_file=paste(sp_nm,"/proj_", proj_nm, "/proj_", sp_bin_file , sep = "") #current_BI_Akepa_bin_ROC_RasterStack          
          try(plot(raster(sp_bin_file)), TRUE)
          dev.off()
        }
      }
    }
    cat('\n',sp_nm,'ensemble projection figures done...')
    cat('\n',sp_nm,'projection graphs done...')
    save.image("temp_workspace4.RData")   #to save workspace
    rm(list=c("spp_info","eval_stats0", "spp_nm0", "clim_surface_to_use", "proj_nm0", "overwrite", 
              "plot_graphs", "local_config_dir","spp_nm", "clim_data_2000", 
              "clim_data_2100", "working_dir", "env_var_files", "csv_dir", "eval_stats"))      
    save.image(workspace_name_out)   #to save workspace
    load("temp_workspace4.RData")
    
    #save.image(workspace_name)
  }else{
    cat('\n',sp_nm,'previously calculated...')
  }
}

