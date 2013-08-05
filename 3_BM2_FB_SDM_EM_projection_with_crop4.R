rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
plot_graphs=1
source(paste0("Y:/PICCC_analysis/code/","directory_registry.r"))
local_config_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m/') #'C:/Users/lfortini/'
#spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Anianiau", "Kauai_Amakihi", "Hawaii_Elepaio", "Palila")
spp_nm=c("Akekee", "Anianiau", "Kauai_Amakihi", "Oahu_Amakihi","Hawaii_Akepa", "Hawaii_Elepaio", "Palila")
baseline_or_future=1 #0 for baseline, 1 for future
overwrite=1 #if 1, will overwrite past results
eval_stats="ROC" #c("ROC", "TSS"), if single, do not place in array 
server=1
models_to_run=c('GBM','MAXENT')
memory = T #keep.in.memory=memory
env_var_files=c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") 

maxentWDtmp = paste("maxentWDtmp_", baseline_or_future, sep = "")
memory.limit(size=24000000)
temp<-paste('Y:/temp/', baseline_or_future, '/', sep='')
dir.create(temp)

if (server==1){
  working_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m_rounded/')
  clim_data_2000=paste0(DR_FB_clim_data,"all_grd/all_baseline/500m_test/")
  clim_data_2100=paste0(DR_FB_clim_data,"all_grd/all_future/500m/")
  clim_data_2000wettest="D:/GIS_Data/REnviroLayers/mixed_data_2000_250mwettest/"
  clim_data_2000driest= "D:/GIS_Data/REnviroLayers/mixed_data_2000_250mdriest/"
  clim_data_2100wettest="D:/GIS_Data/REnviroLayers/mixed_data_2100_250mwettest/"
  clim_data_2100driest= "D:/GIS_Data/REnviroLayers/mixed_data_2100_250mdriest/"  
}else{
  working_dir='C:/Users/lfortini/Forest bird SDM/biomod2/'
  clim_data_2000="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/250m/"
  clim_data_2100="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_future/250m/"
}
####START UNDERHOOD
crop_raster_dir=paste(working_dir, 'map_crop/',sep="")

if (baseline_or_future==1){
  clim_surface_to_use=clim_data_2000 
  proj_nm0='baseline'}
if (baseline_or_future==2){
  clim_surface_to_use=clim_data_2000wettest
  proj_nm0='baseline_wettest'}
if (baseline_or_future==3){
  clim_surface_to_use=clim_data_2000driest
  proj_nm0='baseline_driest'}
if (baseline_or_future==4){
  clim_surface_to_use=clim_data_2100 
  proj_nm0='future'}
if (baseline_or_future==5){
  clim_surface_to_use=clim_data_2100wettest
  proj_nm0='future_wettest'}
if (baseline_or_future==6){
  clim_surface_to_use=clim_data_2100driest
  proj_nm0='future_driest'}
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")

setwd(working_dir)
spp_nm0=spp_nm
eval_stats0=eval_stats
library(biomod2)
library(stringr)
library(colorRamps)
library(rasterVis)
rasterOptions(tmpdir=temp, timer = T, progress = "text", todisk  = T)

#sp_nm="Akepa" #debug
var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}
memory.limit(size=240000)
#sp_nm=spp_nm[1]
spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))

source(paste0(DR_code_S,"Ensemble_SDM/3b_modifications_of_projection_code.r")) #all of fixes to biomod2 code created by AV

for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  workspace_name=paste(sp_nm0,"_FB_EM_fit.RData", sep = "") #set name of file to load workspace data from model run
  load(workspace_name)
  
  plots=paste(working_dir,"/AllEMplots_pmw/",  sep="")
  if (file.exists(plots)==F | overwrite==1){
    dir.create(plots)}
  
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
    cat('\n','using these env files for projection raster:', env_var_files, '\n', 'from dir:', clim_data_dir)
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
    
    if (baseline_or_future==2|3|5|6){
      predictors<-stack((subset(predictors, 1)),
                        (subset(predictors, 2)),
                        (subset(predictors, 3))*4,
                        (subset(predictors, 4)))
      names(predictors)<- var_name
    }
    cat('\n',sp_nm,'projection raster stack created...')
    gc()
    workspace_name_out0=paste(sp_nm,"_FB_all_model_proj_", proj_nm, ".RData", sep = "")
    if (file.exists(workspace_name_out0)==F | overwrite==1){  
      myBiomomodProj_baseline <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = stack(predictors), #error: additional stack fx
        proj.name = proj_nm,
        selected.models = 'all',
        binary.meth = eval_stats,
        compress = 'xz',
        clamping.mask = F, 
        keep.in.memory=memory)
      gc()
      cat('\n',sp_nm,'projection complete...')
      cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn')
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
    #     if (plot_graphs==1){
    #       for (model in models_to_run){
    #         jpeg_name=paste(sp_nm0,"_", model, "_BIN_model", proj_nm, "runs.jpg", sep = "")
    #         if (file.exists(jpeg_name)==F | overwrite==1){
    #           jnk=paste(working_dir, sp_nm, "/proj_", proj_nm, "/", proj_nm, "_", sp_nm, "_bin_ROC_RasterStack" , sep = "")
    #           if (file.exists(jnk)){
    #             sp_bin_file=paste(proj_nm, "_", sp_nm, "_bin_ROC_RasterStack" , sep = "")
    #             jnk=load(paste(sp_nm,"/proj_", proj_nm, "/", sp_bin_file , sep = "")) #current_BI_Akepa_bin_ROC_RasterStack
    #             sp_bin_stack=get(jnk)
    #             sample=c()
    #             try(sample <- raster(sp_bin_stack, layer=paste(sp_nm, "_AllData_Full_",model,".bin" , sep = "")), TRUE)
    #           }else{
    #             sp_bin_file=paste(proj_nm, "_", sp_nm, "_AllData_Full_", model,"_bin_ROC_RasterLayer" , sep = "")
    #             jnk=load(paste(sp_nm,"/proj_", proj_nm, "/", sp_bin_file , sep = "")) #current_BI_Akepa_bin_ROC_RasterStack
    #             sample=get(jnk)
    #           }
    #           
    #           jpeg(jpeg_name,
    #                width = 10, height = 10, units = "in",
    #                pointsize = 12, quality = 90, bg = "white", res = 300)
    #           #par(mfrow=c(1,2))
    #           try(plot(sample), TRUE)
    #           dev.off()
    #         }
    #       }
    #     }
    
    ###################################################
    ### code chunk number 18: Individual model plots
    ###################################################
    # plotting the individual pmw modelling approaches in a single graph for all eval-stats --------    
    if (plot_graphs==1){
      for (ii in 1:length(eval_stats)){         
        for (i in 1:length(models_to_run)){
          #setwd(working_dir)
          
          a <- stack(paste(working_dir, sp_nm0,"/proj_", proj_nm, "/proj_", proj_nm, '_', sp_nm0, "_", eval_stats[ii], "bin.grd", sep = ""))
          
          a<-subset(a, (length(names(a)) - (length(models_to_run)))-1+i)#:length(names(a)))
          nm<-names(a)
          a1<-a
          rcl <- c(NA, 0)
          rcl <- matrix(rcl, ncol=2, byrow=TRUE)
          a <-reclassify(a, rcl)
          names(a)<-nm
          
          b <- stack(paste(working_dir, sp_nm0,"/proj_", proj_nm, "/proj_", proj_nm,"_", sp_nm0, ".grd", sep = ""))
          b<-subset(b, (length(names(b)) - (length(models_to_run)))-1+i)
          
          
          b<-subset(b, length(names(b)))
          c<-a*b
          rcl <- c(0, NA)
          rcl <- matrix(rcl, ncol=2, byrow=TRUE)
          c<-reclassify(c, rcl)
          names(c)<-nm
          
          assign(paste(eval_stats[ii], '_', models_to_run[i], "_scaled_andbinnedEM_pmw", sep = ""), c)
          assign(paste(eval_stats[ii], '_', models_to_run[i], "_binnedEM_pmw", sep = ""), a1)
        }
        
        
        c2<-stack(get(paste(eval_stats[ii], '_', models_to_run[1], "_scaled_andbinnedEM_pmw", sep = "")))
        for (i in 2:length(models_to_run)){
          c2<-addLayer(c2, get(paste(eval_stats[ii], '_', models_to_run[i], "_scaled_andbinnedEM_pmw", sep = "")))
        }
        
        a2<-stack(get(paste(eval_stats[ii], '_', models_to_run[1], "_binnedEM_pmw", sep= "")))
        for (i in 2:length(models_to_run)){
          a2<-addLayer(a2, get(paste(eval_stats[ii], '_', models_to_run[i], "_binnedEM_pmw", sep= "")))
        }
        
        setwd(plots)
        jpeg_name=paste(proj_nm,"_", sp_nm0,"_All_", eval_stats[ii], "Models_Binandscaled_runs_.jpg", sep = "")
        jpeg(jpeg_name, width = 5*length(models_to_run), height = 5, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)  
        
        par(pin = c(4,4), cex = 1, cex.main = 1, cex.axis = 0.8, mfcol=c(1,length(models_to_run)), mgp = c(1, 0.5, 0),
            mar=c(2, 2, 1.5, 0), oma = c(0, 0, 0, 1), bg = "transparent")
        
        gc = c('antiquewhite1', 'transparent')
        #grd <- terrain.colors(255)
        col5 <- colorRampPalette(c('blue', 'sandybrown', 'darkgreen'))
        
        jnk <- subset(a2, 1)
        try(plot(c2, 1,  col = col5(255), useRaster=FALSE, axes = TRUE, addfun=F, interpolate = TRUE, legend = F, add = F, bg = "transparent"),silent=T)
        plot(jnk, col = gc, useRaster=FALSE, axes = F, addfun=F, interpolate = TRUE, legend = F, add = T)
        
        par(mar=c(2, 0, 1.5, 3)) # comment this out if there are three models_to_run
        jnk <- subset(a2, 2)
        try(plot(c2, 2, col = col5(255), useRaster=FALSE, axes = TRUE,  interpolate = TRUE, legend = T, yaxt = 'n', add = F, bg = "transparent"),silent=T)  #addfun=F,
        plot(jnk, col = gc, useRaster=FALSE, axes = F, interpolate = F,  legend = F, add = T)#plot(jnk, col = gc, useRaster=FALSE, axes = F, addfun=F, interpolate = TRUE, yaxt = 'n', legend = F, add = T)
        
        #         par(mar=c(2, 0, 1.5, 3))
        #         jnk <- subset(a2, 3)
        #         plot(c2, 3, col = col5(255), useRaster=FALSE, axes = T, interpolate = F, legend = T, yaxt = 'n', add = F, bg = "transparent")
        #         plot(jnk, col = gc, useRaster=FALSE, axes = F, interpolate = F, legend = F, add = T)
        #         
        legend("bottomright",legend = c("Absent"), fill = gc[1], cex = 0.8)
        dev.off()
      }
    }
    setwd(working_dir)
    cat('/n',sp_nm,'done with individual model plots...')
    
    # subsetting the S4 class object 'myBiomomodProj_baseline'  such that it only uses the two main modelling approaches (GBM and Maxent) 
    # that do not overfit
    
    myBiomodProjection <- LoadProjectionManually(myBiomomodProj_baseline)
    
    cat('\n',sp_nm,'projection graphs done...')
    
    ###################################################
    ### code chunk number 18: EnsembleForecasting_future
    ###################################################
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      projection.output = myBiomodProjection,
      total.consensus = T,
      EM.output = myBiomodEM, 
      binary.meth=eval_stats, 
      keep.in.memory=memory)
    cat('\n',sp_nm,'ensemble projection done...')
    cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn')
    
    ###################################################
    ### code chunk number 19: EnsembleForecasting_loading_res
    ###################################################
    #eval_stats=c("TSS") ###DEBUG###
    #plotting the ensmeble projections per species per projection
    if (plot_graphs==1){
      for (i in 1:length(eval_stats)){
        #setwd(working_dir,)
        jnk=stack(paste(working_dir, sp_nm0, "/proj_", proj_nm, "/proj_", proj_nm, "_", sp_nm0, "_TotalConsensus_EMby", 
                        eval_stats[i], ".grd", sep = ""))
        a<-subset(jnk, length(names(jnk)))
        nm<-names(a)
        
        jnk=stack(paste(working_dir,sp_nm0, "/proj_", proj_nm, "/proj_", proj_nm, "_", sp_nm0, "_TotalConsensus_EMby", 
                        eval_stats[i], "_", eval_stats[i], "bin.grd" , sep = ""))
        b<-subset(jnk, length(names(jnk)))
        
        c<-a*b
        rcl <- c(0, NA)
        rcl <- matrix(rcl, ncol=2, byrow=TRUE)
        c<-reclassify(c, rcl)
        names(c)<-nm
        assign(paste("TotalConsensus_EMScaledandBinnedby_", eval_stats[i], sep=""), c)
        assign(paste("TotalConsensus_EMBinnedby_", eval_stats[i], sep=""), b)
      }
    }
    
    ems_a<-stack(get(paste("TotalConsensus_EMScaledandBinnedby_", eval_stats[1], sep="")))
    if (length(eval_stats)>1){ 
    for (i in 2:length(eval_stats)){
      ems_a<-addLayer(ems_a, get(paste("TotalConsensus_EMScaledandBinnedby_", eval_stats[i], sep="")))
    }}
    
    ems_b<-stack(get(paste("TotalConsensus_EMBinnedby_", eval_stats[1], sep="")))
    if (length(eval_stats)>1){ 
      for (i in 2:length(eval_stats)){
      ems_b<-addLayer(ems_b, get(paste("TotalConsensus_EMBinnedby_", eval_stats[i], sep="")))
    }}
    
    setwd(plots)
    
    jpeg_name=paste(proj_nm,"_", sp_nm0,"_TOTALCONSENSUS_Binandscaled_runs_.jpg", sep = "")
    jpeg(jpeg_name, width = 5*length(eval_stats), height = 5, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)  
    par(pin = c(4,4), cex = 1, cex.main = 1, cex.axis = 0.8, mfcol=c(1,length(eval_stats)), mgp = c(1, 0.5, 0),
        mar=c(2, 2, 1.5, 1), oma = c(0, 0, 0, 1), bg = "transparent")
    
    gc = c('antiquewhite1', 'transparent')
    col5 <- colorRampPalette(c('blue', 'sandybrown', 'darkgreen'))
    jnk <- subset(ems_b, 1)
    try(plot(ems_a[[1]],  col = col5(255), useRaster=FALSE, axes = TRUE, addfun=F, interpolate = TRUE, legend = F, add = F, bg = "transparent"),silent=T)
    plot(jnk, col = gc, useRaster=FALSE, axes = F, addfun=F, interpolate = TRUE, legend = F, add = T)
    
    if (length(eval_stats)>1){ 
      par(mar=c(2, 0, 1.5, 0))
    jnk <- subset(ems_b, 2)
    try(plot(ems_a[[2]], col = col5(255), useRaster=FALSE, axes = TRUE,interpolate = TRUE, legend = F, yaxt = 'n', add = F, bg = "transparent"),silent=T)  # addfun=F, 
    plot(jnk, col = gc, useRaster=FALSE, axes = F, addfun=F, interpolate = TRUE, yaxt = 'n', legend = F, add = T)
    }
    
    if (length(eval_stats)>2){ 
      par(mar=c(2, 0, 1.5, 3.5))
    jnk <- subset(ems_b, 3)
    try(plot(ems_a[[3]], col = col5(255), useRaster=FALSE, axes = T, interpolate = F, legend = T, yaxt = 'n', add = F, bg = "transparent"),silent=T)
    plot(jnk, col = gc, useRaster=FALSE, axes = F, interpolate = F, legend = F, add = T)
    }
    
    legend("bottomright",legend = c("Absent"), fill = gc[1], cex = 0.8)
    dev.off()
    setwd(working_dir)
    if (plot_graphs==1){
      for (eval_stat in eval_stats){
        try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat, sep = "")), TRUE)    
        try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.",eval_stat, sep = "")), TRUE)
        #try(load(paste(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_", sp_nm,"_TotalConsensus_EMby",eval_stat, ".grd", sep = "")), TRUE)
        jpeg_name=paste(sp_nm0,"_", eval_stat,"_ensemble_", proj_nm, "runs.jpg", sep = "")
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        par(mfrow=c(1,2))
        try(plot(get(paste(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat, sep = ""))), TRUE)
        try(plot(get(paste(sp_nm,"_AllData_AllRun_EM.",eval_stat, sep = ""))), TRUE)
        try(plot(raster(paste0(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_", sp_nm,"_TotalConsensus_EMby",eval_stat, ".grd"))), TRUE)
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
          try(plot(raster(paste0(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_", sp_nm,"_TotalConsensus_EMby",eval_stat,"_", eval_stat, "bin.grd"))), TRUE)
          dev.off()
        }
      }
    }
    cat('\n',sp_nm,'ensemble projection figures done...')
    save.image("temp_workspace4.RData")   #to save workspace
    rm(list=c("spp_info","eval_stats0", "spp_nm0", "crop_raster", "clim_surface_to_use", "proj_nm0", "overwrite", 
              "plot_graphs", "local_config_dir","spp_nm", "clim_data_2000", "predictors_temp",
              "clim_data_2100", "working_dir", "env_var_files", "csv_dir", "eval_stats"))      
    save.image(workspace_name_out)   #to save workspace
    load("temp_workspace4.RData")
    removeTmpFiles(h=1)
    cat('\n',sp_nm,'done...')    
    #save.image(workspace_name)
  }else{
    cat('\n',sp_nm,'previously calculated...')
  }
}  



