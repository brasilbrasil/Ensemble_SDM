###USER CONFIGURATION
#see 0_sdm_config file.r

setwd(working_dir)

#loading package libraries
library(biomod2)
library(stringr)
library(colorRamps)
library(rasterVis)
library(tools)

####START UNDERHOOD
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

#sets options for biomod2 fixes in code (if assigned TRUE in config file)
if (apply_biomod2_fixes){
  rasterOptions(tmpdir = dir_for_temp_files, timer = T, progress = "text", todisk  = T) #set options for raster package
  source(paste(codeDir,"Ensemble_SDM/3b_modifications_of_projection_code.r", sep = "/")) #all of fixes to biomod2 code created by AV
}

spp_info = read.csv(paste(csv_dir,'FB_spp_data.csv', sep = "/")) #creates data frame from species info csv file

sp_nm = spp_nm[1] #resets so the first species to run is the first one listed in config file or csv
for (sp_nm in spp_nm){
  sp_nm = as.character(sp_nm) #defines the species name as a character string - not needed if it is already a text name 
  cat('\n',sp_nm,'model projection...') #sign-posting
  workspace_name = paste0(sp_nm,"_FB_EM_fit.RData") #set name of file to load workspace data from model run
  load(workspace_name) #loads workspace from previous ensemble modelling step
  
  plots = paste(working_dir,"AllEMplots_pmw",  sep="/") #assigns location for all plots
  if (file.exists(plots) == F | overwrite == 1){ #run only if plots have not already been collected and the overwrite function is off
    dir.create(plots, showWarnings = FALSE)} #create directory for plots
  
  workspace_name_out = paste0(sp_nm,"_FB_EM_proj_", proj_nm, ".RData") #assigns name to save workspace
  
  if (file.exists(workspace_name_out) == F | overwrite == 1){ #does not run if file already exists and overwrite is off
    #raster_based_env_grid:
    sp_index = which(spp_info[,"Species"] == sp_nm) #assigns index for the line associated with the target species
    raster_res = spp_info[sp_index,"rasterdir"] #assigns the raster resolution to the directory for the species of interest
    cat('\n','using these env files for projection raster:', env_var_files, '\n', 'from dir:', clim_data_dir) #sign-posting
    crop_raster = raster(paste0(crop_raster_dir,"/",raster_res,".grd")) #creates a rasterLayer object from an existing .grd file for ther extent
    predictors = raster(paste(clim_data_dir, env_var_files[1], sep="/")) #creates rasterLayer from the bioclim .grd file
    predictors = crop(predictors,  crop_raster) #crops bioclimate rasterLayer by the raster resolution .grd
    jnk0 = length(env_var_files) #assigns the number of bioclimate variables to a temporary variable for runs
    
    ##creates raster stack from raster layers for all bioclimatic variables of interest
    for (jj in 2:jnk0){ #for all except the first bioclimate variable (which was done above)
      temp = raster(paste(clim_data_dir, env_var_files[jj], sep = "/")) #creates temporary rater layer from bioclimate .grd
      temp = crop(temp,  crop_raster) #crops temporary bioclimate raster layer with the extent raster
      predictors = addLayer(predictors, temp) #adds the new bioclimate raster layer to existing one
    }
    
    var_names <- unlist(file_path_sans_ext(env_var_files)) #creates vector with bioclimatic variable names without the file extension (.grd)
    names(predictors) <- var_names #renames the layers of the raster with the bioclim variable names
    rm("crop_raster" ,"temp", "jnk0") #removes temporary variables
    predictors #returns summary of the predictor rasterStack
    
    ###changing rasterStack for the "wettest" and "driest" scenarios 
    alt_scen = c(2,3,5,6) 
    if (baseline_or_future %in% alt_scen){ 
      predictors<-stack((subset(predictors, 1)),
                        (subset(predictors, 2)),
                        (subset(predictors, 3))*4,#for wet and dry scenarios, multiply bio12 variable (quarterly precipitation) by 4
                        (subset(predictors, 4)))
      names(predictors) <- var_names #assigns names to raster stack again.
    }
    
    if (baseline_or_future == 1){ #for "baseline" runs
      jpeg_name = paste0(sp_nm,"_env_vars_used_for_projection.jpg") #assigns location for map jpeg file
      jpeg(jpeg_name, #settings for map jpeg file
           width = 10, height = 10, units = "in", pointsize = 12, quality = 90, bg = "white", res = 300) 
      plot(predictors, col = rev(terrain.colors(255)), maxpixels = 100000, useRaster = FALSE, axes = TRUE, addfun = NULL, interpolate = TRUE) #builds raster map and sends to jpeg file
      ##RESULTS IN WARNINGS RELATED TO "INTERPOLATE IS NOT A GRAPHICAL PARAMETER" BUT WORK
      dev.off() #turns off jpeg device
    }
    cat('\n',sp_nm,'projection raster stack created...') #sign-posting
    gc() #reclaims memory that is no longer used and returns summary of memory usage
    workspace_name_out0 = paste0(sp_nm,"_FB_all_model_proj_", proj_nm, ".RData") #sets location to save R data
    
    if (file.exists(workspace_name_out0) == F | overwrite == 1){ #does not run if RData file already exists and overwrite is off  
      myBiomodProj_baseline <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut, #results from previous model step
        new.env = predictors, #new environment to project onto (in case of baseline it is not new)
        proj.name = proj_nm, #name for save folder
        selected.models = 'all', #whether all of a subset of the models should be used
        binary.meth = eval_stats, #vector of evaluation statistics to use for projection to presence/ absence
        compress = 'xz', #compression format for files
        build.clamping.mask = clampingMask, #whether a clamping mask should be saved or not
        keep.in.memory = memory) #whether or not the clamping mask should be saved to hard disk
      gc() #reclaims memory that is no longer used.
      cat('\n',sp_nm,'projection complete...') #sign-posting
      cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn') #sign-posting
      save("myBiomodProj_baseline", file = workspace_name_out0)   #save projection workspace      
    }else{
      load(workspace_name_out0) #loads existing workspace
      cat('\n',sp_nm,'projection of individual models loaded from past run...') #sign-posting
    }
    ###################################################
    ### code chunk number 17: projection_current_plot
    ###################################################
    
    if (plot_graphs == T){ #set in the config file
      jnk = length(myBiomodProj_baseline@models.projected) #assigns number of projection models from output
      jnk = jnk/(length(models_to_run)+1) #calculates the number of projections per model type +1
      if (jnk<26){ 
        jpeg_nam = paste0(sp_nm,"_", proj_nm, "_", "runs.jpg") #assigns location for map jpeg
        jpeg(jpeg_name, #settings for the jpeg map
             width = 10, height = 10, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        par(mfrow = c(1,2)) #subsequent figures are drawn in an array of 1 by 2
        try(plot(myBiomodProj_baseline, str.grep = "Full"), TRUE) #attempts to plot results of full models but saves error if not
        dev.off() #turns device off - saves plot as jpeg
      }
    }
    
#     if (plot_graphs == T){ ###this has to be fixed- in situations where there is only enough data for a single PA, no others are run
#       for (model in models_to_run){ #runs for each model
#         jpeg_name = paste0(sp_nm,"_", model, "_BIN_model", proj_nm, "runs.jpg")
#         if (file.exists(jpeg_name) == F | overwrite == 1){ #does not run if file exists and overwrite is off
#           sample=c()
#           sp_bin_file = paste0(proj_nm, "_", sp_nm, "_bin_ROC_RasterStack")
#           sp_bin_dir = paste0(sp_nm,"/proj_", proj_nm, "/", sp_bin_file)
#           if (file.exists(sp_bin_dir)){
#             jnk = load(sp_bin_file) #current_BI_Akepa_bin_ROC_RasterStack
#             sp_bin_stack = get(jnk)
#             try(sample <- raster(sp_bin_stack, layer = paste0(sp_nm, "_AllData_Full_",model,".bin")), TRUE)
#           }
#           sp_bin_file = paste0(proj_nm, "_", sp_nm, "_AllData_Full_", model,"_bin_ROC_RasterLayer")
#           sp_bin_file = paste0(sp_nm,"/proj_", proj_nm, "/", sp_bin_file)          
#           if (file.exists(sp_bin_file)){
#             jnk = load(sp_bin_file) #current_BI_Akepa_bin_ROC_RasterStack
#             sample = get(jnk)
#           }
#           
#           sp_bin_file=paste(proj_nm, "_", sp_nm, "_ROCbin.grd", sep = "")
#           sp_bin_file=paste(sp_nm,"/proj_", proj_nm, "/proj_", sp_bin_file , sep = "") #current_BI_Akepa_bin_ROC_RasterStack
#           if (file.exists(sp_bin_file)){
#             sample=stack(sp_bin_file)
#             names=names(sample)
#             lyr_name=paste0(sp_nm,"_PA",PA.nb.rep,"_Full_", model)
#             jnk=which(names==lyr_name)
#             sample=raster(sample,jnk)
#           }
#           
#           jpeg(jpeg_name,
#                width = 10, height = 10, units = "in",
#                pointsize = 12, quality = 90, bg = "white", res = 300)
#           try(plot(sample), TRUE)
#           dev.off()
#           
#         }
#       }
#     }
    
    ###################################################
    ### code chunk number 18: Individual model plots
    ###################################################
    # plotting the individual pmw modelling approaches in a single graph for all eval-stats --------    
    if (plot_graphs == T){ #set in config file
      for (ii in 1:length(eval_stats)){ #for each evaluation statistic        
        for (i in 1:length(models_to_run)){ #for each model type
          sp_nm = str_replace_all(sp_nm,"_",".") #replaces all "_" with "." to account for biomod2 naming
          binGrdDir <- paste0(working_dir, "/", sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, '_', sp_nm, "_", eval_stats[ii], "bin.grd") #points to location of bin.grd file
          binGrdStack <- stack(binGrdDir) #creates a raster stack from the bin.grd file
          #not sure what next line is attempting
          binGrdSub <- subset(binGrdStack, (length(names(binGrdStack)) - (length(models_to_run)))-1+i) #selects a layer of the rasterStack
          binGrdNames <- names(binGrdSub)
          ###IRC STOPPED HERE###
          a1 <- binGrdSub
          rclassVect <- c(NA, 0)
          rclassMatrix <- matrix(rcl, ncol=2, byrow=TRUE)
          binGrdSub <- reclassify(binGrdSub, rclassMatrix)
          names(binGrdSub) <- binGrdNames
          
          b <- stack(paste(working_dir, sp_nm,"/proj_", proj_nm, "/proj_", proj_nm,"_", sp_nm, ".grd", sep = ""))
          b <- subset(b, (length(names(b)) - (length(models_to_run)))-1+i)
          
          
          b <- subset(b, length(names(b)))
          c <- a*b
          rcl <- c(0, NA)
          rcl <- matrix(rcl, ncol=2, byrow=TRUE)
          c <- reclassify(c, rcl)
          names(c) <- binGrdNames
          
          assign(paste0(eval_stats[ii], '_', models_to_run[i], "_scaled_andbinnedEM_pmw"), c)
          assign(paste0(eval_stats[ii], '_', models_to_run[i], "_binnedEM_pmw"), a1)
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
        jpeg_name=paste(proj_nm,"_", sp_nm,"_All_", eval_stats[ii], "Models_Binandscaled_runs_.jpg", sep = "")
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
    
    # subsetting the S4 class object 'myBiomodProj_baseline'  such that it only uses the two main modelling approaches (GBM and Maxent) 
    # that do not overfit
    
    if (apply_biomod2_fixes){
      myBiomodProjection <- LoadProjectionManually(myBiomodProj_baseline)
    }else{
      myBiomodProjection <- myBiomodProj_baseline
    }
    
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
    #cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn')
    
    ###################################################
    ### code chunk number 19: EnsembleForecasting_loading_res
    ###################################################
    #eval_stats=c("TSS") ###DEBUG###
    #plotting the ensmeble projections per species per projection
    if (plot_graphs==1){
      for (i in 1:length(eval_stats)){
        #setwd(working_dir,)
        jnk=stack(paste(working_dir, sp_nm, "/proj_", proj_nm, "/proj_", proj_nm, "_", sp_nm, "_TotalConsensus_EMby", 
                        eval_stats[i], ".grd", sep = ""))
        a<-subset(jnk, length(names(jnk)))
        nm<-names(a)
        
        jnk=stack(paste(working_dir,sp_nm, "/proj_", proj_nm, "/proj_", proj_nm, "_", sp_nm, "_TotalConsensus_EMby", 
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
      
      jpeg_name=paste(proj_nm,"_", sp_nm,"_TOTALCONSENSUS_Binandscaled_runs_.jpg", sep = "")
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
    }  
    
    setwd(working_dir)
    if (plot_graphs==1){
      for (eval_stat in eval_stats){
        try(load(paste0(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat)), TRUE)    
        try(load(paste0(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.",eval_stat)), TRUE)
        jpeg_name=paste0(sp_nm,"_", eval_stat,"_ensemble_", proj_nm, "runs.jpg")
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        par(mfrow=c(1,2))
        try(plot(get(paste0(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat))), TRUE)
        try(plot(get(paste0(sp_nm,"_AllData_AllRun_EM.",eval_stat))), TRUE)
        sp_bin_file=paste(proj_nm, "_", sp_nm, "_TotalConsensus_EMby", eval_stat,".grd", sep = "")
        sp_bin_file=paste(sp_nm,"/proj_", proj_nm, "/proj_", sp_bin_file , sep = "") #current_BI_Akepa_bin_ROC_RasterStack          
        try(plot(raster(sp_bin_file)), TRUE)
        dev.off()
        eval_stat0=eval_stat
        for (eval_stat in eval_stats){
          try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat0,".bin.",eval_stat, sep = "")), TRUE)    
          try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.",eval_stat0,".bin.",eval_stat, sep = "")), TRUE)
          jpeg_name=paste(sp_nm,"_", eval_stat0,"_ensemble_", proj_nm, "_bin_",eval_stat,"runs.jpg", sep = "")
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
    save("myBiomodProj_baseline", "myBiomodEF", file=workspace_name_out)   #save workspace
    
    removeTmpFiles(h=1)
    cat('\n',sp_nm,'done...')    
  }else{
    cat('\n',sp_nm,'previously calculated...')
  }
}  



