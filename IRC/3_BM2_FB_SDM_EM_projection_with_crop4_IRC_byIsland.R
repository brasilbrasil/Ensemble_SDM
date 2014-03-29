###USER CONFIGURATION
#see 0_sdm_config file.r

setwd(working_dir)

#loading package libraries
library(biomod2)
library(stringr)
library(colorRamps)
library(rasterVis)
library(tools)
library(ncdf)
library(gbm) #needed for projection of "gbm" data (i.e. "BIOMOD_Projection(...)" function)

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
  source(paste(codeDir,"3b_modifications_of_projection_code.r", sep = "/")) #all of fixes to biomod2 code created by AV
}

spp_info = read.csv(paste(csv_dir,'FB_spp_data.csv', sep = "/")) #creates data frame from species info csv file

sp_nm = spp_nm[3] #resets so the first species to run is the first one listed in config file or csv
for (sp_nm in spp_nm){
  sp_nm = as.character(sp_nm) #defines the species name as a character string - not needed if it is already a text name 
  sp_dir = str_replace_all(sp_nm,"_", ".") #replaces "_" with "." in sp_nm
  cat('\n',sp_nm,'model projection...') #sign-posting
  workspace_name = paste0(sp_nm,"_FB_EM_fit.RData") #set name of file to load workspace data from model run
  load(workspace_name) #loads workspace from previous ensemble modelling step
  plots = paste0(working_dir, "AllEMplots_pmw") #assigns location for all plots
  
  if (file.exists(plots) == F | overwriteData == T) { #run only if plots have not already been collected and the overwrite function is off
    dir.create(plots, showWarnings = FALSE)
  } #create directory for plots
  
  workspace_name_out = paste0(sp_nm,"_FB_EM_proj_", proj_nm, ".RData") #assigns name to save workspace
  
  if (file.exists(workspace_name_out) == F | overwriteData == T){ #does not run if file already exists and overwrite is off
    #raster_based_env_grid:
    sp_index = which(spp_info[,"Species"] == sp_nm) #assigns index for the line associated with the target species
    raster_res = spp_info[sp_index, "rasterdir"] #assigns the raster resolution to the directory for the species of interest
    cat('\n','using these env files for projection raster:', env_var_files, 
        '\n', 'from dir:', clim_data_dir) #sign-posting
    crop_raster = raster(paste0(crop_raster_dir,raster_res,".grd")) #creates a rasterLayer object from an existing .grd file for ther extent
    predictors = raster(paste0(clim_data_dir, env_var_files[1])) #creates rasterLayer from the bioclim .grd file
    predictors = crop(predictors,  crop_raster) #crops bioclimate rasterLayer by the raster resolution .grd
    jnk0 = length(env_var_files) #assigns the number of bioclimate variables to a temporary variable for runs
    
    ##creates raster stack from raster layers for all bioclimatic variables of interest
    for (jj in 2:jnk0){ #for all except the first bioclimate variable (which was done above)
      temp = raster(paste0(clim_data_dir, env_var_files[jj])) #creates temporary rater layer from bioclimate .grd
      temp = crop(temp,  crop_raster) #crops temporary bioclimate raster layer with the extent raster
      predictors = addLayer(predictors, temp) #adds the new bioclimate raster layer to existing one
    }
    
    var_names <- unlist(file_path_sans_ext(env_var_files)) #creates vector with bioclimatic variable names without the file extension (.grd)
    names(predictors) <- var_names #renames the layers of the raster with the bioclim variable names
    rm("crop_raster" ,"temp", "jnk0") #removes temporary variables
    predictors #returns summary of the predictor rasterStack

    
    ###island by island code
    
    # Defining the extent of the different islands
    Kauai = c(-159.82,-159.26, 21.84, 22.25)
    Oahu = c(-158.32, -157.62,  21.22, 21.73)
    Molokai = c(-157.34, -156.69, 21.03, 21.25)
    Lanai = c(-157.08, -156.78, 20.70, 20.92)
    Maui= c(-156.8, -155.53, 20.46, 21.05)
    Hawaii = c(-156.10,-154.74, 18.87, 20.30)
    Kahoolawe = c(-156.8, -156.51, 20.46, 20.62)
    
    #Identify which islands the species is found on
    sp_row <- which(spp_info[,"Species"] == sp_nm) #returns the row number for the species of interest - same as sp_index above
    spIslands <- spp_info[sp_row,(6:length(names(spp_info)))] #returns dataframe indicating whether or not the species is found on each of 6 main Hawaiian islands 
    spIslandNames <- names(spIslands)[spIslands > 0] #names of islands where the species is found
    
    #add Kahoolawe to species island list if Maui is included - so that it can later be cut out
    if ("Maui" %in% spIslandNames) {
      spIslandNames <- append(spIslandNames, "Kahoolawe")
    }

    # Cutting out each island
    for (i in 1:length(spIslandNames)){
      e = extent(get(spIslandNames[i]))
      Isras = stack(crop(predictors, e, snap = 'in'))
      names(Isras)<- var_names
      assign(spIslandNames[i], Isras)
    }
    
    # for Maui need to cut out Kahoolawe, but because of the extent issues one has to first reclass, merge and then reclass the merged Kahoo to NA
    if ("Maui" %in% spIslandNames) {
      rcl <- c(0, 10000, -1)
      rcl <- matrix(rcl, ncol=3, byrow=TRUE)
      a = reclassify(Kahoolawe, rcl)  # defining the Kahoolawe raster values as -1 so that they can be distinguished once merged with Maui
      Maui = try(merge(a, Maui), T)
      
      rcl <- c(-1, NA)
      rcl <- matrix(rcl, ncol=2, byrow=TRUE)
      Maui = try(reclassify(Maui, rcl), T)   
      Maui = stack(Maui)
      names(Maui)<- var_names
      spIslandNames <- spIslandNames[which(spIslandNames != "Kahoolawe")] #remove Kahoolawe from the list of islands with species
    }
    
    ###START OF CODE TO RUN PROJECTIONS ONE ISLAND AT A TIME###
    spIsland = spIslandNames[1] #for testing - run this whole loop for testing on a species.
    for (spIsland in spIslandNames){
      workspace_name_out0 = paste0(sp_nm, spIsland, "_FB_all_model_proj_", proj_nm, ".RData") #sets location to save R data
      if (file.exists(workspace_name_out0) == F | overwriteData == T){  #does not run if RData file already exists and overwrite is off  
    
        predictors <- get(spIsland)

        ###changing raster for the "wettest" and "driest" scenarios 
        alt_scen = c(2,3,5,6)
        if (baseline_or_future %in% alt_scen){
          predictors<-stack((subset(predictors, 1)),
                            (subset(predictors, 2)),
                            (subset(predictors, 3))*4,
                            (subset(predictors, 4)))
          names(predictors)<- var_name
        }
        
        if (baseline_or_future == 1){ #for "baseline" runs
          jpeg_name = paste0(sp_nm, "_", spIsland, "_env_vars_used_for_projection.jpg") #assigns location for map jpeg file
          jpeg(jpeg_name, #settings for map jpeg file
               width = 10, height = 10, units = "in", pointsize = 12, 
               quality = 90, bg = "white", res = 300) 
          plot(predictors, col = rev(terrain.colors(255)), maxpixels = 100000, 
               useRaster = useRasterDef, axes = TRUE, addfun = NULL, 
               interpolate = interpolateDef) #builds raster map and sends to jpeg file  
          dev.off() #saves plot to jpeg
        }
        cat('\n', sp_nm, "_", spIsland,'projection raster stack created...') #sign-posting
        gc() #reclaims memory that is no longer used and returns summary of memory usage
        
        if (length(spIslandNames) > 1) {
          projection_name = paste0(proj_nm, "_", spIsland)
        } else {
          projection_name = proj_nm
        }
        
        myBiomodProj_baseline <- BIOMOD_Projection(
          modeling.output = myBiomodModelOut, #results from previous model step
          new.env = predictors, #new environment to project onto (in case of baseline it is not new)
          proj.name = projection_name, #name for save folder
          selected.models = remaining_models, #whether all of a subset of the models should be used
          binary.meth = eval_stats, #vector of evaluation statistics to use for projection to presence/ absence
          compress = 'xz', #compression format for files
          build.clamping.mask = clampingMask, #whether a clamping mask should be saved or not
          keep.in.memory = memory) #whether or not the clamping mask should be saved to hard disk
        gc() #reclaims memory that is no longer used.
        cat('\n',sp_nm,'projection complete...') #sign-posting
        cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn') #sign-posting
        save("myBiomodProj_baseline", file = workspace_name_out0)   #save projection workspace  
        
        ###################################################
        ### code chunk number 17: projection_current_plot
        ###################################################
        
        #     if (plot_graphs == T){ #set in the config file
        #       jnk = length(myBiomodProj_baseline@models.projected) #gets number of projection models from output
        #       jnk = jnk/(length(models_to_run)+1) #calculates the number of projections per model type - +1 because there is a combined model type as well
        #       if (jnk<26){ 
        #         jpeg_nam = paste0(sp_nm,"_", proj_nm, "_", "runs.jpg") #assigns location for map jpeg
        #         jpeg(jpeg_name, #settings for the jpeg map
        #              width = 10, height = 10, units = "in",
        #              pointsize = 12, quality = 90, bg = "white", res = 300)
        #         par(mfrow = c(1,2)) #subsequent figures are drawn in an array of 1 by 2
        #         try(plot(myBiomodProj_baseline, str.grep = "Full"), TRUE) #attempts to plot results of full models but saves error if not
        #         dev.off() #turns device off - saves plot as jpeg
        #       }
        #     }
        
        #     if (plot_graphs == T){ ###this has to be fixed- in situations where there is only enough data for a single pseudo absence, no others are run
        #       for (model in models_to_run){ #runs for each model
        #         jpeg_name = paste0(sp_nm,"_", model, "_BIN_model", proj_nm, "runs.jpg")
        #         if (file.exists(jpeg_name) == F | overwriteData == T){ #does not run if file exists and overwrite is off
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
        #    if (plot_graphs == T){ #set in config file
        #      for (ii in 1:length(eval_stats)){ #for each evaluation statistic        
        #        for (i in 1:length(models_to_run)){ #for each model type
        #          sp_nm = str_replace_all(sp_nm,"_",".") #replaces all "_" with "." to account for biomod2 naming
        #          spEvalGrdDir <- paste0(working_dir, "/", sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, '_', sp_nm, "_", eval_stats[ii], "bin.grd") #points to location of bin.grd file
        #          spEvalGrdStack <- stack(spEvalGrdDir) #creates a raster stack from the bin.grd file
        #          spEvalGrdSub <- subset(spEvalGrdStack, (length(names(spEvalGrdStack)) - (length(models_to_run)))-1+i) #selects a layer of the rasterStack
        #          spEvalGrdNames <- names(spEvalGrdSub) #assigns names of the rasterStack Layer to a new vector
        #          rclassVect <- c(NA, 0) #creates relassification vector to change 'NA' to '0'
        #          rclassMatrix <- matrix(rclassVect, ncol=2, byrow=TRUE) #creates reclassification matrix from the reclassification vector
        #          spEvalGrdReclass <- reclassify(spEvalGrdSub, rclassMatrix) #reclassifies all 'NA' as '0' in bin.grd raster layer file
        #          names(spEvalGrdReclass) <- spEvalGrdNames #gives name from original raster layer to reclassified layer.
        #          spGrdDir <- paste0(working_dir, "/", sp_nm,"/proj_", proj_nm, "/proj_", proj_nm,"_", sp_nm, ".grd") #assigns file location for species grd file
        #          spGrdStack <- stack(spGrdDir) #loads species grd file as a raster stack
        #          spGrdSub <- subset(spGrdStack, (length(names(spGrdStack)) - (length(models_to_run)))-1+i) #selects a layer of the raster stack
        #          spGrdSub2 <- subset(spGrdSub, length(names(spGrdSub))) #selects a layer of the previous raster layer???
        #          spGrdComb <- spEvalGrdReclass * spGrdSub2 #combines the species raster layer and the species evaluation raster layer
        #          reclassVect2 <- c(0, NA) #creates reclassification vector to turn "0" into "NA"
        #          reclassMatrix2 <- matrix(reclassVect2, ncol=2, byrow=TRUE) #creates reclassification matrix from vector
        #          spGrdCombReclass <- reclassify(spGrdComb, reclassMatrix2) #reclassifies raster layer, turning "0" onto "NA"
        #          names(spGrdCombReclass) <- spEvalGrdNames #assigns original layer name to reclassified layer
        #          scaledBinEMStackNm <- paste0(eval_stats[ii], '_', models_to_run[i], "_scaled_andbinnedEM_pmw") #filename for scaled and binned EM 
        #          assign(scaledBinEMStackNm, spGrdCombReclass) #assigns combined raster layer to a filename
        #          binnedEMStackNm <- paste0(eval_stats[ii], '_', models_to_run[i], "_binnedEM_pmw") #filename for binned EM
        #          assign(binnedEMStackNm, spEvalGrdSub) #assigns evaulation raster layer to a filename
        #        }
        #        
        #        scaledBinEMStack <- stack(get(paste0(eval_stats[ii], '_', models_to_run[1], "_scaled_andbinnedEM_pmw"))) #creates a raster stack from the first raster layer created in the loop above
        #        for (i in 2:length(models_to_run)){ #for all model types
        #          scaledBinEMStack <- addLayer(scaledBinEMStack, get(scaledBinEMStackNm)) #combine raster layers from all other model types into raster stack
        #        }
        #       
        #        binnedEMStack <- stack(get(paste0(eval_stats[ii], '_', models_to_run[1], "_binnedEM_pmw"))) #creates a raster stack from the first raster layer created in the loop above
        #        for (i in 2:length(models_to_run)){ #for all model types
        #          binnedEMStack <- addLayer(binnedEMStack, get(binnedEMStackNm)) #combine raster layers from all other model types into raster stack
        #        }
        
        #        setwd(plots) #set working directory 
        #        jpeg_name = paste0(proj_nm,"_", sp_nm,"_All_", eval_stats[ii], 
        #             "Models_BinandScaled_runs_.jpg") #assigns name for jpeg
        #        jpeg(jpeg_name, width = 5*length(models_to_run), height = 5, 
        #             units = "in", pointsize = 12, quality = 90, bg = "white", 
        #             res = 300) #creates jpeg 
        
        #        par(pin = c(4,4), cex = 1, cex.main = 1, cex.axis = 0.8, 
        #            mfcol=c(1,length(models_to_run)), mgp = c(1, 0.5, 0),
        #            mar=c(2, 2, 1.5, 0), oma = c(0, 0, 0, 1), 
        #            bg = "transparent") #set graphical parameter
        
        #        gc = c('antiquewhite1', 'transparent')
        #        #grd <- terrain.colors(255)
        #        col5 <- colorRampPalette(c('blue', 'sandybrown', 'darkgreen'))
        
        #        jnk <- subset(binnedEMStack, 1) #junk variable assigned to first layer of raster stack
        #        try(plot(scaledBinEMStack, 1,  col = col5(255), 
        #                 useRaster = useRasterDef, axes = TRUE, addfun = F, 
        #                 interpolate = interpolateDef, legend = F, add = F, 
        #                 bg = "transparent"), silent = T) #plot first layer of scaled and binned raster stack with error recovery
        #        plot(jnk, col = gc, useRaster = useRasterDef, axes = F, 
        #             addfun=F, interpolate = interpolateDef, legend = F, 
        #             add = T) #plot first layer of binned raster stack 
        
        #        par(mar = c(2, 0, 1.5, 3)) #change graphical parameters *comment this out if 3 models to run
        #        jnk <- subset(binnedEMStack, 2) #junk variable assigned to second layer of binned raster stack
        #        try(plot(scaledBinEMStack, 2, col = col5(255), 
        #                 useRaster = useRasterDef, axes = TRUE,  
        #                 interpolate = interpolateDef, legend = T, yaxt = 'n', 
        #                 add = F, bg = "transparent"), 
        #            silent = T)  #plot second layer of scaled binned raster stack
        #        plot(jnk, col = gc, useRaster = useRasterDef, axes = F, 
        #             interpolate = interpolateDef,  legend = F, add = T)# "yaxt = 'n'" - plot second layer of binned raster stack
        
        #         par(mar = c(2, 0, 1.5, 3)) #change graphical parameters
        #         jnk <- subset(binnedEMStack, 3) #assign junk variable to 
        #         plot(scaledBinEMStack, 3, col = col5(255), useRaster = useRasterDef, axes = T, interpolate = interpolateDef, legend = T, yaxt = 'n', add = F, bg = "transparent") #plot third layer of binned raster stack
        #         plot(jnk, col = gc, useRaster = useRasterDef, axes = F, interpolate = interpolateDef, legend = F, add = T) #plot third layer of scaled and binned raster stack
        #         
        #        legend("bottomright",legend = c("Absent"), fill = gc[1], cex = 0.8) #creating legend for saved plot
        #        dev.off() #save plot to file
        #      }
        #    }
        #    setwd(working_dir) #return working directory to default
        #    cat('\n',sp_nm,'done with individual model plots...') #sign-posting
        
      } else {
        load(workspace_name_out0) #loads existing workspace
        cat('\n',sp_nm,'projection of individual models loaded from past run...') #sign-posting
      }
      
      if (apply_biomod2_fixes){ #parameter set in config file
        myBiomodProjection <- LoadProjectionManually(myBiomodProj_baseline) #function set in another module?
        
        } else {
          myBiomodProjection <- myBiomodProj_baseline
        }
      
      cat('\n',sp_nm,'projection graphs done...') #sign-posting
      
      workspace_name_out1 = paste0(sp_nm, spIsland, "_FB_all_model_proj_EF", proj_nm, ".RData") #sets location to save R data
      if (file.exists(workspace_name_out1) == F | overwriteData == T){  #does not run if RData file already exists and overwrite is off  
        ###################################################
        ### code chunk number 19: EnsembleForecasting
        ###################################################
        myBiomodEF <- BIOMOD_EnsembleForecasting( #Ensemble projections of species
          projection.output = myBiomodProjection,
          total.consensus = T, #setting em.by to all to combine all models
          EM.output = myBiomodEM, #from module 2 BIOMOD_EnsembleModeling output
          binary.meth = eval_stats, #names of evaluation metrics - defined in config module
          keep.in.memory = memory)
        cat('\n',sp_nm,'ensemble projection done...') #sign-posting
        #cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn') #returns memory used
        
        ###################################################
        ### code chunk number 20: EnsembleForecasting_loading_res
        ###################################################
        
        #plotting the ensemble projections per species per projection
        #    if (plot_graphs == T){ #set in config file
        #      for (i in 1:length(eval_stats)){ #for each evaluation statistic
        #        totalConsDir1 <- paste0(working_dir, "/", sp_nm, "/proj_", proj_nm, 
        #                                  "/proj_", proj_nm, "_", sp_nm, 
        #                                  "_TotalConsensus_EMby", eval_stats[i], 
        #                                  ".grd") #sets location of total consensus ensemble model .grd file
        #        totalConsStack1 = stack(totalConsDir1) #creates a raster stack from the .grd file
        #        totalConsSub1 <- subset(totalConsStack1, length(names(totalConsStack1))) #creates a raster layer from the .pmw in the raster stack
        #        
        #        #WHAT IS THE DIFFERENCE BETWEEN THIS STACK AND PREVIOUS?
        #        totalConsDir2 <- paste0(working_dir, "/", sp_nm, "/proj_", proj_nm, 
        #                                "/proj_", proj_nm, "_", sp_nm, 
        #                                "_TotalConsensus_EMby", eval_stats[i], "_", 
        #                                eval_stats[i], "bin.grd") #sets location of total consensus .grd file
        #        totalConsStack2 = stack(totalConsDir2) #creates raster stack from .grd file
        #        totalConsSub2 <- subset(totalConsStack2, length(names(totalConsStack2))) #creates raster layer from .pmw file
        #        
        #        totalConsComb <- totalConsSub1 * totalConsSub2 #combines the two raster layers from the .pmw files
        #        totalConsCombReclass <- reclassify(totalConsComb, reclassMatrix2) #reclassifies any values of 0 into NA to get rid of island outline
        
        #        names(totalConsCombReclass) <- names(totalConsSub1) #assigns layer name from the original .pmw file to new reclassified layer
        #        assign(paste0("TotalConsensus_EMScaledandBinnedby_", eval_stats[i]), totalConsComb) #assigns the combined raster layers to a character string
        #        assign(paste0("TotalConsensus_EMBinnedby_", eval_stats[i]), totalConsSub2) #assigns the raster layer from the 2nd .pmw file to a character string
        #      }
        
        #      emsScaledBinStack <- stack(get(paste0("TotalConsensus_EMScaledandBinnedby_", eval_stats[1])))
        #      emsBinStack <- stack(get(paste0("TotalConsensus_EMBinnedby_", eval_stats[1])))
        #      if(length(eval_stats)>1){
        #        for (i in 2:length(eval_stats)){
        #          emsScaledBinStack <- addLayer(emsScaledBinStack, get(paste("TotalConsensus_EMScaledandBinnedby_", eval_stats[i], sep="")))
        #          emsBinStack <- addLayer(emsBinStack, get(paste("TotalConsensus_EMBinnedby_", eval_stats[i], sep="")))
        #      }}
        
        
        #      setwd(plots)
        
        #      jpeg_name=paste(proj_nm,"_", sp_nm,"_TOTALCONSENSUS_Binandscaled_runs_.jpg", sep = "")
        #      jpeg(jpeg_name, width = 5*length(eval_stats), height = 5, units = "in",
        #           pointsize = 12, quality = 90, bg = "white", res = 300)  
        #      par(pin = c(4,4), cex = 1, cex.main = 1, cex.axis = 0.8, mfcol=c(1,length(eval_stats)), mgp = c(1, 0.5, 0),
        #          mar=c(2, 2, 1.5, 1), oma = c(0, 0, 0, 1), bg = "transparent")
        
        #      gc = c('antiquewhite1', 'transparent')
        #      col5 <- colorRampPalette(c('blue', 'sandybrown', 'darkgreen'))
        #      jnk <- subset(emsBinStack, 1)
        #      try(plot(emsScaledBinStack[[1]],  col = col5(255), useRaster=useRasterDef, axes = TRUE, addfun=F, interpolate = interpolateDef, legend = F, add = F, bg = "transparent"),silent=T)
        #      plot(jnk, col = gc, useRaster=useRasterDef, axes = F, addfun=F, interpolate = interpolateDef, legend = F, add = T)
        
        #      if (length(eval_stats)>1){ 
        #        par(mar=c(2, 0, 1.5, 0))
        #        jnk <- subset(emsBinStack, 2)
        #        try(plot(emsScaledBinStack[[2]], col = col5(255), useRaster=useRasterDef, axes = TRUE,interpolate = interpolateDef, legend = F, yaxt = 'n', add = F, bg = "transparent"),silent=T)  # addfun=F, 
        #        plot(jnk, col = gc, useRaster=useRasterDef, axes = F, addfun=F, interpolate = interpolateDef, yaxt = 'n', legend = F, add = T)
        #      }
        
        #      if (length(eval_stats)>2){ 
        #        par(mar=c(2, 0, 1.5, 3.5))
        #        jnk <- subset(emsBinStack, 3)
        #        try(plot(emsScaledBinStack[[3]], col = col5(255), useRaster=useRasterDef, axes = T, interpolate = interpolateDef, legend = T, yaxt = 'n', add = F, bg = "transparent"),silent=T)
        #        plot(jnk, col = gc, useRaster=useRasterDef, axes = F, interpolate = interpolateDef, legend = F, add = T)
        #      }
        
        #     legend("bottomright",legend = c("Absent"), fill = gc[1], cex = 0.8)
        #      dev.off()
        #    }  
        
        #    setwd(working_dir)
        #    if (plot_graphs==1){
        #      for (eval_stat in eval_stats){
        #        try(load(paste0(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat)), TRUE)    
        #        try(load(paste0(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.",eval_stat)), TRUE)
        #        jpeg_name=paste0(sp_nm,"_", eval_stat,"_ensemble_", proj_nm, "runs.jpg")
        #        jpeg(jpeg_name,
        #             width = 10, height = 8, units = "in",
        #             pointsize = 12, quality = 90, bg = "white", res = 300)
        #        par(mfrow=c(1,2))
        #        try(plot(get(paste0(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat))), TRUE)
        #        try(plot(get(paste0(sp_nm,"_AllData_AllRun_EM.",eval_stat))), TRUE)
        #        sp_bin_file=paste(proj_nm, "_", sp_nm, "_TotalConsensus_EMby", eval_stat,".grd", sep = "")
        #        sp_bin_file=paste(sp_nm,"/proj_", proj_nm, "/proj_", sp_bin_file , sep = "") #current_BI_Akepa_bin_ROC_RasterStack          
        #        try(plot(raster(sp_bin_file)), TRUE)
        #        dev.off()
        #        eval_stat0=eval_stat
        #        for (eval_stat in eval_stats){
        #          try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat0,".bin.",eval_stat, sep = "")), TRUE)    
        #          try(load(paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.",eval_stat0,".bin.",eval_stat, sep = "")), TRUE)
        #          jpeg_name=paste(sp_nm,"_", eval_stat0,"_ensemble_", proj_nm, "_bin_",eval_stat,"runs.jpg", sep = "")
        #          jpeg(jpeg_name,
        #               width = 10, height = 8, units = "in",
        #               pointsize = 12, quality = 90, bg = "white", res = 300)
        #par(mfrow=c(1,2))
        #          try(plot(get(paste(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat0,".bin.",eval_stat, sep = ""))), TRUE)
        #          try(plot(get(paste(sp_nm,"_AllData_AllRun_EM.",eval_stat0,".bin.",eval_stat, sep = ""))), TRUE)
        #          sp_bin_file=paste(proj_nm, "_", sp_nm, "_TotalConsensus_EMby", eval_stat0,"_",eval_stat, "bin.grd", sep = "")
        #          sp_bin_file=paste(sp_nm,"/proj_", proj_nm, "/proj_", sp_bin_file , sep = "") #current_BI_Akepa_bin_ROC_RasterStack          
        #          try(plot(raster(sp_bin_file)), TRUE)
        #          dev.off()
        #        }
        #      }
        #    }
        #    cat('\n',sp_nm,'ensemble projection figures done...')
        save("myBiomodProjection", "myBiomodEF", file = workspace_name_out1)   #save workspace
        
        removeTmpFiles(h=1)
        cat('\n',sp_nm, "_", spIsland,'done...')
      } else {
        cat('\n',sp_nm, "_", spIsland,'previously done...')
      }
        
    }
    ###END OF CODE TO PROJECT ONE ISLAND AT A TIME###
    

    ####START OF CODE TO MERGE RASTERS FROM EACH ISLAND INTO RASTER ONE FOR ALL ISLANDS###
    combinedDir = paste0(working_dir, sp_dir, "/proj_", proj_nm) #names the directory where the combined rasters will be placed
    dir.create(combinedDir, showWarnings = FALSE) #creates the named directory above
    island1Dir = paste0(working_dir, sp_dir, "/proj_", proj_nm, "_", spIslandNames[1], "/") #points to the name of the 1st island directory
    fileList <- list.files(island1Dir, pattern = "*.gri") #lists all files within the 1st island directory
    for (file in fileList) { #for each file in the file list
      rasterStack <- stack(paste0(island1Dir, file)) #create a rasterstack from the *.gri file in the list
      rasterStackNames = names(rasterStack) #assigns names of the original rasterStack to a vector to rename later
      if (length(spIslandNames) > 1) { #if there is more than one island covered by the extent of the species
        for (islandNum in 2:length(spIslandNames)) { #for each of the islands covered by the data except the 1st island
          newIslandName = spIslandNames[islandNum] #find the name of the island
          tempIslandDir = str_replace_all(island1Dir, spIslandNames[1], newIslandName) #finds the directory for the island by replacing 1st island name with the new island name 
          tempFile = str_replace_all(file, spIslandNames[1], newIslandName) #finds the filename for the island by replacing the 1st island name with the new island name
          tempStack <- stack(paste0(tempIslandDir, tempFile)) #creates a rasterstack from the *.gri file for the new island
          rasterStack = stack(merge(rasterStack, tempStack)) #merges the new island data (stack) with the previous island data and creates a RasteStack with the result (need stack() because otherwise a RasterBrick is created)
        }
        names(rasterStack) <- rasterStackNames #assigns the names of the original rasterStack to the combined one
      }
      combinedFileName = str_replace_all(file, paste0("_", spIslandNames[1]), "") #creates a filename for the merged rasterstack by removing the island name
      combinedFileLoc = paste0(combinedDir, "/", combinedFileName) #creates a file location for the file with data for all islands
      writeRaster(rasterStack, combinedFileLoc, overwrite = TRUE) #saves the RasterStack to a file and location indicated above  
    }   
    ####END OF CODE TO MERGE RASTERS FROM EACH ISLAND INTO ONE FOR ALL ISLANDS###
        
  } else {
    cat('\n',sp_nm,'previously calculated...')
  }
}  


    
            
    



