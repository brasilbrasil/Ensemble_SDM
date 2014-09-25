###USER CONFIGURATION
#see 0_sdm_config file.r

setwd(working_dir)
require(snowfall)

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
if (baseline_or_future == 7){
  clim_data_dir = clim_data_2100rev
  proj_nm = 'future_rev'}


spp_info = read.csv(paste(csv_dir,'FB_spp_data.csv', sep = "/")) #creates data frame from species info csv file

sp_nm = spp_nm[1] #resets so the first species to run is the first one listed in config file or csv
sp_parallel_run=function(sp_nm){
  #loading package libraries
  library(biomod2)
  library(stringr)
  library(colorRamps)
  library(rasterVis)
  library(tools)
  library(ncdf)
  library(gbm) #needed for projection of "gbm" data (i.e. "BIOMOD_Projection(...)" function)
  
  #sets options for biomod2 fixes in code (if assigned TRUE in config file)
  if (apply_biomod2_fixes){
    rasterOptions(tmpdir = dir_for_temp_files, timer = T, progress = "text", todisk  = T) #set options for raster package
    source(paste(codeDir,"3b_modifications_of_projection_code.r", sep = "/")) #all of fixes to biomod2 code created by AV
  }

  sp_nm = as.character(sp_nm) #defines the species name as a character string - not needed if it is already a text name 
  sp_dir = str_replace_all(sp_nm,"_", ".") #replaces "_" with "." in sp_nm
  sink(file(paste0(working_dir,sp_dir,"/",sp_dir,Sys.Date(),"_proj_log.txt"), open="wt"))#######NEW
  
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
        
        save("myBiomodProjection", "myBiomodEF", file = workspace_name_out1)   #save workspace
        
        removeTmpFiles(h=1)
        cat('\n',sp_nm, "_", spIsland,'done...')
      } else {
        cat('\n',sp_nm, "_", spIsland,'previously done...')
      }   
    }
    ###END OF CODE TO PROJECT ONE ISLAND AT A TIME###
    

    ####START OF CODE TO MERGE RASTERS FROM EACH ISLAND INTO RASTER ONE FOR ALL ISLANDS###
    if (length(spIslandNames) > 1) { #only need to recombine if the species is found on more than one island - otherwise the combined raster is already created above
      combinedDir = paste0(working_dir, sp_dir, "/proj_", proj_nm) #names the directory where the combined rasters will be placed
      dir.create(combinedDir, showWarnings = FALSE) #creates the named directory above
      island1Dir = paste0(working_dir, sp_dir, "/proj_", proj_nm, "_", spIslandNames[1], "/") #points to the name of the 1st island directory
      fileList <- list.files(island1Dir, pattern = "*.gri") #lists all files within the 1st island directory
      for (file in fileList) { #for each file in the file list
        rasterStack <- stack(paste0(island1Dir, file)) #create a rasterstack from the *.gri file in the list
        rasterStackNames = names(rasterStack) #assigns names of the original rasterStack to a vector to rename later
        for (islandNum in 2:length(spIslandNames)) { #for each of the islands covered by the data except the 1st island
          newIslandName = spIslandNames[islandNum] #find the name of the island
          tempIslandDir = str_replace_all(island1Dir, spIslandNames[1], newIslandName) #finds the directory for the island by replacing 1st island name with the new island name 
          tempFile = str_replace_all(file, spIslandNames[1], newIslandName) #finds the filename for the island by replacing the 1st island name with the new island name
          tempStack <- stack(paste0(tempIslandDir, tempFile)) #creates a rasterstack from the *.gri file for the new island
          rasterStack = stack(merge(rasterStack, tempStack)) #merges the new island data (stack) with the previous island data and creates a RasteStack with the result (need stack() because otherwise a RasterBrick is created)
        }
        names(rasterStack) <- rasterStackNames #assigns the names of the original rasterStack to the combined one
        combinedFileName = str_replace_all(file, paste0("_", spIslandNames[1]), "") #creates a filename for the merged rasterstack by removing the island name
        combinedFileLoc = paste0(combinedDir, "/", combinedFileName) #creates a file location for the file with data for all islands
        writeRaster(rasterStack, combinedFileLoc, overwrite = TRUE) #saves the RasterStack to a file and location indicated above 
      }
    }
     
    ####END OF CODE TO MERGE RASTERS FROM EACH ISLAND INTO ONE FOR ALL ISLANDS###
        
  } else {
    cat('\n',sp_nm,'previously calculated...')
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



    
            
    



