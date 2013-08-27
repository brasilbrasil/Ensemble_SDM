###USER CONFIGURATION
#see 0_sdm_config.r file

###START UNDERHOOD
setwd(working_dir) #sets the working directory

#Load required libraries
library(biomod2)
library(raster)
library(randomForest)
library(dismo)
library(mda)
library(stringr)
library(tools)

###not in FWS code (copy necessary files)
#this loop copies the necessary data to run the models into the working directory
dirs = list.dirs(necessary_run_data, full.names = FALSE, recursive = TRUE)
for (dir in dirs){
  layers <- list.files(dir, pattern=NULL, full.names=FALSE, include.dirs = FALSE)
  for (layer in layers){
    layer_full_nm = paste(dir,layer, sep="/")
    if (file.info(layer_full_nm)$isdir==FALSE){
      out_dir_nm = str_replace(dir, necessary_run_data, working_dir)
      dir.create(out_dir_nm, showWarnings = FALSE, recursive = TRUE, mode = "0777")
      out_lyr_nm = str_replace(layer_full_nm, necessary_run_data, working_dir)
      if (file.exists(out_lyr_nm)==FALSE){
        cat('\n','found ', layer, 'in ', dir)
        file.copy(layer_full_nm, out_lyr_nm, overwrite = TRUE, recursive = TRUE,
                  copy.mode = TRUE)
        cat('\n','saved as ', out_lyr_nm)
      }
    }
  }
}
cat('\n','Copying of necessary files is complete')

#assigns where to find species information (e.g. raster directory, raster set, endemic...) in a csv file )
spp_info = read.csv(paste0(csv_dir,'/FB_spp_data.csv'))

#creates vector with bioclimatic variable names without the file extension (.grd)
var_name <- unlist(file_path_sans_ext(env_var_files))

sp_nm = spp_nm[1]
n_abs_removed = c()
for (sp_nm in spp_nm){
  sp_nm = as.character(sp_nm) #converts sp_nm into character variable (in case the species are numbered)
  sp_dir = str_replace_all(sp_nm,"_", ".") #replaces "_" with "." in sp_nm
  dir.create(sp_dir, showWarnings = FALSE) #creates new directory with species name
  
  cat('\n',sp_nm,'model fitting...') #sign-posting
  FileName00 <- paste0(sp_nm, "_VariImp.csv") ###not in FWS code - allows for overwrite capacity
  if (file.exists(FileName00) == FALSE | overwrite == 1){ #check if analysis for species already done or overwrite requested in config    
    # Start the clock!
    ptm0 <- proc.time()
    workspace_name = paste0(sp_nm,"_FB_modelfitting.RData") #set name of file to save all workspace data after model run
    
    #######Loading datasets#######
    
    ##raster_based_env_grid:
    cat('\n','loading rasters...') #sign-posting
    
    sp_index = which(spp_info[,"Species"] == sp_nm) #finds which line of species .csv file has information needed
    raster_res = paste0("/", spp_info[sp_index,"rasterdir"]) #finds which raster directory should be used for the species based on the .csv file
    crop_raster = raster(paste0(crop_raster_dir,raster_res,".grd")) #assigns cropped raster associated with sp_nm to "crop_raster" variable
    predictors = raster(paste0(fitting_clim_data_dir, "/", env_var_files[1])) #assigns bioclimate raster to "predictors" variable
    predictors = crop(predictors,  crop_raster) #crops predictor grid using crop_raster
    jnk0 = length(env_var_files) #creates variable with no. of bioclimate variables to use
    for (jj in 2:jnk0){ #this loop adds the rest of the rest of the bioclimate variables to the "predictors" raster
      temp = raster(paste0(fitting_clim_data_dir, "/", env_var_files[jj]))
      temp = crop(temp,  crop_raster)
      predictors = addLayer(predictors, temp)
    }
    names(predictors) <- var_name #assigns names to bioclimate raster stack
    rm("jnk0","jj","crop_raster" ,"temp") #removes temporary variables
    
    jpeg_name = paste0(sp_nm,"_env_vars_used.jpg") #names jpeg file to be created
    jpeg(jpeg_name, #creates blank jpeg file in working directory
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(predictors, col=rev(terrain.colors(255)), maxpixels = 100000, useRaster = FALSE, axes = TRUE, addfun = NULL) #
      #ERROR - 24 warnings to do with "interpolate"  - check with warnings()  
    dev.off()
    
    ####species point data
    cat('\n','loading species data...') #sign-posting
    mySpeciesOcc = read.csv(paste0(csv_dir,"/",sp_nm,'_pres_abs.csv')) #FB_data_points4_PAandA
    
    #presence (and absence) data handling)
    mySpeciesOcc = cbind(mySpeciesOcc[,2:3],pa=mySpeciesOcc[,1]) #extracts just the x, y, and PA data from the species csv
    head(mySpeciesOcc) #check header of PA data
    if (!include_Abs){ ##NEW
      mySpeciesOcc=mySpeciesOcc[mySpeciesOcc$pa==1,]
    }
    
    ##pseudo-absence handling
    cat('\n','defining candidate PA points...')
    PA_XY = mySpeciesOcc[,1:2] #extracts the lat and long for all PA points
    mySREresp <- reclassify(subset(predictors,1,drop=TRUE), c(-Inf,Inf,0)) #builds a raster layer based on environmental rasters for response variable
    mySREresp[cellFromXY(mySREresp,PA_XY)] <- 1 #assigns all shared cells in "bioclim" and "PA_XY" to "1"
    Act_abs = dim(mySpeciesOcc[mySpeciesOcc$pa==0,])[1] #calculates number of real absences

    #this loop makes different runs depending on whether pseudo absences should be assigned outside climate envelope or randomly
    if (PAs_outside_CE){
      sre_tail = 0.025
      sp_CE = sre(Response = mySREresp, Explanatory = predictors,NewData = predictors, Quant = sre_tail) #Calculates surface range envelope for distribution removing a percentage of extremes (based on sre_tail
      #calculate density of points within sre
      n_PandA = sum(as.matrix(mySREresp), na.rm=T)*(1-sre_tail*2) #this corrects value for density calc #counts number of cells where species occurence data available
      CE_cells = sum(as.matrix(sp_CE), na.rm=T) #counts number of cells where presence or absence predicted in climate envelope
      CE_point_density = round(n_PandA/CE_cells, digits = 4) #density of points with real data within the climate envelope
      
      #creates raster of all cells outside CE
      neg_sp_CE = sp_CE == 0 
      
      #calculate desired number of pseudo absence points outside CE based on PandA density within CE
      neg_CE_cells = sum(as.matrix(neg_sp_CE), na.rm=T) #calculates no. cells outside CE
      n_PseudoAbs_pts = round(neg_CE_cells * dens_PAs_outside_CE * CE_point_density) + Act_abs #accounts for the actual absences for calculating pseudo absences      
      PseudoAbs_cand_pts = rasterToPoints(neg_sp_CE, fun = function(x){x == 1}) #Creates matrix of candidate points (x,y,layer)
      
      plot(mySREresp) #plots all cells with data
      plot(sp_CE) #plots climate envelope
      plot(neg_sp_CE) #plots areas outside climate envelope

      # next section assigns pseudo absences anywhere (not limited to CE)       
    }else{
      neg_mySREresp = mySREresp == 0 #creates raster of areas outside those with known data (presence and absence)
      plot(neg_mySREresp) #plots raster outside known data (presence and absence)
      PseudoAbs_cand_pts = rasterToPoints(neg_mySREresp, fun=function(x){x==1}) # Creates matrix of candidate pseudo absence points (x,y,layer)
      if (candidatePA.per.PA==0){
        n_PseudoAbs_pts = PA.nb.absences  #assigns number of pseudo absence points as indicated in config code - SHOULD WE ACCOUNT FOR ACTUAL ABSENCES?
      }else{
        n_PseudoAbs_pts = round(dim(PseudoAbs_cand_pts)[1]/candidatePA.per.PA)        
      }
    }  
    PseudoAbs_cand_pts = as.data.frame(PseudoAbs_cand_pts[,1:2]) #extracts only the geographic information for the candidate pseudo absence points
    head(PseudoAbs_cand_pts) #returns the first lines for the pseudo absence candidate pts 
    dim(PseudoAbs_cand_pts) #returns dimensions of the pseudo absence points (#pts and 2 rows)
    PseudoAbs_cand_pts_noNA=PseudoAbs_cand_pts[complete.cases(PseudoAbs_cand_pts),] #creates new data frame removing any pts with missing geographic information
    PseudoAbs_cand_pts_noNA=cbind(PseudoAbs_cand_pts_noNA,pa=rep('NA', dim(PseudoAbs_cand_pts_noNA)[1],1)) #adds "pa" data column to dataframe and assigns "NA" to all rows for that column
    #IRC: Is this next line necessary - seems like it already has a header with lower case x and y
    names(PseudoAbs_cand_pts_noNA)=c('X', 'Y', 'pa') #reassigns names to capitalizes X and Y
    head(PseudoAbs_cand_pts_noNA) #returns the first lines of the data frame
    
    #merge data with pseudoabsence
    mySpeciesOcc_w_Pseudo <- data.frame(rbind(mySpeciesOcc, PseudoAbs_cand_pts_noNA)) #creates new data frame with real data and pseudo absence candidate points combined
    
    #### EXTRACTION OF ENV DATA FOR POINT DATA
    cat('\n','extracting env vars to points...')
    relBioclimData <- extract(predictors, mySpeciesOcc_w_Pseudo[,1:2], cellnumbers = T) #creates new matrix with relavent bioclim variables and the cell numbers for the real data and candidate pseudo absence points
    XY_PresAbsPA_Bioclim = data.frame(cbind(mySpeciesOcc_w_Pseudo,relBioclimData)) #creates data frame with bioclimate data for each of the points with real data and pseudo absence candidates
    head(XY_PresAbsPA_Bioclim) #returns first lines of the dataset
    XY_PresAbsPA_Bioclim_noNA = XY_PresAbsPA_Bioclim[complete.cases(XY_PresAbsPA_Bioclim[4:dim(XY_PresAbsPA_Bioclim)[2]]),] #creates new data frame with any rows with missing data removed
    head(XY_PresAbsPA_Bioclim_noNA) #returns first lines of new data frame
    
    ### Select, Count and Remove presence duplicate points in cells 
    jnk = c(1,0, NA) #temporary vector with unique values for "pa" column
    XY_PresAbsPA_Bioclim_sort <- XY_PresAbsPA_Bioclim_noNA[order(match(XY_PresAbsPA_Bioclim_noNA$pa,jnk)),] #sorting so if duplicates, PA removed before abs, abs removed before pres
    rm(jnk)    
    
    dupEntries <- duplicated(XY_PresAbsPA_Bioclim_sort$cells) # Identifies duplicates in cell column 
    n_dups = length(dupEntries[dupEntries == TRUE]) #calculates number of duplicate entries
    cat('\n','Out of ', length(dupEntries), "points, ",n_dups, "were removed because they were within the same raster cell for", sp_nm) #sign-posting
    PresAbsPA_noDup <- XY_PresAbsPA_Bioclim_sort[!dupEntries, -4] #creates new dataframe with duplicates and cell column removed
    head(PresAbsPA_noDup)
    
    #write.table(PresAbsPA_noDup, file = paste0(sp_nm,"_loc_data_table.csv"), sep = ",", col.names = NA) #this generates a table way to large
    sp_loc_summary = table(PresAbsPA_noDup$pa, useNA = "ifany")
    sp_loc_summary = as.data.frame(sp_loc_summary)
    levels(sp_loc_summary$Var1) = c(0, 1, NA, "n to select")
    jnk = c("n to select", n_PA_points)#3 is code for n of pseudo absences to select per run
    sp_loc_summary = rbind(sp_loc_summary, jnk)
    names(sp_loc_summary) = c("Data_type", "Freq")
    write.table(sp_loc_summary, file = paste0(sp_nm,"_loc_data_summary.csv"), sep = ",", col.names = NA)

    
    ###PLOTS MAPS OF DATA (PRESENCE, ABSENCE, AND PSEUDO ABSENCE)
    jpeg_name2=paste0(sp_nm,"_loc_data_used.jpg") #assigning filename to jpeg
    jpeg(jpeg_name2,
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300) #creating blank jpeg file
    plot(seq((min(PresAbsPA_noDup[,1])-0.1),(max(PresAbsPA_noDup[,1])+0.1),by=((max(PresAbsPA_noDup[,1])+0.1)-(min(PresAbsPA_noDup[,1])-0.1))/5), 
         seq((min(PresAbsPA_noDup[,2])-0.1),(max(PresAbsPA_noDup[,2])+0.1),by=((max(PresAbsPA_noDup[,2])+0.1)-(min(PresAbsPA_noDup[,2])-0.1))/5), 
         type = "n", xlab="Lon", ylab="Lat")# setting up coord. system
    points(x=PresAbsPA_noDup[is.na(PresAbsPA_noDup[,3]),1], y=PresAbsPA_noDup[is.na(PresAbsPA_noDup[,3]),2], type = "p", col = "grey", pch=20,cex = 0.7)
    points(x=PresAbsPA_noDup[PresAbsPA_noDup[,3]==0,1], y=PresAbsPA_noDup[PresAbsPA_noDup[,3]==0,2], type = "p", col = "red", pch=20,cex = 0.7)
    points(x=PresAbsPA_noDup[PresAbsPA_noDup[,3]==1,1], y=PresAbsPA_noDup[PresAbsPA_noDup[,3]==1,2], type = "p", col = "blue", pch=20,cex = 0.7)
    dev.off()
    
    
    ###defining the variables used by biomod2
    cat('\n','biomod model config...') #sign-posting
    myRespXY = PresAbsPA_noDup[,1:2] #new data frame with the x and y coordinates only
    myResp <- data.frame(Sp_Bio = PresAbsPA_noDup[,3]) #new data frame with only one column for whether it is a presence (1), absence (0), or PA (NA)
    myResp[myResp =='NA'] = NA #re-assigns all "NA" (text) to real <NA> 
    myExpl = PresAbsPA_noDup[,4:dim(PresAbsPA_noDup)[2]] #new data frame with only the explanatory variables
    
    myBiomodData <- BIOMOD_FormatingData(
      resp.var = myResp,
      expl.var = myExpl,
      resp.xy = myRespXY,
      resp.name = sp_nm,
      PA.nb.rep = PA.nb.rep,
      PA.nb.absences = n_PseudoAbs_pts,
      PA.strategy = PA.strategy,
      PA.dist.min = PA.dist.min) #checking formatting of biomod 2 variables
    
    #This plotting methods takes way too long!!!  (but it is useful since it plots PAs selected)
    if (plot_graphs == 1 & PA.nb.rep < 9){  
    jpeg_name3=paste0(sp_nm,"_loc_data_used2.jpg") #assigning location for jpeg file
    jpeg(jpeg_name3,
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300) #creating blank jpeg file
    plot(myBiomodData) #plots the psuedo absences, real data, and candidate points used
    dev.off()
    }
    
    #memory.limit(size=4095)
    myBiomodOption <- BIOMOD_ModelingOptions(
      GBM = list( distribution = 'bernoulli', interaction.depth = 7,  shrinkage = 0.001, bag.fraction = 0.5, train.fraction = 1, n.trees = 100,
                  cv.folds = 10),
      MARS = list( degree = 2,
                   penalty = 2,
                   thresh = 0.001,
                   prune = TRUE),
      RF = list(do.classif = TRUE, ntree = 100, mtry = 'default', max.nodes=10, corr.bias = T), 
      MAXENT = list(maximumiterations = 100, visible = F, linear = TRUE, quadratic = TRUE,
                    product = TRUE, threshold = TRUE, hinge = TRUE, lq2lqptthreshold = 80, l2lqthreshold = 10,
                    hingethreshold = 15, beta_threshold = -1, beta_categorical = -1, beta_lqp = -1, 
                    beta_hinge = -1,defaultprevalence = 0.5)
    )
    
    cat('\n','fitting...') #sign-posting
    
    #CREATING ENSEMBLE MODEL
    myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
                                        models = models_to_run, #from config file
                                        models.options = myBiomodOption, #from above
                                        NbRunEval = NbRunEval, #from config file
                                        DataSplit = 80,
                                        Yweights = NULL, 
                                        VarImport = 10,
                                        do.full.models = do.full.models,
                                        models.eval.meth = eval_stats, #c('TSS','ROC', 'KAPPA'),  #options set in config file
                                        SaveObj = TRUE,
                                        rescal.all.models = TRUE)
    
    ## Output the biomod models
    myBiomodModelOut #returns summary of model runs
    
    # output model evaluation metrics
    myBiomodModelEval <- getModelsEvaluations(myBiomodModelOut) #creates an array with model evaluation results for all models   
    dimnames(myBiomodModelEval) #returns the header names for the model evaluation array 
    
    # Outputting the validation metrics for all tests
    if ("TSS" %in% eval_stats){
      myBiomodModelEval["TSS","Testing.data",,,] #returns variable importances based on TSS for each model
      Spp_TSS<- data.frame(myBiomodModelEval["TSS","Testing.data",,,]) #creates data frame with TSS variable importances
      FileName<-paste0(sp_nm, "_TSS.csv") #assigns file path for results
      write.table(Spp_TSS, file = FileName, sep=",", col.names=NA) #creates csv file with TSS variable importances
    }
    
    if ("ROC" %in% eval_stats){
      myBiomodModelEval["ROC","Testing.data",,,]
      Spp_ROC<- data.frame(myBiomodModelEval["ROC","Testing.data",,,])
      FileName<-paste0(sp_nm, "_ROC.csv")
      write.table(Spp_ROC, file = FileName, sep=",", col.names=NA)
    }
    if ("KAPAA" %in% eval_stats){
      myBiomodModelEval["KAPPA","Testing.data",,,]
      Spp_KAP<- data.frame(myBiomodModelEval["KAPPA","Testing.data",,,])
      FileName<-paste0(sp_nm, "_KAP.csv")
      write.table(Spp_KAP, file = FileName, sep=",", col.names=NA)
    }
    ## getting the variable importance ##
    getModelsVarImport(myBiomodModelOut) #returns an array with model variable importances (i.e bio1, bio7, etc)
    Spp_VariImp <- data.frame(getModelsVarImport(myBiomodModelOut)) #creates data frame with model variable importances
    write.table(Spp_VariImp, file = FileName00, sep=",", col.names=NA) #creates csv of variable importances with name assigned above.
    
    save("myBiomodModelOut", file=workspace_name)   #save workspace      
     
    ptm1 = proc.time() - ptm0 #calculates time it took to run all code
    jnk = as.numeric(ptm1[3]) #assigns temporary variable to the numeric value of the time elapsed
    jnk = jnk/3600 #converts elapsed time into hours
    cat('\n','It took ', jnk, "hours to model", sp_nm) #sign-posting
  }else{
    cat('\n','fitting for ',sp_nm,'already done...') #sign-posting in case file for variable importance has already been created (indicating this species has already been run)  
  }    
}