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
        file.copy(layer_full_nm, out_lyr_nm, overwrite = TRUE, recursive = FALSE,
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

for (sp_nm in spp_nm){
  sp_nm = as.character(sp_nm) #converts sp_nm into character variable (in case the species are numbered)
  sp_dir = str_replace_all(sp_nm,"_", ".") #replaces "_" with "." in sp_nm
  dir.create(sp_dir, showWarnings = FALSE) #creates new directory with species name
  
  cat('\n',sp_nm,'model fitting...') #sign-posting
  FileName00 <- paste0(sp_nm, "_VariImp.csv") ###not in FWS code - allows for overwrite capacity
  if (file.exists(FileName00)==FALSE | overwrite==1){ #check if analysis for species already done or overwrite requested in config    
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
    rm("crop_raster" ,"temp") #removes temporary variables
    
    jpeg_name = paste0(sp_nm,"_env_vars_used.jpg") #names jpeg file to be created
    jpeg(jpeg_name, #creates blank jpeg file in working directory
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300)
    dev.off() #shuts down current device  - may need to use this before next line or it may not plot
    plot(predictors, col=rev(terrain.colors(255)), maxpixels = 100000, useRaster = FALSE, axes = TRUE, addfun = NULL) #NOT WORKING
      #ERROR - 24 warnings to do with "interpolate" check with warnings()  
    dev.off()
    
    ####species point data
    cat('\n','loading species data...') #sign-posting
    mySpeciesOcc = read.csv(paste0(csv_dir,"/",sp_nm,'_pres_abs.csv')) #FB_data_points4_PAandA
    
    #presence (and absence) data handling)
    mySpeciesOcc = cbind(mySpeciesOcc[,2:3],pa=mySpeciesOcc[,1]) #extracts just the x, y, and PA data from the species csv
    head(mySpeciesOcc) #check header of PA data
    
    ##pseudo-absence handling
    cat('\n','defining candidate PA points...')
    PA_XY = mySpeciesOcc[,1:2] #extracts the lat and long for all PA points
    mySREresp <- reclassify(subset(predictors,1,drop=TRUE), c(-Inf,Inf,0)) #builds a raster layer based on environmental rasters for response variable
    mySREresp[cellFromXY(mySREresp,PA_XY)] <- 1 #assigns all shared cells in "bioclim" and "PA_XY" to "1"
    Act_abs = dim(mySpeciesOcc[mySpeciesOcc$pa==0,])[1] #calculates number of real absences

    #this loop makes different runs depending on whether Pseudo Absences should be assigned outside climate envelope or randomly
    if (PseudoAbs_outside_CE){      
      sp_CE = sre(Response = mySREresp,Explanatory = predictors,NewData = predictors,Quant = 0.025) #Calculates surface range envelope for distribution removing 2.5% of extremes 
      #calculate density of points within sre
      n_PandA = sum(as.matrix(mySREresp), na.rm=T) #counts number of cells where species occurence data available
      CE_cells = sum(as.matrix(sp_CE), na.rm=T) #counts number of cells where presence or absence predicted in climate envelope
      CE_point_density = round(n_PandA/CE_cells, digits = 4) #density of points with real data within the climate envelope
      
      #creates raster of all cells outside CE
      neg_sp_CE = sp_CE == 0 
      
      #calculate desired number of pseudo absence points outside CE based on PandA density within CE
      neg_CE_cells = sum(as.matrix(neg_sp_CE), na.rm=T) #calculates no. cells outside CE
      n_PseudoAbs_pts = round(neg_CE_cells*CE_point_density)+Act_abs #accounts for the actual absences for calculating pseudo absences      
      PseudoAbs_cand_pts = rasterToPoints(neg_sp_CE, fun=function(x){x==1}) #Creates matrix of candidate points (x,y,layer)
      
      plot(mySREresp) #plots all cells with data
      plot(sp_CE) #plots climate envelope
      plot(neg_sp_CE) #plots areas outside climate envelope

      # next section assigns pseudo absences anywhere (not limited to CE)       
    }else{
      neg_mySREresp = mySREresp == 0 #creates raster of areas outside those with known data (presence and absence)
      plot(neg_mySREresp) #plots raster outside known data (presence and absence)
      PseudoAbs_cand_pts = rasterToPoints(neg_mySREresp, fun=function(x){x==1}) # Creates matrix of candidate pseudo absence points (x,y,layer)
      n_PseudoAbs_pts = PA.nb.absences + Act_abs #assigns number of pseudo absence points as indicated in config code, accounting for # actual absences
    }  
    PseudoAbs_cand_pts = as.data.frame(PseudoAbs_cand_pts[,1:2]) #extracts only the geographic information for the candidate pseudo absence points
    head(PseudoAbs_cand_pts) #checks the header for the pseudo absence 
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
    
    dups3 <- duplicated(XY_PresAbsPA_Bioclim_sort[, c('cells')]) # Identifies duplicates in cell column 
    n_dups=length(dups3[dups3==TRUE])
    cat('\n','out of ', length(dups3), "points, ",n_dups, "were removed because they were within the same raster cell for", sp_nm)
    mySpeciesOcc_w_Pseudo<-XY_PresAbsPA_Bioclim_sort[!dups3, ] 
    #n_PandA=dim(XY_PresAbsPA_Bioclim_sort)[1]
    
    mySpeciesOcc_w_Pseudo<-mySpeciesOcc_w_Pseudo[,-4] # This drops the cell column from the data frame
    
    head(mySpeciesOcc_w_Pseudo)
    tail(mySpeciesOcc_w_Pseudo)
    
    ###not in FWS code (points map)
    jpeg_name2=paste(sp_nm,"_loc_data_used.jpg", sep = "")
    jpeg(jpeg_name2,
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(seq((min(mySpeciesOcc_w_Pseudo[,1])-0.1),(max(mySpeciesOcc_w_Pseudo[,1])+0.1),by=((max(mySpeciesOcc_w_Pseudo[,1])+0.1)-(min(mySpeciesOcc_w_Pseudo[,1])-0.1))/5), 
         seq((min(mySpeciesOcc_w_Pseudo[,2])-0.1),(max(mySpeciesOcc_w_Pseudo[,2])+0.1),by=((max(mySpeciesOcc_w_Pseudo[,2])+0.1)-(min(mySpeciesOcc_w_Pseudo[,2])-0.1))/5), 
         type = "n", xlab="Lon", ylab="Lat")# setting up coord. system
    points(x=mySpeciesOcc_w_Pseudo[mySpeciesOcc_w_Pseudo[,3]=='NA',1], y=mySpeciesOcc_w_Pseudo[mySpeciesOcc_w_Pseudo[,3]=='NA',2], type = "p", col = "grey", pch=20,cex = 0.7)
    points(x=mySpeciesOcc_w_Pseudo[mySpeciesOcc_w_Pseudo[,3]==0,1], y=mySpeciesOcc_w_Pseudo[mySpeciesOcc_w_Pseudo[,3]==0,2], type = "p", col = "red", pch=20,cex = 0.7)
    points(x=mySpeciesOcc_w_Pseudo[mySpeciesOcc_w_Pseudo[,3]==1,1], y=mySpeciesOcc_w_Pseudo[mySpeciesOcc_w_Pseudo[,3]==1,2], type = "p", col = "blue", pch=20,cex = 0.7)

    dev.off()
    
    
    ###defining the variables used by biomod2
    cat('\n','biomod model config...')
    myRespName = sp_nm # Insert Species Name Here
    myRespXY = mySpeciesOcc_w_Pseudo[,1:2]
    myResp<-data.frame(Sp_Bio=mySpeciesOcc_w_Pseudo[,3])
    myResp[myResp=='NA']=NA
    #unique(myResp)
    #head(myResp)
    
    
    jnk=dim(mySpeciesOcc_w_Pseudo)[2]
    myBiomodData <- BIOMOD_FormatingData(
      resp.var = myResp,
      expl.var = mySpeciesOcc_w_Pseudo[,4:jnk], # Modify based on number of variables 
      resp.xy = myRespXY,
      resp.name = myRespName,
      PA.nb.rep=PA.nb.rep,
      PA.nb.absences = n_PseudoAbs_pts,
      PA.strategy = PA.strategy,
      PA.dist.min = PA.dist.min)
    #This plotting methods takes way too long!!!  (but it is useful since it plots PAs selected)
    if (plot_graphs==1 & PA.nb.rep<9){  
    jpeg_name3=paste(sp_nm,"_loc_data_used2.jpg", sep = "")
    jpeg(jpeg_name3,
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(myBiomodData)
    dev.off()
    }
    
    memory.limit(size=4095)
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
    
    rm("predictors", "xybackg", "PseudoAbs_cand_pts", "dups2", "jnk", "jnk1", "jnk2") 
    
    cat('\n','fitting...')
    
    myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
                                        models = models_to_run, models.options = myBiomodOption,
                                        NbRunEval=NbRunEval,
                                        DataSplit=80,
                                        Yweights=NULL, 
                                        VarImport=10,
                                        do.full.models=T,
                                        models.eval.meth = eval_stats, #c('TSS','ROC', 'KAPPA'),
                                        SaveObj = TRUE,
                                        rescal.all.models = TRUE)
    
    ## Output the biomod models
    myBiomodModelOut
    
    # output model evaluation metrics
    myBiomodModelEval <- getModelsEvaluations(myBiomodModelOut)    
    dimnames(myBiomodModelEval)
    
    # Outputting the validation metrics for all tests
    if ("TSS" %in% eval_stats){
      myBiomodModelEval["TSS","Testing.data",,,]
      Spp_TSS<- data.frame(myBiomodModelEval["TSS","Testing.data",,,])
      FileName<-paste(sp_nm, "_TSS.csv")
      write.table(Spp_TSS, file = FileName, sep=",", col.names=NA)
    }
    
    if ("ROC" %in% eval_stats){
      myBiomodModelEval["ROC","Testing.data",,,]
      Spp_ROC<- data.frame(myBiomodModelEval["ROC","Testing.data",,,])
      FileName<-paste(sp_nm, "_ROC.csv")
      write.table(Spp_ROC, file = FileName, sep=",", col.names=NA)
    }
    if ("KAPAA" %in% eval_stats){
      myBiomodModelEval["KAPPA","Testing.data",,,]
      Spp_KAP<- data.frame(myBiomodModelEval["KAPPA","Testing.data",,,])
      FileName<-paste(sp_nm, "_KAP.csv")
      write.table(Spp_KAP, file = FileName, sep=",", col.names=NA)
    }
    ## getting the variable importance ##
    getModelsVarImport(myBiomodModelOut)
    Spp_VariImp<- data.frame(getModelsVarImport(myBiomodModelOut))
    #FileName<-paste(sp_nm, "_VariImp.csv")
    write.table(Spp_VariImp, file = FileName00, sep=",", col.names=NA)
    
    save.image("temp_workspace1.RData")   #to save workspace
    rm(list=c("sp_nm","local_config_dir", "spp_nm", "models_to_run", "working_dir", 
              "fitting_clim_data_dir", "env_var_files", "csv_dir", "spp_info", "var_name",
              "eval_stats0", "spp_nm0", "clim_surface_to_use", "proj_nm0", "overwrite", 
              "plot_graphs", "local_config_dir","spp_nm", "clim_data_2000", 
              "clim_data_2100", "working_dir", "csv_dir", "eval_stats",  "crop_raster", "necessary_run_data"))      
    save.image(workspace_name)   #save workspace
    load("temp_workspace1.RData")        
    
    
    ptm1=proc.time() - ptm0
    jnk=as.numeric(ptm1[3])
    jnk=jnk/3600
    cat('\n','It took ', jnk, "hours to model", sp_nm)
  }else{
    cat('\n','fitting for ',sp_nm,'already done...')  
  }    
}