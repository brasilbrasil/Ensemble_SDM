###USER CONFIGURATION
#see 0_sdm_config.r file

###START UNDERHOOD
setwd(working_dir) #sets the working directory
base_working_dir=working_dir #?

#Load required libraries
library(biomod2)
library(raster)
library(randomForest)
library(dismo)
library(mda)
library(stringr)


###not in FWS code (copy necessary files)
#this loop copies the necessary data to run the models into the working directory
dirs=list.dirs(necessary_run_data, full.names = FALSE, recursive = TRUE)
for (dir in dirs){
  layers<-list.files(dir, pattern=NULL, full.names=FALSE, include.dirs = FALSE)
  for (layer in layers){
    layer_full_nm=paste(dir,layer, sep="/")
    if (file.info(layer_full_nm)$isdir==FALSE){
      out_dir_nm=str_replace(dir, necessary_run_data, working_dir)
      dir.create(out_dir_nm, showWarnings = FALSE, recursive = TRUE, mode = "0777")
      out_lyr_nm=str_replace(layer_full_nm, necessary_run_data, working_dir)
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

spp_info=read.csv(paste0(csv_dir,'/FB_spp_data.csv'))


var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}

sp_nm=spp_nm[1]
n_abs_removed=c()
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)
  
  sp_nm_temp=str_replace_all(sp_nm,"_", ".")
  sp_dir=paste0(sp_nm_temp,"/") ###not in FWS code (dir creation)
  dir.create(sp_dir, showWarnings = FALSE)
  ##copy the maxent jar file into the species subdirectory
  #file.copy("maxent.jar", paste0(sp_dir,"maxent.jar"), overwrite = TRUE, recursive = TRUE,
  #          copy.mode = TRUE)
  
  cat('\n',sp_nm,'model fitting...')
  FileName00<-paste(sp_nm, "_VariImp.csv") ###not in FWS code (orverwrite capacity)
  if (file.exists(FileName00)==F | overwrite==1){ #check to see if the analysis for this species was already done    
    # Start the clock!
    ptm0 <- proc.time()
    workspace_name=paste(sp_nm,"_FB_modelfitting.RData", sep = "") #set name of file to save all workspace data after model run
    
    #######Loading datasets#######
    
    ##raster_based_env_grid:
    cat('\n','loading rasters...')
    
    sp_index=which(spp_info[,"Species"]==sp_nm)
    raster_res= spp_info[sp_index,"rasterdir"]
    clim_data_dir=fitting_clim_data_dir 
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
    
    jpeg_name=paste(sp_nm,"_env_vars_used.jpg", sep = "")
    jpeg(jpeg_name,
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(predictors, col=rev(terrain.colors(255)), maxpixels=100000, useRaster=FALSE, axes = TRUE, addfun=NULL, Interpolate = TRUE)
    dev.off()
    
    ####species point data
    cat('\n','loading species data...')
    mySpeciesOcc=read.csv(paste(csv_dir,sp_nm,'_pres_abs.csv', sep = "")) #FB_data_points4_PAandA
    
    #presence (and absence) data handling)
    mySpeciesOcc=cbind(mySpeciesOcc[,2:3],pa=mySpeciesOcc[,1])
    head(mySpeciesOcc)
    
    
    ##pseudo-absence handling
    cat('\n','defining candidate PA points...')
    if (PAs_outside_CE){
      P_and_A=mySpeciesOcc[,1:2]
      mySREresp <- reclassify(subset(predictors,1,drop=TRUE), c(-Inf,Inf,0))
      mySREresp[cellFromXY(mySREresp,P_and_A)] <- 1
      sp_CE=sre(Response=mySREresp,Explanatory=predictors,NewData=predictors,Quant=0.025)
      #calculate density of points within sre
      n_PandA=sum(as.matrix(mySREresp), na.rm=T)
      CE_cells=sum(as.matrix(sp_CE), na.rm=T)
      CE_point_density=round(CE_cells/n_PandA)
      
      #create raster outside CE
      neg_sp_CE=sp_CE==0
      
      #calculate desired number of PA points based on PandA density within CE
      neg_CE_cells=sum(as.matrix(neg_sp_CE), na.rm=T)
      jnk=dim(mySpeciesOcc[mySpeciesOcc$pa==0,])[1]
      n_PA_points=round(neg_CE_cells/CE_point_density)+jnk
      PA_candidate_points=rasterToPoints(neg_sp_CE, fun=function(x){x==1})
      
      plot(mySREresp)
      plot(sp_CE)
      plot(neg_sp_CE)
            
    }else{
      Ps=mySpeciesOcc[,1:2]
      mySREresp <- reclassify(subset(predictors,1,drop=TRUE), c(-Inf,Inf,0))
      mySREresp[cellFromXY(mySREresp,Ps)] <- 1
      mySREresp=mySREresp==0
      plot(mySREresp)
      PA_candidate_points=rasterToPoints(mySREresp, fun=function(x){x==1})
      n_PA_points=PA.nb.absences
    }  
    PA_candidate_points=as.data.frame(PA_candidate_points[,1:2])
    head(PA_candidate_points)
    dim(PA_candidate_points)
    PA_candidate_points_noNA=PA_candidate_points[complete.cases(PA_candidate_points),] #removes rows with NAs
    PA_candidate_points_noNA=cbind(PA_candidate_points_noNA,pa=rep('NA', dim(PA_candidate_points_noNA)[1],1))
    names(PA_candidate_points_noNA)=c('X', 'Y', 'pa') 
    head(PA_candidate_points_noNA)
    
    #merge data with pseudoabsence
    mySpeciesOcc<-data.frame(rbind(mySpeciesOcc, PA_candidate_points_noNA))
    
    #### EXTRACTION OF ENV DATA FOR POINT DATA
    cat('\n','extracting env vars to points...')
    XY_pres_extr<-extract(predictors, mySpeciesOcc[,1:2], cellnumbers=T) ###NEW:This creates a new column call "cell" with the cell numbers from the rasterstack ) 
    XY_pres_extr=data.frame(cbind(mySpeciesOcc,XY_pres_extr)) ###NEW CHANGE
    #XY_pres_extr<-cbind(mySpeciesOcc, XY_pres_extr)
    head(XY_pres_extr)
    XY_pres_extrnoNA=XY_pres_extr[complete.cases(XY_pres_extr[4:dim(XY_pres_extr)[2]]),] #removes rows with NAs
    head(XY_pres_extrnoNA) 
    #tail(XY_pres_extrnoNA) 
    
    ### NEW: Select, Count and Remove presence duplicate points in cells 
    jnk=c(1,0, NA)
    jnk=order(match(XY_pres_extrnoNA$pa, jnk))
    #jnk=order(XY_pres_extrnoNA$pa, decreasing=T)
    XY_pres_extrnoNA=XY_pres_extrnoNA[jnk,] #sorting so if duplicates, PA removed before abs, abs removed before Pres
    
    dups3<- duplicated(XY_pres_extrnoNA[, c('cells')]) # Identifies duplicates in cell column 
    n_dups=length(dups3[dups3==TRUE])
    cat('\n','out of ', length(dups3), "points, ",n_dups, "were removed because they were within the same raster cell for", sp_nm)
    mySpeciesOcc<-XY_pres_extrnoNA[!dups3, ] 
    #n_PandA=dim(XY_pres_extrnoNA)[1]
    
    mySpeciesOcc<-mySpeciesOcc[,-4] # This drops the cell column from the data frame
    
    head(mySpeciesOcc)
    tail(mySpeciesOcc)
    
    ###not in FWS code (points map)
    jpeg_name=paste(sp_nm,"_loc_data_used.jpg", sep = "")
    jpeg(jpeg_name,
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(seq((min(mySpeciesOcc[,1])-0.1),(max(mySpeciesOcc[,1])+0.1),by=((max(mySpeciesOcc[,1])+0.1)-(min(mySpeciesOcc[,1])-0.1))/5), 
         seq((min(mySpeciesOcc[,2])-0.1),(max(mySpeciesOcc[,2])+0.1),by=((max(mySpeciesOcc[,2])+0.1)-(min(mySpeciesOcc[,2])-0.1))/5), 
         type = "n", xlab="Lon", ylab="Lat")# setting up coord. system
    points(x=mySpeciesOcc[mySpeciesOcc[,3]=='NA',1], y=mySpeciesOcc[mySpeciesOcc[,3]=='NA',2], type = "p", col = "grey", pch=20,cex = 0.7)
    points(x=mySpeciesOcc[mySpeciesOcc[,3]==0,1], y=mySpeciesOcc[mySpeciesOcc[,3]==0,2], type = "p", col = "red", pch=20,cex = 0.7)
    points(x=mySpeciesOcc[mySpeciesOcc[,3]==1,1], y=mySpeciesOcc[mySpeciesOcc[,3]==1,2], type = "p", col = "blue", pch=20,cex = 0.7)

    dev.off()
    
    
    ###defining the variables used by biomod2
    cat('\n','biomod model config...')
    myRespName = sp_nm # Insert Species Name Here
    myRespXY = mySpeciesOcc[,1:2]
    myResp<-data.frame(Sp_Bio=mySpeciesOcc[,3])
    myResp[myResp=='NA']=NA
    #unique(myResp)
    #head(myResp)
    
    
    jnk=dim(mySpeciesOcc)[2]
    myBiomodData <- BIOMOD_FormatingData(
      resp.var = myResp,
      expl.var = mySpeciesOcc[,4:jnk], # Modify based on number of variables 
      resp.xy = myRespXY,
      resp.name = myRespName,
      PA.nb.rep=PA.nb.rep,
      PA.nb.absences = n_PA_points,
      PA.strategy = PA.strategy,
      PA.dist.min = PA.dist.min)
    #This plotting methods takes way too long!!!  (but it is useful since it plots PAs selected)
    if (plot_graphs==1 & PA.nb.rep<9){  
    jpeg_name=paste(sp_nm,"_loc_data_used2.jpg", sep = "")
    jpeg(jpeg_name,
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
    
    rm("predictors", "xybackg", "PA_candidate_points", "dups2", "jnk", "jnk1", "jnk2") 
    
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