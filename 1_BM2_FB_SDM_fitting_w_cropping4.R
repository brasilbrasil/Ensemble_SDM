rm(list = ls()) #remove all past worksheet variables
#for each species modeled, have csv of presence data in working directory for the species named speciesname_Ps.csv formated with 3 cols: x,y,pa where pa = 1
#after running the code for whichever many species, copy results (species output folder and workspace file) to a new directory, along with the maxent.jar file
#use the projection code to project the distribution model on different environmental surfaces (do not forget to change the working directory)

###USER CONFIGURATION
#local_config_dir='C:/Users/lfortini/'
local_config_dir='Y:/FB_analysis/FB_SDM/biomod2/' #'C:/Users/lfortini/'
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Kauai_Amakihi")#, "Oahu_Amakihi", "Apapane")# "Akekeke", "Akikiki", "Anianiau", "Kauai_Amakihi", "Kauai_Elepaio", "Puaiohi", "Oahu_Amakihi", "Oahu_Elepaio", "Apapane", "Iiwi",)   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
server=1
remove_PA_abs=TRUE
overwrite=0

models_to_run=c('GBM','RF','MAXENT')
if (server==1){
  working_dir='Y:/FB_analysis/FB_SDM/biomod2/response_curves_500/'
  clim_data_dir0="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/100m/" 
  necessary_run_data='Y:/FB_analysis/FB_SDM/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/test/'
  necessary_run_data='C:/Users/lfortini/Data/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
  clim_data_dir0="C:/Users/lfortini/Data/SDM_env_data/all_grd/all_baseline/100m/"
}

env_var_files=c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") 
crop_raster_dir=paste(working_dir, 'map_crop/',sep="")
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")


###START UNDERHOOD
setwd(working_dir)
library(biomod2)
library(raster)
library(rJava)
library(randomForest)
library(dismo)
library(mda)
library(stringr)

dirs=list.dirs(necessary_run_data, full.names = FALSE, recursive = TRUE)
for (dir in dirs){
  layers<-list.files(dir, pattern=NULL, full.names=FALSE, include.dirs = FALSE)
  for (layer in layers){
    layer_full_nm=paste(dir,layer, sep="/")
    if (file.info(layer_full_nm)$isdir==FALSE){
      out_dir_nm=str_replace(dir, necessary_run_data, working_dir)
      dir.create(out_dir_nm, showWarnings = FALSE, recursive = TRUE, mode = "0777")
      out_lyr_nm=str_replace(layer_full_nm, necessary_run_data, working_dir)
      #out_lyr_nm=str_replace(out_lyr_nm, bl_lr, output_filenm)
      if (file.exists(out_lyr_nm)==F){
        cat('\n','found ', layer, 'in ', dir)
        file.copy(layer_full_nm, out_lyr_nm, overwrite = TRUE, recursive = TRUE,
                  copy.mode = TRUE)
        cat('\n','saved as ', out_lyr_nm)
      }
    }
  }
}
spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))


var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}

#sp_nm=spp_nm[1]
#spp_nm=spp_nm[11:length(spp_nm)]
#spp_nm=spp_nm[1:10]
#spp_nm=rev(spp_nm)
n_abs_removed=c()
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)
  cat('\n',sp_nm,'modeling...')
  FileName00<-paste(sp_nm, "_VariImp.csv")
  if (file.exists(FileName00)==F | overwrite==1){    
    # Start the clock!
    ptm0 <- proc.time()
    workspace_name=paste(sp_nm,"_FB_run.RData", sep = "") #set name of file to save all workspace data after model run
    
    #######Loading datasets#######
    mySpeciesOcc=read.csv(paste(csv_dir,sp_nm,'_pres_abs.csv', sep = "")) #FB_data_points4_PAandA
    
    # Select, Count and Remove presence Duplicates
    jnk=dim(mySpeciesOcc)[1]
    dups2<- duplicated(mySpeciesOcc[, c('X','Y')])
    sum(dups2)
    mySpeciesOcc<-mySpeciesOcc[!dups2, ]
    jnk1=dim(mySpeciesOcc)[1]
    mySpeciesOcc=mySpeciesOcc[mySpeciesOcc[,"pa"]==1,] #get rid of absences
    jnk2=dim(mySpeciesOcc)[1]
    head(mySpeciesOcc)
    cat('\n','removed ', jnk-jnk1, "duplicates for", sp_nm)
    cat('\n','removed ', jnk1-jnk2, "absence records for", sp_nm)
    
    ##raster_based_env_grid:
    sp_index=which(spp_info[,"Species"]==sp_nm)
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
    
    # Ploting predictors may take a substantial amount of time, depending on the file type (i.e. ascii or grd) and resolution)
    # Irritatingly, if the resolution is too high they may not plot...  this also may be a function of the plotting package used
    # i.e. whether 'useRaster' is true or false.  On windows 2008 it doesn't work under 'TRUE'.
    
    jpeg_name=paste(sp_nm,"_env_vars_used.jpg", sep = "")
    jpeg(jpeg_name,
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(predictors, col=rev(terrain.colors(255)), maxpixels=100000, useRaster=FALSE, axes = TRUE, addfun=NULL, Interpolate = TRUE)
    dev.off()
    
    ####Generate 10000 random background pts with good env data
    xybackg<-randomPoints(predictors, n=20000) # Creates 10,000 background/absence points
    colnames(xybackg)=c('X', 'Y')
    XYabackg <- c(rep(0, nrow(xybackg)))
    XYabackg <- data.frame(cbind(xybackg, pa=XYabackg))
    head(XYabackg)
    
    XYabackg_extr<-extract(predictors, XYabackg[,1:2], cellnumbers=T)
    XYabackg_extr<-cbind(XYabackg, XYabackg_extr)
    head(XYabackg_extr)
    dim(XYabackg_extr)
    XYabackg_extrnoNA=XYabackg_extr[complete.cases(XYabackg_extr),] #removes rows with NAs
    
    jnk=nrow(XYabackg_extrnoNA)
    if (jnk>10000){
      XYabackg_extrnoNA=XYabackg_extrnoNA[1:10000,]
    }else{
      cat('\n','could only generate', jnk, "random background points for", sp_nm)
    } #pick only 10k good points, will give error if not enough good points area available
    head(XYabackg_extrnoNA)
    tail(XYabackg_extrnoNA)
    
    #### EXTRACTION OF ENV DATA FOR PRESENCE DATA
    XY_pres_extr<-extract(predictors, mySpeciesOcc[,2:3], cellnumbers=T) ###NEW:This creates a new column call "cell" with the cell numbers from the rasterstack ) 
    XY_pres_extr=data.frame(cbind(mySpeciesOcc[,2:3], pa= c(rep(1, nrow(mySpeciesOcc))),XY_pres_extr)) ###NEW CHANGE
    #XY_pres_extr<-cbind(mySpeciesOcc, XY_pres_extr)
    head(XY_pres_extr)
    XY_pres_extrnoNA=XY_pres_extr[complete.cases(XY_pres_extr),] #removes rows with NAs
    head(XY_pres_extrnoNA) 
    tail(XY_pres_extrnoNA) 
    
    ### NEW: Select, Count and Remove presence duplicate points in cells 
    dups3<- duplicated(XY_pres_extrnoNA[, 'cells']) # Identifies duplicates in cell column 
    #sum(dups3) 
    XY_pres_extrnoNA<-XY_pres_extrnoNA[!dups3, ] 
    #XY_pres_extrnoNA<-XY_pres_extrnoNA[,-4] # This drops the cell column from the data frame
    
    ####combining the presence and pseudoabsence background points
    if (remove_PA_abs){
      jnk= XYabackg_extrnoNA[, 'cells'] %in% XY_pres_extrnoNA[, 'cells']  ####DEBUG DEBUG DEBUG
      XYabackg_extrnoNA=XYabackg_extrnoNA[jnk=='FALSE',]
      jnk=length(which(jnk==TRUE))
      cat('\n','removed', jnk, "random background points that overlaped with presence for", sp_nm)
      n_abs_removed=cbind(n_abs_removed,jnk)
    }
    mySpeciesOcc<-data.frame(rbind(XY_pres_extrnoNA, XYabackg_extrnoNA))
    mySpeciesOcc<-mySpeciesOcc[,-4] # This drops the cell column from the data frame
    
    head(mySpeciesOcc)
    tail(mySpeciesOcc)
    
    jpeg_name=paste(sp_nm,"_loc_data_used.jpg", sep = "")
    jpeg(jpeg_name,
         width = 10, height = 10, units = "in",pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(seq((min(mySpeciesOcc[,1])-0.1),(max(mySpeciesOcc[,1])+0.1),by=((max(mySpeciesOcc[,1])+0.1)-(min(mySpeciesOcc[,1])-0.1))/5), 
         seq((min(mySpeciesOcc[,2])-0.1),(max(mySpeciesOcc[,2])+0.1),by=((max(mySpeciesOcc[,2])+0.1)-(min(mySpeciesOcc[,2])-0.1))/5), 
         type = "n", xlab="Lon", ylab="Lat")# setting up coord. system
    points(x=mySpeciesOcc[mySpeciesOcc[,3]==0,1], y=mySpeciesOcc[mySpeciesOcc[,3]==0,2], type = "p", col = "red", pch=20,cex = 0.7)
    points(x=mySpeciesOcc[mySpeciesOcc[,3]==1,1], y=mySpeciesOcc[mySpeciesOcc[,3]==1,2], type = "p", col = "blue", pch=20,cex = 0.7)
    #a=rasterFromXYZ(mySpeciesOcc[,1:3])
    #plot(a,  title=title_temp)
    dev.off()
    
    
    ###defining the variables used by biomod2
    myRespName = sp_nm # Insert Species Name Here
    myRespXY = mySpeciesOcc[,1:2]
    myResp<-data.frame(Sp_Bio=mySpeciesOcc[,3])
    head(myResp)
    
    jnk0=length(env_var_files)
    jnk=4+jnk0-1
    myBiomodData <- BIOMOD_FormatingData(
      resp.var = myResp,
      expl.var = mySpeciesOcc[,4:jnk], # Modify based on number of variables 
      resp.xy = myRespXY,
      resp.name = myRespName,
      PA.nb.rep = 0)
    
    
    memory.limit(size=4095)
    myBiomodOption <- BIOMOD_ModelingOptions(
      GBM = list( distribution = 'bernoulli', interaction.depth = 7,  shrinkage = 0.001, bag.fraction = 0.5, train.fraction = 1, n.trees = 100,
                  cv.folds = 10),
      MARS = list( degree = 2,
                   penalty = 2,
                   thresh = 0.001,
                   prune = TRUE),
      RF = list(do.classif = TRUE, ntree = 100, mtry = 'default'), 
      MAXENT = list(maximumiterations = 100, visible = TRUE, linear = TRUE, quadratic = TRUE,
                    product = TRUE, threshold = TRUE, hinge = TRUE, lq2lqptthreshold = 80, l2lqthreshold = 10,
                    hingethreshold = 15, beta_threshold = -1, beta_categorical = -1, beta_lqp = -1, 
                    beta_hinge = -1,defaultprevalence = 0.5)
    )
    
    rm("predictors", "xybackg", "XYabackg_extr", "dups2", "jnk", "jnk1", "jnk2") 
    
    # attempting to change the java and R default parameters to 12 gb of ram and to
    # omit rows with NA
    # options(java.parameters = "-Xmx12g" )  # Modify this based on ram available
    # options(na.action=na.omit)  ## Is this redundant. I thought NAs were removed.
    
    ## Modelling ##
    myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
                                        models = models_to_run, models.options = myBiomodOption,
                                        NbRunEval=500,
                                        DataSplit=80,
                                        Yweights=NULL, 
                                        VarImport=10,
                                        models.eval.meth = 'ROC', #c('TSS','ROC', 'KAPPA'),
                                        SaveObj = TRUE,
                                        rescal.all.models = TRUE)
    
    ## Output the biomod models
    myBiomodModelOut
    
    # output model evaluation metrics
    myBiomodModelEval <- getModelsEvaluations(myBiomodModelOut)    
    dimnames(myBiomodModelEval)
    
    #     # Outputting the validation metrics for all tests
    #     myBiomodModelEval["TSS","Testing.data",,,]
    #     Spp_TSS<- data.frame(myBiomodModelEval["TSS","Testing.data",,,])
    #     FileName<-paste(sp_nm, "_TSS.csv")
    #     write.table(Spp_TSS, file = FileName, sep=",", col.names=NA)
    
    myBiomodModelEval["ROC","Testing.data",,,]
    Spp_ROC<- data.frame(myBiomodModelEval["ROC","Testing.data",,,])
    FileName<-paste(sp_nm, "_ROC.csv")
    write.table(Spp_ROC, file = FileName, sep=",", col.names=NA)
    
    #     myBiomodModelEval["KAPPA","Testing.data",,,]
    #     Spp_KAP<- data.frame(myBiomodModelEval["KAPPA","Testing.data",,,])
    #     FileName<-paste(sp_nm, "_KAP.csv")
    #     write.table(Spp_KAP, file = FileName, sep=",", col.names=NA)
    
    ## getting the variable importance ##
    getModelsVarImport(myBiomodModelOut)
    Spp_VariImp<- data.frame(getModelsVarImport(myBiomodModelOut))
    #FileName<-paste(sp_nm, "_VariImp.csv")
    write.table(Spp_VariImp, file = FileName00, sep=",", col.names=NA)
    
    save.image("temp_workspace1.RData")   #to save workspace
    rm(list=c("sp_nm","local_config_dir", "spp_nm", "models_to_run", "working_dir", 
              "clim_data_dir0", "env_var_files", "csv_dir", "spp_info", "var_name",
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