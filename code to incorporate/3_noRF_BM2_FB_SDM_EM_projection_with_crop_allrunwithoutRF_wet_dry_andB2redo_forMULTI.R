rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
# install.packages(c("biomod2", "stringr", "colorRamps", "rasterVis", "raster", "dismo"), 
#                  dependencies = c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances"))
options(java.parameters="-Xmx100g", expressions= 500000)

library(biomod2)
library(stringr)
library(colorRamps)
library(rasterVis)
library(raster)
library(dismo)

plot_graphs=1

# !!!!R version 3.0 seems to have a flaw in the path.files function where it cannot handle a "paste" function, whereas earlier versions can
# !!!!Check lines 219-220 and 274 and see if you can modify them for multi-iteration processing.
#-------------------------------------------------------------------------------------------------
# !!!! Getting the following error "protect(): protection stack overflow" is common with this script, the following will only help if you use 
# this code through the command line (e.g. not Rstudio)
# --max-ppsize=500000
#-------------------------------------------------------------------------------------------------

spp_nm=c("Iiwi")

notserver = T # set to "F" if a server is being used to develop this model and "T" if not
models_to_run=c('GBM','MAXENT')#,'RF')
eval_stats=c("ROC", "TSS", "KAPPA")
working_dir='G:/RWorkDir_Multi2/'
memory = T #keep.in.memory=memory
baseline_or_future=1 #1 for baseline, 2 for baseline_wettest, 3 for baseline_driest, 4 for future, 5 for future_wettest, 6 for future_driest
overwrite=1 #if 1, will overwrite past results
eval_stats=c("ROC","TSS","KAPPA") 
maxentWDtmp = paste("maxentWDtmp_", baseline_or_future, sep = "")
memory.limit(size=24000000)
temp<-paste('G:/temp_RWorkDir_Multi2_', baseline_or_future, '/', sep='')
dir.create(temp)
rasterOptions(tmpdir=temp, timer = T, progress = "text", todisk  = T)
file.copy(paste(working_dir, 'maxent.jar', sep = ""), paste(working_dir, 'maxent', baseline_or_future, '.jar', sep = ""), copy.mode = TRUE)
# If an error occurs on a specific Island, remove the islands (up to that Island in the sequence) 
# Do not remove Kahoolawe b/c will removed in the code.
# you can also run iterations of the code using different islands to reduce the stack overflow issue and get the outputs faster
# if this is done you may get another error associated with the base files overlapping writing over oneanother, but that can be fixed 
Islands <- c('Kauai', 'Oahu', 'Molokai','Lanai','Maui','Kahoolawe','Hawaii')#

clim_data_2000="G:/REnviroLayers/mixed_data_2000_250m/"
clim_data_2000wettest="G:/REnviroLayers/mixed_data_2000_250mwettest/"
clim_data_2000driest= "G:/REnviroLayers/mixed_data_2000_250mdriest/"
clim_data_2100="G:/REnviroLayers/mixed_data_2100_250m/"
clim_data_2100wettest="G:/REnviroLayers/mixed_data_2100_250mwettest/"
clim_data_2100driest= "G:/REnviroLayers/mixed_data_2100_250mdriest/"
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

env_var_files=c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") ###DEBUG!!
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")

####START UNDERHOOD
setwd(working_dir)
spp_nm0=spp_nm
eval_stats0=eval_stats

#sp_nm="Akepa" #debug
var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}
memory.limit(size=240000)
#sp_nm=spp_nm[1]
spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))

LoadProjectionManually <- function(bm_proj){
  if(bm_proj@proj@inMemory){
    cat("\n\tprojection already loaded!")
  } else{
    filesToLoad <- list.files(path=sub("/individual_projections","", bm_proj@proj@link), full.names=T)
    toMatch <- c('.grd$','.img$')
    filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)  
    if(length(filesToLoad)){
      bm_proj@proj@val <- raster::stack(filesToLoad[1])
      bm_proj@proj@inMemory=TRUE
      
    } else {
      filesToLoad <- list.files(path=bm_proj@proj@link, full.names=T)
      toMatch <- c('.grd$','.img$')
      filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
      toMatch <- bm_proj@models.projected
      filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
      bm_proj@proj@val <- raster::stack(filesToLoad)
      toMatch <- c(bm_proj@proj@link,".img$",'.grd$', .Platform$file.sep)
      names(bm_proj@proj@val) <- gsub(pattern=paste(toMatch,collapse="|"), "", filesToLoad)
      bm_proj@proj@inMemory=TRUE
    }    
  }
  return(bm_proj)
}

# Modified the basic biomod2 functions such that a secondary 
# temp. file will be made per modelling iteration = allows for multiple 
# iterations of the same code to run consectutively.
#----------------------------------------------------------------------------------------------------------
# See above comment concerning the path.files function


.testnull <-
  function(object, Prev = 0.5 , dat){
    
    if( is.finite(object$deviance) & is.finite(object$null.deviance)){
      if(object$deviance != object$null.deviance){
        if(inherits(dat,'Raster')){
          pred <- predict(dat, model=object, type='response')
        } else{
          pred <- predict(object, dat, type="response")
        }
      }
    }
    
    if(!exists('pred')){
      if(inherits(dat,'Raster')){
        pred <- subset(dat,1,drop=TRUE)
        if(Prev < 0.5) pred <- reclassify(x=pred, rcl=c(-Inf,Inf,0))
        if(Prev >= 0.5) pred <- reclassify(x=pred, rcl=c(-Inf,Inf,1))
      } else{
        if(Prev < 0.5) pred <- rep(0, nrow(dat))
        if(Prev >= 0.5) pred <- rep(1, nrow(dat))      
      }      
    }    
    return(pred)
  }

setGeneric(".Prepare.Maxent.Proj.WorkDir", 
           def = function(Data, proj.name, ...){
             standardGeneric( ".Prepare.Maxent.Proj.WorkDir" )
           } )


setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='data.frame'),
          def = function(Data, xy, species.name =".", proj.name=".", silent=FALSE){
            if(!silent) cat('\n\t\tCreating Maxent Temp Proj Data...')
            if(is.null(xy)) xy <- matrix(1,nrow=nrow(Data), ncol=2, dimnames=list(NULL, c("X","Y")))
            
            #             if(is.null(proj.name))proj.name <- format(Sys.time(), "%s")
            dir.create(file.path(species.name,proj.name,maxentWDtmp,'Proj'), showWarnings=FALSE, recursive=TRUE)
            
            # Proj Data
            Proj_swd <- cbind(rep("proj",nrow(xy)),xy,Data)
            colnames(Proj_swd)  <- c("proj","X","Y",colnames(Data))
            write.table(Proj_swd, file=file.path(species.name,proj.name,maxentWDtmp, 'Proj_swd.csv'), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
          })


setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='RasterStack'),
          def = function(Data, species.name =".",proj.name=".", silent=FALSE){
            if(!silent) cat('\n\t\tCreating Maxent Temp Proj Data...')
            
            #             if(is.null(proj.name))proj.name <- colnames(Data)[1]
            dir.create(file.path(species.name,proj.name,maxentWDtmp,'Proj'), showWarnings=FALSE, recursive=TRUE)
            
            # Proj Data
            for(l in names(Data)){
              if(! file.exists(file.path(species.name,proj.name,maxentWDtmp,'Proj',paste(l,'.asc',sep='')))){
                if(!silent) cat("\n\t\t\t>",l ,"\t:\t" )
                if(grepl(".asc", filename(raster::subset(Data,l,drop=TRUE)) ) ){
                  if(!silent) cat("coping ascii file")
                  file.copy(filename(raster::subset(Data,l,drop=TRUE)), file.path(species.name,proj.name,maxentWDtmp, 'Proj' ,paste(l,'.asc',sep='')))
                } else{
                  if(!silent) cat("creating ascii file")
                  writeRaster(raster::subset(Data,l,drop=TRUE), filename=file.path(species.name,proj.name,maxentWDtmp, 'Proj' ,paste(l,'.asc',sep='')),
                              format='ascii', overwrite=TRUE)        
                }                
              } else{
                if(!silent) cat("\n", file.path(species.name,proj.name,maxentWDtmp,'', paste(l,'.asc',sep='')),'already created !')
              }              
            }
          })

setClass('MAXENT_biomod2_model',
         representation(model_output_dir = 'character'),
         contains = 'biomod2_model',
         prototype(model_class = 'MAXENT'),
         validity = function(object){
           # check model class
           #            if(sum(! ( c("randomForest.formula", "randomForest") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'MAXENT_biomod2_model'),
          function(object, newdata, ...){            
            args <- list(...)            
            if(inherits(newdata, 'Raster')){            
              return(.predict.MAXENT_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.MAXENT_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }          
          })


.predict.MAXENT_biomod2_model.RasterStack <- function(object, newdata,  ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  rm_tmp_files <- args$rm_tmp_files
  temp_workdir <- args$temp_workdir  
  if (is.null(temp_workdir)) temp_workdir <- maxentWDtmp
  if (is.null(rm_tmp_files)) rm_tmp_files <- TRUE
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE  
  .Prepare.Maxent.Proj.WorkDir(Data = newdata, proj.name = file.path(object@resp_name,temp_workdir))  
  cat("\n\t\tRunning Maxent...")
 # for some reason using a paste (e.g. maxent1) in the file.path does not work
  maxent1 <-paste("maxent", baseline_or_future, ".jar", sep='')
  system(command=paste("java -cp ", file.path(object@model_options$path_to_maxent.jar, "maxent.jar"), 
                       " density.Project \"", 
                       file.path(object@model_output_dir, sub("_MAXENT",".lambdas",object@model_name, fixed=T)),"\" ",
                       file.path(object@resp_name, temp_workdir, maxentWDtmp,"Proj"), " ",
                       file.path(object@resp_name, temp_workdir, maxentWDtmp, "projMaxent.grd") , 
                       " doclamp=false visible=false autorun nowarnings notooltips", sep=""), wait = TRUE)   
  cat("\n\t\tReading Maxent outputs...")
  proj <- raster(file.path(object@resp_name, temp_workdir , maxentWDtmp,"projMaxent.grd"))  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }  
  if(on_0_1000) proj <- round(proj*1000)  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  } else if(!inMemory(proj)){
    proj <- readAll(proj) # to prevent from tmp files removing
  }  
  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
      unlink(x=file.path(object@resp_name, temp_workdir), recursive=TRUE, force=TRUE )
    }
  }  
  return(proj)
}

.predict.MAXENT_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  temp_workdir <- args$temp_workdir
  rm_tmp_files <- args$rm_tmp_files
  xy <- args$xy 
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(temp_workdir)) temp_workdir <- maxentWDtmp
  if (is.null(rm_tmp_files)) rm_tmp_files <- TRUE
  
  #   if( is.null(xy) ){
  #     if( sum(c('x','y') %in% colnames(newdata) ) == 2 ){
  #       coor_col <- c( which(colnames(newdata) == 'x'), which(colnames(newdata) == 'y') )
  #       xy <- newdata[,coor_col]
  #       newdata <- newdata[,- coor_col]
  #     } else { 
  #       xy <- data.frame(x=rep(0,nrow(newdata)), y=rep(0,nrow(newdata)))
  #     }
  #   }
  
  ## no xy needed for models projections
  xy <- NULL
  .Prepare.Maxent.Proj.WorkDir(Data = as.data.frame(newdata), xy = xy , proj.name = file.path(object@resp_name,temp_workdir)) 
  cat("\n\t\tRunning Maxent...")
  system(command=paste("java -cp ", file.path(object@model_options$path_to_maxent.jar, "maxent.jar"),
                       " density.Project \"", 
                       file.path(object@model_output_dir, sub("_MAXENT",".lambdas",object@model_name, fixed=T)),"\" ",
                       file.path(object@resp_name, temp_workdir, maxentWDtmp,"Proj_swd.csv"), " ",
                       file.path(object@resp_name, temp_workdir, maxentWDtmp, "projMaxent.asc")," doclamp=false", sep=""), wait = TRUE)   
  cat("\n\t\tReading Maxent outputs...")
  proj <- as.numeric(read.asciigrid(file.path(object@resp_name, temp_workdir , maxentWDtmp, "projMaxent.asc"))@data[,1])  
  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
      unlink(file.path(object@resp_name, temp_workdir),recursive=TRUE,force=TRUE)
    }
  }  
  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }  
  if(on_0_1000) proj <- round(proj*1000) 
  return(proj)  
}

for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('/n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  setwd(working_dir)
  workspace_name=paste(sp_nm,"_FB_EM_fit.RData", sep = "") #set name of file to load workspace data from model run
  if (file.exists(workspace_name)==F){ 
    workspace_name=paste(sp_nm,"_FB_run.RData", sep = "")}
  load(workspace_name)
  plots=paste(working_dir,"AllEMplots_pmw/",  sep="")
  
  if (file.exists(plots)==F | overwrite==1){
    dir.create(plots)}
  
  #model run specific variables that must not be saved to workspace
  spp_nm=spp_nm0
  #sp_nm0=sp_nm
  eval_stats=eval_stats0  
#   clim_data_dir0=paste(clim_surface_to_use,raster_res,"/grd/",sep="")
  
  proj_nm=proj_nm0 
  
  sp_nm0=str_replace_all(sp_nm,"_", ".")
  workspace_name_out=paste(sp_nm0,"_FB_EM_proj_", proj_nm, ".RData", sep = "")
  #jnk=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.TSS", sep = "")
  #jnk2=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.TSS", sep = "")
  
  predictors = raster(paste(clim_surface_to_use, "Apapane_Iiwi/grd/", env_var_files[1], sep=""))
  for (jj in 2:jnk0){
    temp=raster(paste(clim_surface_to_use, "Apapane_Iiwi/grd/", env_var_files[jj], sep=""))
    predictors = addLayer(predictors, temp)
  }
  names(predictors)<- var_name
  
  # Defineing the extent of the different islands each Trupanea lives on 
  Kauai = c(-159.82,-159.26, 21.84, 22.25)
  Oahu = c(-158.32, -157.62,  21.22, 21.73)
  Molokai = c(-157.34, -156.69, 21.03, 21.25)
  Lanai = c(-157.08, -156.78, 20.70, 20.92)
  Maui= c(-156.8, -155.53, 20.46, 21.05)
  Hawaii = c(-156.10,-154.74, 18.87, 20.30)
  Kahoolawe = c(-156.8, -156.51, 20.46, 20.62)
  
  # Cutting out each island
  for (i in 1:length(Islands)){
    e = extent(get(Islands[i]))
    Isras = crop(predictors, e, snap = 'in')
    assign(Islands[i], Isras)
  }
  
  # for Maui I need to cut out Kahoolawe, but because of the extent issues one has to first reclass, merge and then reclass the merged Kahoo to NA
  rcl <- c(0, 10000, -1)
  rcl <- matrix(rcl, ncol=3, byrow=TRUE)
  a=reclassify(Kahoolawe, rcl)  # defining the Kahoolawe raster values as -1 so that they can be distinguished once merged with Maui
  Maui = try(merge(a, Maui), T)
  
  rcl <- c(-1, NA)
  rcl <- matrix(rcl, ncol=2, byrow=TRUE)
  Maui=try(reclassify(Maui, rcl), T) 
  
  sp_row<- which(spp_info[,"Species"]==sp_nm)
  spIsland<- spp_info[sp_row,(6:length(names(spp_info)))]

  jnk<-seq(1:length(names(spp_info))-1)*spIsland
  jnkx<-jnk[jnk>0]
  spIsland2<- names(spp_info)[c(jnkx+5)]

  SpIslandX<- which(spIsland2 %in% Islands[-(which(Islands =="Kahoolawe"))])
  spIsland2<-names(spIsland[SpIslandX])
  workspace_name_out0=paste(sp_nm,"_FB_all_model_proj_", proj_nm, ".RData", sep = "")
  if (file.exists(workspace_name_out0)==F | overwrite==1){  
    for (spIsl in spIsland2){
      predictors<-get(spIsl)
      if (baseline_or_future==2|3|5|6){
        predictors<-stack((subset(predictors, 1)),
                          (subset(predictors, 2)),
                          (subset(predictors, 3))*4,
                          (subset(predictors, 4)))
        names(predictors)<- var_name
      }
      
      cat('/n',sp_nm,'projection raster stack created...')
      gc()
      
      # If the analysis has already been completed to get the graphics run everything from this line above 
      #  and then skip to the graphics/plots seciton (cntrl+alt+b)
      
      myBiomomodProj_baseline <- BIOMOD_Projection(
        modeling.output = myBiomodModelOut,
        new.env = stack(predictors), #error: additional stack fx
        proj.name = proj_nm,
        selected.models = 'all',
        binary.meth = eval_stats,
        compress = 'xz',
        clamping.mask = T, 
        keep.in.memory=memory,
        overwite = T)
      gc()
      cat('/n',sp_nm,'projection complete...')
      cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn')
      save.image("temp_workspace3.RData")   #to save workspace
      #       try(rm(list=c("spp_info", "eval_stats0", "spp_nm0", "clim_surface_to_use", "proj_nm0", "overwrite", 
      #                 "plot_graphs", "local_config_dir","spp_nm", "clim_data_2000", 
      #                 "clim_data_2100", "working_dir", "env_var_files", "csv_dir", "eval_stats")))      
      save.image(workspace_name_out0)   #to save workspace
      load("temp_workspace3.RData")
      
      cat('/n',sp_nm,'projection of individual models loaded from past run...')
      
      sp_nm0=str_replace_all(sp_nm,"_", ".")
      
      for (ii in 1:length(eval_stats)){         
        for (i in 1:length(models_to_run)){
          setwd(working_dir)
          sp_nm0=str_replace_all(sp_nm,"_", ".")
          a <- stack(paste(working_dir, sp_nm0,"/proj_", proj_nm, "/proj_", proj_nm, '_', sp_nm0, "_", eval_stats[ii], "bin.grd", sep = ""))
          
          a<-subset(a, (length(names(a)) - (length(models_to_run)))+i)#:length(names(a)))
          nm<-names(a)
          a1<-a
          rcl <- c(NA, 0)
          rcl <- matrix(rcl, ncol=2, byrow=TRUE)
          a <-reclassify(a, rcl)
          names(a)<-nm
          
          b <- stack(paste(working_dir, sp_nm0,"/proj_", proj_nm, "/proj_", proj_nm,"_", sp_nm0, ".grd", sep = ""))
          b<-subset(b, (length(names(b)) - (length(models_to_run)))+i)
          
          b<-subset(b, length(names(b)))
          c<-a*b
          rcl <- c(0, NA)
          rcl <- matrix(rcl, ncol=2, byrow=TRUE)
          c<-reclassify(c, rcl)
          names(c)<-nm
          
          assign(paste(eval_stats[ii], '_', models_to_run[i], "_scaled_andbinnedEM_pmw", sep = ""), c)
          assign(paste(eval_stats[ii], '_', models_to_run[i], "_binnedEM_pmw", sep = ""), a1)
        }
        
        # If there is a problem with this part of the code it may be due to the names section = added and would take some time to test
        
        c2<-stack(get(paste(eval_stats[ii], '_', models_to_run[1], "_scaled_andbinnedEM_pmw", sep = "")))
        names(c2)<-paste(eval_stats[ii], '_', models_to_run[1], "_scaled_andbinnedEM_pmw", sep = "")
        for (i in 2:length(models_to_run)){
          jnkc<-raster(get(paste(eval_stats[ii], '_', models_to_run[i], "_scaled_andbinnedEM_pmw", sep = "")))
          names(jnkc)<-paste(eval_stats[ii], '_', models_to_run[i], "_scaled_andbinnedEM_pmw", sep = "")
          c2<-addLayer(c2, jnkc)
        }
        
        a2<-stack(get(paste(eval_stats[ii], '_', models_to_run[1], "_binnedEM_pmw", sep= "")))
        names(a2)<-paste(eval_stats[ii], '_', models_to_run[1], "_binnedEM_pmw", sep= "")
        for (i in 2:length(models_to_run)){
          jnkc<-raster(get(paste(eval_stats[ii], '_', models_to_run[i], "_binnedEM_pmw", sep= "")))
          names(jnkc)<-paste(eval_stats[ii], '_', models_to_run[i], "_binnedEM_pmw", sep= "")
          a2<-addLayer(a2, get(paste(eval_stats[ii], '_', models_to_run[i], "_binnedEM_pmw", sep= "")))
        }
        
        writeRaster(c2, paste(sp_nm, "/proj_", proj_nm0, "/", names(c2),
                              spIsl, ".grd", sep = ""), bylayer=T ,overwrite=T) #defualt format is .grd
        writeRaster(a2, paste(sp_nm, "/proj_",proj_nm0, "/",names(a2),
                              spIsl, ".grd", sep = ""), bylayer=T ,overwrite=T) #defualt format is .grd
        
      }
      
      # subsetting the S4 class object 'myBiomomodProj_baseline'  such that it only uses the two main modelling approaches (GBM and Maxent) 
      # that do not overfit         
      myBiomodProjection <- LoadProjectionManually(myBiomomodProj_baseline)
      
      ###################################################
      ### code chunk number 18: EnsembleForecasting_future
      ###################################################
      myBiomodEF <- BIOMOD_EnsembleForecasting(
        projection.output = myBiomodProjection,
        total.consensus = T,
        EM.output = myBiomodEM, 
        binary.meth=eval_stats, 
        keep.in.memory=memory,
        overwrite=T) ###DEBUG###
      cat('/n',sp_nm,'ensemble projection done...')
      cat('point 1 mem', memory.size(), memory.size(max=TRUE), 'nn')
      
      ###################################################
      ### code chunk number 19: EnsembleForecasting_loading_res
      ###################################################
      #eval_stats=eval_stats ###DEBUG###
           
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
                
      ems_a<-stack(get(paste("TotalConsensus_EMScaledandBinnedby_", eval_stats[1], sep="")))
      for (i in 2:length(eval_stats)){
        ems_a<-addLayer(ems_a, get(paste("TotalConsensus_EMScaledandBinnedby_", eval_stats[i], sep="")))
      }          
      ems_b<-stack(get(paste("TotalConsensus_EMBinnedby_", eval_stats[1], sep="")))
      for (i in 2:length(eval_stats)){
        ems_b<-addLayer(ems_b, get(paste("TotalConsensus_EMBinnedby_", eval_stats[i], sep="")))
      }          
      writeRaster(ems_a, paste(sp_nm, "/proj_", proj_nm0, "/", names(ems_a),"_Tot_Consensus_", 
                               spIsl, ".grd", sep = ""), bylayer=T ,overwrite=T) #defualt format is .grd
      writeRaster(ems_b, paste(sp_nm, "/proj_",proj_nm0, "/",names(ems_b),"_Tot_Consensus_", 
                               spIsl, ".grd", sep = ""), bylayer=T ,overwrite=T) #defualt format is .grd
    }

    # If the analysis has already been completed run everything below this point (cntrl + alt + E) for the plots
    #-----------------------------------------------------------------------
    # reading in the rasters from the saved .grd files and merging them
    
    for (jj in 1:length(eval_stats)){     
      island_mod = raster(list.files(paste(working_dir, sp_nm0,"/proj_", proj_nm, "/", sep = ""), 
                                     pattern = paste(eval_stats[jj], '_ef.pmw_Tot_Consensus_', names(spIsland)[1], ".grd", sep=""),
                                     full.names=T)) 
      
      for (spIsla in names(spIsland)){
        temp=raster(list.files(paste(working_dir, sp_nm0,"/proj_", proj_nm, "/", sep = ""), 
                               pattern = paste(eval_stats[jj], '_ef.pmw_Tot_Consensus_', spIsla, ".grd", sep=""),
                               full.names=T))
        island_mod = merge(island_mod, temp)
      }
      assign(paste(eval_stats[jj], "binandscale_Consensus_", sp_nm0, sep=""), island_mod)
    }
    
    for (jj in 1:length(eval_stats)){     
      island_mod = raster(list.files(paste(working_dir, sp_nm0,"/proj_", proj_nm, "/", sep = ""), 
                                     pattern = paste(eval_stats[jj], '_ef.pmw_', eval_stats[jj], 'bin_Tot_Consensus_', names(spIsland)[1], ".grd", sep=""),
                                     full.names=T)) 
      
      for (spIsla in names(spIsland)){
        temp=raster(list.files(paste(working_dir, sp_nm0,"/proj_", proj_nm, "/", sep = ""), 
                               pattern = paste(eval_stats[jj], '_ef.pmw_', eval_stats[jj], 'bin_Tot_Consensus_', spIsla, ".grd", sep=""),
                               full.names=T))
        island_mod = merge(island_mod, temp)
      }
      assign(paste(eval_stats[jj], "binned_Consensus_", sp_nm0, sep=""), island_mod)
    }
    plot(KAPPAbinned_Consensus_Iiwi)
    plot(KAPPAbinandscale_Consensus_Iiwi)
    
    setwd(plots)
    
    #make sure to turn "server" to "T" (e.g. useRaster = T) if Windows server or XP is being used to model/develop these graphics
    
    gc = c('antiquewhite1', 'transparent')
    col5 <- colorRampPalette(c('blue', 'sandybrown', 'darkgreen'))
    
    jpeg_name=paste(proj_nm,"_", sp_nm0,"_TOTALCONSENSUS_Binandscaled_runs_.jpg", sep = "")
    jpeg(jpeg_name, width = 5*length(eval_stats), height = 5, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)  
    par(pin = c(4,4), cex = 1, cex.main = 1, cex.axis = 0.8, mfcol=c(1,length(eval_stats)), mgp = c(1, 0.5, 0),
        mar=c(2, 2, 1.5, 1), oma = c(0, 0, 0, 1), bg = "transparent")
    
    ras1<-get(paste("ROCbinandscale_Consensus_", sp_nm0, sep = ""))
    names(ras1)<-paste("Binned_by_ROC_EM_Consensus_for_", sp_nm0, sep = "")
    plot(ras1, col=col5(255), useRaster = notserver, axes = TRUE, addfun=F, 
             interpolate = TRUE, legend = F, add = F, bg = "transparent", main = names(ras1))  
    ras2<-get(paste("ROCbinned_Consensus_", sp_nm0, sep = ""))
    plot(ras2, col = gc, useRaster = notserver, axes = F, addfun=F, interpolate = TRUE, legend = F, add = T) 
    
    par(mar=c(2, 0, 1.5, 0))
    ras1<-get(paste("TSSbinandscale_Consensus_", sp_nm0, sep = ""))
    names(ras1)<-paste("Binned_by_TSS_EM_Consensus_for_", sp_nm0, sep = "")
    plot(ras1, col = col5(255), useRaster = notserver, axes = TRUE,interpolate = TRUE, 
             legend = F, yaxt = 'n', add = F, bg = "transparent", main = names(ras1))
    ras2<-get(paste("TSSbinned_Consensus_", sp_nm0, sep = ""))
    plot(ras2, col = gc, useRaster = notserver, axes = F, addfun=F, interpolate = TRUE, yaxt = 'n', legend = F, add = T)
    
    par(mar=c(2, 0, 1.5, 3.5))
    ras1<-get(paste("KAPPAbinandscale_Consensus_", sp_nm0, sep = ""))
    names(ras1)<-paste("Binned_by_KAPPA_EM_Consensus_for_", sp_nm0, sep = "")
    plot(ras1, col = col5(255), useRaster = notserver, axes = T, interpolate = F, 
             legend = T, yaxt = 'n', add = F, bg = "transparent", main = names(ras1))
    ras2<-get(paste("KAPPAbinned_Consensus_", sp_nm0, sep = ""))
    plot(ras2, col = gc, useRaster = notserver, axes = F, interpolate = F, legend = F, add = T)
    
    legend("bottomright",legend = c("Absent"), fill = gc[1], cex = 0.8)
    dev.off()
    
    cat('/n',sp_nm,'consensus projection graphs done...')
  }
  save.image(workspace_name_out0)
  removeTmpFiles(h=1)
}

setwd(working_dir)
cat('/n',sp_nm,'ensemble projection figures done...')
save.image("temp_workspace4.RData")   #to save workspace
# try(rm(list=c("spp_info","eval_stats0", "spp_nm0", "clim_surface_to_use", "proj_nm0", "overwrite", 
#               "plot_graphs", "local_config_dir","spp_nm", "clim_data_2000", 
#               "clim_data_2100", "working_dir", "env_var_files", "csv_dir", "eval_stats")))      
#    #to save workspace
#  load("temp_workspace4.RData")}                                 

  # }else{
  #    cat('/n',sp_nm,'previously calculated...')
 