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
# iterations of the same code to run consectutively
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
                       file.path(object@resp_name, temp_workdir, maxentWDtmp, "projMaxent.asc") , 
                       " doclamp=false", sep=""), wait = TRUE)  
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
