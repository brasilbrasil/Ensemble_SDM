rm(list = ls()) #remove all past worksheet variables
#for each species modeled, have csv of presence data in working directory for the species named speciesname_Ps.csv formated with 3 cols: x,y,pa where pa = 1
#after running the code for whichever many species, copy results (species output folder and workspace file) to a new directory, along with the maxent.jar file
#use the projection code to project the distribution model on different environmental surfaces (do not forget to change the working directory)

###USER CONFIGURATION
local_config_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/' #'C:/Users/lfortini/'
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Akikiki","Akekee", "Anianiau", "Kauai_Amakihi", "Kauai_Elepaio", "Puaiohi", "Oahu_Amakihi", "Oahu_Elepaio", "Hawaii_Creeper", "Hawaii_Amakihi","Omao", "Hawaii_Akepa","Palila","Hawaii_Elepaio","Akohekohe","Maui_Alauahio","Maui_Parrotbill","Elepaio", "Amakihi","Apapane", "Iiwi")   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
spp_nm=c("Akekee")   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
proj_name="_1_bio1_12"
server=1
min_n_points=3
overwrite=1
apply_veg_mask=FALSE

if (server==1){
  working_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/'
  clim_data_dir0="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/100m/" 
  clim_data_dir0f="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_future/100m/" 
  
  necessary_run_data='Y:/FB analysis/FB SDM/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/test/'
  clim_data_dir0="C:/Users/lfortini/Data/SDM_env_data/all_grd/all_baseline/100m/"
  clim_data_dir0f="C:/Users/lfortini/Data/SDM_env_data/all_grd/all_future/100m/"
  necessary_run_data='C:/Users/lfortini/Data/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
}

#env_var_files=c("bio1.grd", "bio12.grd") 
#env_var_files=c("bio5.grd", "bio6.grd", "bio12.grd") 
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

if (apply_veg_mask){
  desc="_masked"
}else{
  desc=""    
}
#sp_nm=spp_nm[1]
for (sp_nm in spp_nm){
  
  sp_nm=as.character(sp_nm)
  cat('\n',sp_nm,'modeling...')
  sp_str=sp_nm
  
  out_raster_name=paste(working_dir, "SREs/", sp_str, proj_name, desc,"_rect_resp_envelopes.tif", sep="")
  if (file.exists(out_raster_name)==F | overwrite==1){
    #######Loading datasets#######
    sp_data=read.csv(paste(csv_dir,sp_nm,'_pres_abs.csv', sep = "")) #FB_data_points4_PAandA
    
    # Select, Count and Remove presence Duplicates
    jnk=dim(sp_data)[1]
    dups2<- duplicated(sp_data[, c('X','Y')])
    sum(dups2)
    sp_data<-sp_data[!dups2, ]
    jnk1=dim(sp_data)[1]
    sp_data=sp_data[sp_data[,"pa"]==1,] #get rid of absences
    jnk2=dim(sp_data)[1]
    head(sp_data)
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
    predictors
    
    
    #future
    clim_data_dir=clim_data_dir0f 
    predictors_fut = raster( paste(clim_data_dir, env_var_files[1], sep=""))
    predictors_fut=crop(predictors_fut,  crop_raster)
    for (jj in 2:jnk0){
      temp=raster(paste(clim_data_dir, env_var_files[jj], sep=""))
      temp=crop(temp,  crop_raster)
      predictors_fut = addLayer(predictors_fut, temp)
    }
    names(predictors_fut)<- var_name
    rm("crop_raster" ,"temp") 
    predictors_fut
    
    
    
    cat('\n',sp_nm,'modeling...')
    sp_data=data.frame(cbind(x=sp_data$X, y=sp_data$Y, pres=sp_data$pa))
    n_points=dim(sp_data)[1]
    if (n_points>=min_n_points){    
      Response_var=sp_data
      explanatory=extract(predictors, Response_var[,1:2])
      pred <- sre(Response_var, explanatory, predictors, Quant = 0.025)
      SRE_raster_present=subset(pred, 3)
      
      tifname=paste("SREs/", sp_str, proj_name,  desc,"_SRE_baseline.tif", sep="")
      writeRaster(SRE_raster_present, tifname, format="GTiff", overwrite=TRUE)
      
      
      #future
      pred <- sre(Response_var, explanatory, predictors_fut, Quant = 0.025)
      SRE_raster_future=subset(pred, 3)
      
      tifname=paste("SREs/", sp_str, proj_name,  desc,"_SRE_future.tif", sep="")
      writeRaster(SRE_raster_future, tifname, format="GTiff", overwrite=TRUE)
      
      
      jnk=SRE_raster_future*10
      BIN_dif=SRE_raster_present+jnk
      m  =  c(9.9,  10.1,  3, 10.9, 11.1, 2)
      rclmat  =  matrix(m,  ncol=3,  byrow=TRUE)
      resp_zone  =  reclassify(BIN_dif,  rclmat)
      
      if (apply_veg_mask){
        habitat = raster( paste(clim_data_dir, "veg_areas.grd", sep=""))
        habitat=crop(habitat, resp_zone)
        #habitat=habitat==0
        resp_zone=resp_zone*habitat
      }
      
      mypalette_numbers=c(0, 1, 2, 3)
      mypalette=c("Grey", "Red", "Green", "Yellow")
      #resp_zone_names0=c("Micro refugia", "Tolerate", "Migrate")
      resp_zone_names0=c("Lost", "Overlap", "Gained")
      
      jnk=unique(resp_zone)
      zones_present=jnk[jnk>0]
      zones_present=zones_present[zones_present<=3]
      resp_zone_colors=mypalette[zones_present+1]
      mypalette_numbers_selected=mypalette[jnk+1] #CHANGED
      resp_zone_names=resp_zone_names0[zones_present]
      
      jpeg_name=paste("SREs/", sp_str, proj_name, desc,"_SRE_baseline.jpg", sep = "")
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(SRE_raster_present,  col=mypalette, legend=F)
      title(paste(sp_nm, " Baseline SRE (n=", n_points, ")", sep=""))
      dev.off()
      
      jpeg_name=paste("SREs/", sp_str, proj_name, desc,"_SRE_future.jpg", sep = "")
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(SRE_raster_future,  col=mypalette, legend=F)
      title(paste(sp_nm, " Future SRE (n=", n_points, ")", sep=""))
      dev.off()
      
      
      jpeg_name=paste("SREs/", sp_str, proj_name, desc,"_response_zones.jpg", sep = "")
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(resp_zone,  col=mypalette_numbers_selected, legend=F)
      #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
      legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
      title(paste(sp_nm, " response zones (n=", n_points, ")", sep=""))
      dev.off()
      
      writeRaster(resp_zone, out_raster_name, format="GTiff", overwrite=TRUE) 
    }
    
  }else{
    cat('\n', sp_nm, " already calculated")
  }
  
}
