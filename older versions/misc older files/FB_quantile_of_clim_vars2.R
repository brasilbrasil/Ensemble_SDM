rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
local_config_dir='C:/Users/lfortini/'
#local_config_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/' #'C:/Users/lfortini/'
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Akikiki","Akekee", "Anianiau", "Kauai_Amakihi", "Kauai_Elepaio", "Puaiohi", "Oahu_Amakihi", "Oahu_Elepaio", "Hawaii_Creeper", "Hawaii_Amakihi","Omao", "Hawaii_Akepa","Palila","Hawaii_Elepaio","Akohekohe","Maui_Alauahio","Maui_Parrotbill","Elepaio", "Amakihi","Apapane", "Iiwi")   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
#spp_nm=c("Akekee", "Akikiki", "Anianiau")#, "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
server=0
overwrite=1
exclude_abs=FALSE
proj_name="_presabs_95pctQuant"

if (server==1){
  working_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/'
  clim_data_dir0="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/100m/" 
  necessary_run_data='Y:/FB analysis/FB SDM/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/response_curves/'
  clim_data_dir0="C:/Users/lfortini/Data/SDM_env_data/all_grd/all_baseline/100m/"
  necessary_run_data='C:/Users/lfortini/Data/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
}

env_var_files=c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") 
crop_raster_dir=paste(working_dir, 'map_crop/',sep="")
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")

###START UNDERHOOD
setwd(working_dir)
library(biomod2)
library(raster)
library(randomForest)
library(dismo)
library(mda)
library(stringr)

all_quantiles=c()
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
dir.create("quantiles/", showWarnings=FALSE)

var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}

sp_nm = spp_nm[1]
for (sp_nm in spp_nm){  
  sp_nm=as.character(sp_nm)
  cat('\n',sp_nm,'modeling...')
  sp_str=sp_nm
  
  FileName<-paste("quantiles/", sp_nm, proj_name,".csv", sep="")
  if (file.exists(FileName)==F | overwrite==1){
    #######Loading datasets#######
    sp_data=read.csv(paste(csv_dir,sp_nm,'_pres_abs.csv', sep = "")) #FB_data_points4_PAandA
    
    # Select, Count and Remove presence Duplicates
    jnk=dim(sp_data)[1]
    dups2<- duplicated(sp_data[, c('X','Y')])
    sum(dups2)
    sp_data<-sp_data[!dups2, ]
    jnk1=dim(sp_data)[1]
    if (exclude_abs){
      sp_data=sp_data[sp_data[,"pa"]==1,] #get rid of absences      
    }
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
    
    cat('\n',sp_nm,'modeling...')
    sp_data=data.frame(cbind(x=sp_data$X, y=sp_data$Y, pres=sp_data$pa))
    n_points=dim(sp_data)[1]
    Response_var=sp_data
    explanatory=extract(predictors, Response_var[,1:2])
    Quant = 0.025
    extrem.cond <- t(apply(as.data.frame(explanatory), 2, quantile, probs = c(0 + Quant, 1 - Quant), 
                           na.rm = TRUE))          
    extrem.cond=as.data.frame(extrem.cond)
    temp_names=names(extrem.cond)
    for (temp_name in temp_names){
      temp_names[which(temp_names==temp_name)]=paste(temp_name,sp_nm, sep="_")
    }
    names(extrem.cond)=temp_names  
    write.table(extrem.cond, file = FileName, sep=",", col.names=NA)
    
  }else{
    cat('\n', sp_nm, " already calculated")
  }
  if (sp_nm==spp_nm[1]){
    all_quantiles=c(extrem.cond)
  }else{
    all_quantiles=cbind(all_quantiles,extrem.cond)    
  }
}
FileName<-paste("quantiles/all", proj_name,".csv", sep="")
write.table(t(all_quantiles), file = FileName, sep=",", col.names=NA)
