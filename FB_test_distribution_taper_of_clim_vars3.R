
###USER CONFIGURATION
local_config_dir='C:/Users/lfortini/'
#local_config_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/' #'C:/Users/lfortini/'
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Akikiki","Akekee", "Anianiau", "Kauai_Amakihi", "Kauai_Elepaio", "Puaiohi", "Oahu_Amakihi", "Oahu_Elepaio", "Hawaii_Creeper", "Hawaii_Amakihi","Omao", "Hawaii_Akepa","Palila","Hawaii_Elepaio","Akohekohe","Maui_Alauahio","Maui_Parrotbill","Elepaio", "Amakihi","Apapane", "Iiwi")   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
#spp_nm=c("Akekee", "Akikiki", "Kauai_Amakihi")#, "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
server=0
overwrite=1
exclude_abs=FALSE
proj_name="Pres_Abs"
Quant = 0.025
Quant_val=(1-Quant*2)*100

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
  
  FileName<-paste("quantiles/Q", Quant_val, "_", proj_name, "_", sp_nm,"_taper_off_test.csv", sep="")
  if (file.exists(FileName)==F | overwrite==1){
    #######Loading datasets#######
    sp_data=read.csv(paste(csv_dir,sp_nm,'_pres_abs.csv', sep = "")) #FB_data_points4_PAandA
    
    # Select, Count and Remove presence Duplicates
    jnk=dim(sp_data)[1]
    dups2<- duplicated(sp_data[, c('X','Y')])
    sum(dups2)
    sp_data<-sp_data[!dups2, ]
    jnk1=dim(sp_data)[1]
#     if (exclude_abs){
#       sp_data=sp_data[sp_data[,"pa"]==1,] #get rid of absences      
#     }
    jnk2=dim(sp_data)[1]
    pres_total=sum(sp_data$pa)
    abs_total=jnk2-pres_total
    prop_pres=pres_total/jnk2
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
    extrem.cond <- t(apply(as.data.frame(explanatory), 2, quantile, probs = c(0 + Quant, 1 - Quant), 
                           na.rm = TRUE))          
    extrem.cond=as.data.frame(extrem.cond)
    extrem.cond=cbind(species=matrix(sp_nm,dim(extrem.cond)[1],1),extrem.cond)
    extrem.cond=cbind(type=matrix(proj_name,dim(extrem.cond)[1],1),extrem.cond)
    extrem.cond=cbind(var=var_name,extrem.cond)
    row.names(extrem.cond) <- NULL 
    all_var_res=c()
    var=var_name[1]
    for (var in var_name){
      temp_explanatory=explanatory[,var]
      temp_response_var=Response_var[,3]
      jnk_low=extrem.cond[which(extrem.cond$var==var),4]
      jnk_high=extrem.cond[which(extrem.cond$var==var),5]

      ##low extreme presence
      jnk=which(temp_explanatory<jnk_low)
      temp_response_var1=temp_response_var[jnk]
      Low_extr_total=length(temp_response_var1)
      Low_extr_pres=sum(temp_response_var1)
      Low_extr_prop_pres=Low_extr_pres/Low_extr_total

      ##low extreme presence
      jnk=which(temp_explanatory>jnk_high)
      temp_response_var1=temp_response_var[jnk]
      High_extr_total=length(temp_response_var1)
      High_extr_pres=sum(temp_response_var1)
      High_extr_prop_pres=High_extr_pres/High_extr_total

      var_res=c(pres_total, abs_total, prop_pres,Low_extr_total,Low_extr_pres,
        Low_extr_prop_pres, High_extr_total, High_extr_pres, High_extr_prop_pres,
                Low_extr_prop_pres/prop_pres, High_extr_prop_pres/prop_pres)
      all_var_res=rbind(all_var_res,var_res)        
    }
    
    all_var_res=as.data.frame(all_var_res)
    all_var_res=cbind(species=matrix(sp_nm,dim(all_var_res)[1],1),all_var_res)
    all_var_res=cbind(var=var_name,all_var_res)
    var_names=c("Var", "Species", "pres_total", "abs_total", "prop_pres","Low_extr_total",
                "Low_extr_pres", "Low_extr_prop_pres", "High_extr_total", 
                "High_extr_pres", "High_extr_prop_pres", "Prop_pres_in_Low_extreme",
                "Prop_pres_in_High_extreme")
    names(all_var_res)=var_names
    #FileName=
    write.table(all_var_res, file = FileName, sep=",", col.names=NA)
    
  }else{
    cat('\n', sp_nm, " already calculated")
  }
  if (sp_nm==spp_nm[1]){
    all_quantiles=all_var_res
  }else{
    all_quantiles=rbind(all_quantiles,all_var_res)    
  }
}
#names(all_quantiles)[4:5]=c("QL","QU")
FileName<-paste("quantiles/Q", Quant_val, "_", proj_name, "_all_spp_taper_test.csv", sep="")
write.table((all_quantiles), file = FileName, sep=",", col.names=NA)
