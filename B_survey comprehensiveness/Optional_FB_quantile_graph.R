rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
source(paste0("C:/Users/lfortini/","directory_registry.r"))
spp_nm=c('Akekee', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi', 'Kauai_Amakihi', 'Hawaii_Amakihi', 'Apapane', 'Amakihi', 'Elepaio', 'Iiwi')
server=1
overwrite=1
exclude_abs=F
exclude_pres=T
proj_name="Abs"
Quant = 0.025

Quant_val=(1-Quant*2)*100

if (server==1){
#   working_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/'
#   clim_data_dir0="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/100m/" 
  working_dir='Y:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_less_PAs/'
  clim_data_dir0="Y:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg/bioclims_abs/all_baseline/250m/"
}else{
  working_dir='D:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_220runs/'
  clim_data_dir0="D:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg/bioclims_abs/all_baseline/250m/"
}

env_var_files=c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif") 
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
spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))
dir.create("quantiles/", showWarnings=FALSE)

var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}
