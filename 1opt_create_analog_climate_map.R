rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
plot_graphs=1
source(paste0("Y:/PICCC_analysis/code/","directory_registry.r"))
local_config_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m/') #'C:/Users/lfortini/'
#spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Anianiau", "Kauai_Amakihi", "Hawaii_Elepaio", "Palila")
baseline_or_future=1 #0 for baseline, 1 for future
overwrite=1 #if 1, will overwrite past results
server=1
env_var_files=c("bio1.grd", "bio7.grd", "bio12.grd", "bio15.grd") 

if (server==1){
  working_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m_rounded/')
  clim_data_2000=paste0(DR_FB_clim_data,"all_grd/all_baseline/500m_test/")
  clim_data_2100=paste0(DR_FB_clim_data,"all_grd/all_future/500m/")
  clim_data_2000wettest="D:/GIS_Data/REnviroLayers/mixed_data_2000_250mwettest/"
  clim_data_2000driest= "D:/GIS_Data/REnviroLayers/mixed_data_2000_250mdriest/"
  clim_data_2100wettest="D:/GIS_Data/REnviroLayers/mixed_data_2100_250mwettest/"
  clim_data_2100driest= "D:/GIS_Data/REnviroLayers/mixed_data_2100_250mdriest/"  
}else{
  working_dir='C:/Users/lfortini/Forest bird SDM/biomod2/'
  clim_data_2000="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/250m/"
  clim_data_2100="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_future/250m/"
}
####START UNDERHOOD
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
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")

setwd(working_dir)
library(biomod2)
library(stringr)

#sp_nm="Akepa" #debug
var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}
memory.limit(size=240000)
#sp_nm=spp_nm[1]
spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))


####create analog climate map

all_clim_data=c("clim_data_2000", "clim_surface_to_use")
clim_stacks=c("biovars2000", "biovars2100")
crop_raster=raster(paste(crop_raster_dir,raster_res,".grd",sep=""))
for (dfd in 1:length(all_clim_data)){
  clim_data=all_clim_data[dfd]
  clim_data_dir=get(clim_data)
  predictors_temp = raster( paste(clim_data_dir, env_var_files[1], sep=""))
  predictors_temp=crop(predictors_temp,  crop_raster)
  for (jj in 2:jnk0){
    temp=raster(paste(clim_data_dir, env_var_files[jj], sep=""))
    temp=crop(temp,  crop_raster)
    predictors_temp = addLayer(predictors_temp, temp)
  }
  names(predictors_temp)<- var_name
  assign(clim_stacks[dfd], predictors_temp)
}

jnk=raster(paste(clim_data_dir, "bio12.grd", sep=""))
jnk=crop(jnk,  crop_raster)
Response_var  =  reclassify(jnk,  c(0.1,Inf,1))

pred <- sre(Response_var, biovars2000, biovars2100, Quant = 0.001)
analog_climates2100=subset(pred, 1)
rm("crop_raster" ,"temp")         

raster_name=paste(sp_nm0,"_analog_climates_",proj_nm0, sep="")
tif_name=paste(raster_name,".tif", sep="")
jpeg_name=paste(raster_name,".jpg", sep="")
writeRaster(analog_climates2100, tif_name, format="GTiff", overwrite=TRUE)

plot_vars=c("analog_climates2100")

for (plot_var in plot_vars){
  jpeg_name=paste(proj_nm0,"_",plot_var, ".jpg", sep = "")
  jpeg(jpeg_name,
       width = 10, height = 10, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  plot(get(plot_var))
  dev.off()
}
cat('\n',sp_nm,'analog climate raster created...')
