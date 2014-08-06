#current_T
#future_T
#D:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg/bioclims_abs/all_baseline/500m
library(raster)
clim_data_dir0=paste0("D:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg/bioclims_abs/all_baseline/500m/")  
clim_data_dir1=paste0("D:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg/bioclims_abs/all_future/500m/")    
files=c("bio1.tif","bio12.tif")
future_img=stack(paste0(clim_data_dir1,files))
baseline_img=stack(paste0(clim_data_dir0,files))
#maskLayer=raster(paste0(clim_data_dir1,files[1]))
#maskLayer=maskLayer<17
maskLayer=raster("Y:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_less_PAs/output_rasters/spp_ensembles_HiRelSpp/current_habitat_spp_em_current_bin.tif")
maskLayer=maskLayer>-1

future_img=mask(future_img,maskLayer, maskvalue=0)
baseline_img=mask(baseline_img,maskLayer, maskvalue=0)


delta_img=future_img-baseline_img
summary_list=cellStats(delta_img, stat='mean', na.rm=TRUE)
#2.6 (all spp, clim only)
#2.55 (hiRelSpp, habitat clip)
#2.5 (all islands)

#plot(delta_img)
drying_areas=delta_img<0
drying_img=(0-delta_img)*drying_areas
wetter_areas=delta_img>0
wetter_img=delta_img*(delta_img>0)
plot(drying_areas)
plot(drying_img)
plot(wetter_areas)
plot(wetter_img)
