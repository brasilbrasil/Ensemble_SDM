rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))
###USER CONFIGURATION
spp_nm = c('Akekee', 'Hawaii_Amakihi', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi','Kauai_Amakihi', 'Iiwi', 'Apapane')#'Amakihi', 'Elepaio', 
project_name='finalmodel_P_PA_oldcode_less_PAs'
ensemble_type="ef.pmw"
eval_stats=c('ROC') 
habitat_overlay=F

working_dir=paste0(resultsDir,project_name,'/')
clim_data_dir=paste0(bioclimData2013Dir,"all_baseline/500m/")
overwrite=0 #if 1, will overwrite past results
current_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/current_veg_mask/"
projected_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/projected_veg_mask/"

####START UNDERHOOD
setwd(working_dir)
library(biomod2)
library(stringr)
dir.create("tables/", showWarnings=F)
vars=c("Eval stat", "Species", "baseline_area", "future_area", "% change", "lost", "kept", "gained", "baseline_suitability", "future_suitability")
island_mask=raster(paste0(clim_data_dir,"bio1.tif"))
island_mask[!is.na(island_mask)]=1

extend_raster=function(raster_lyr){
  raster_lyr[raster_lyr==NA]=0
  raster_lyr=extend(raster_lyr,island_mask,value=0)
  expanded_raster_lyr=raster_lyr*island_mask
  return(expanded_raster_lyr)  
}

sp_nm=spp_nm[1]
eval_stat = eval_stats[1]

#current suitability
#future suitability
#masked current suitability
#masked future suitability
#current bin
#future bin
#suitability_delta
#lost_range
#kept_range
#gained_range
jnk=island_mask
jnk[jnk==1]=0
spp_em_masked_current_suitability=jnk
spp_em_masked_future_suitability=jnk
spp_em_masked_current_suitability_CV=jnk
spp_em_masked_future_suitability_CV=jnk
spp_em_current_suitability=jnk
spp_em_future_suitability=jnk
spp_em_current_bin=jnk
spp_em_future_bin=jnk
spp_em_suitability_delta=jnk
spp_em_lost_range=jnk
spp_em_kept_range=jnk
spp_em_gained_range=jnk

sp_nm = spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  sp_nm=str_replace_all(sp_nm,"_", ".")
  
  response_raster=raster(paste('output_rasters/main/', sp_nm0,"_response_zones_",eval_stat, "_", ensemble_type, ".tif", sep = ""))  
  baseline_masked_suitability=raster(paste('output_rasters/', sp_nm0,"_", "clipped_suitability_","baseline","_",eval_stat,"_",ensemble_type, ".tif", sep = ""))
  future_masked_suitability=raster(paste('output_rasters/', sp_nm0,"_", "clipped_suitability_","future","_",eval_stat,"_",ensemble_type, ".tif", sep = ""))
  baseline_suitability=raster(paste('output_rasters/', sp_nm0,"_", "suitability_","baseline","_",eval_stat,"_",ensemble_type, ".tif", sep = ""))
  future_suitability=raster(paste('output_rasters/', sp_nm0,"_", "suitability_","future","_",eval_stat,"_",ensemble_type, ".tif", sep = ""))
  baseline_suitability_CV=raster(paste0('output_rasters/', sp_nm0,"_", "suitability_CV_","baseline","_",eval_stat,"_",ensemble_type, ".tif"))
  future_suitability_CV=raster(paste0('output_rasters/', sp_nm0,"_", "suitability_CV_","future","_",eval_stat,"_",ensemble_type, ".tif", sep = ""))
  current_bin=raster(paste('output_rasters/', sp_nm0,"_", "BIN_","baseline","_",eval_stat,"_",ensemble_type, ".tif", sep = ""))
  future_bin=raster(paste('output_rasters/', sp_nm0,"_", "BIN_","future","_",eval_stat,"_",ensemble_type, ".tif", sep = ""))
  suitability_delta=raster(paste('output_rasters/', sp_nm0,"_", "suitability_change_",eval_stat,"_",ensemble_type, ".tif", sep = ""))
  lost_range=response_raster==1
  kept_range=response_raster==2
  gained_range=response_raster==3

  if (habitat_overlay){
    #with current veg overlay
    Response_var=paste(current_biome_distribution_dir, sp_nm0,"_current_veg_mask.tif", sep="")
    Response_var=raster(Response_var)
    res_ratio=round(res(response_raster)[1]/res(Response_var)[1])
    if (res_ratio>1){
      Response_var=aggregate(Response_var,  fact=res_ratio,  fun=max)          
    }
    Response_var=crop(Response_var,response_raster)
    Response_var1=resample(Response_var,response_raster,  method="ngb")

    response_raster=response_raster*Response_var1  
    baseline_masked_suitability=baseline_masked_suitability*Response_var1
    future_masked_suitability=future_masked_suitability*Response_var1
    baseline_suitability=baseline_suitability*Response_var1
    future_suitability=future_suitability*Response_var1
    baseline_suitability_CV=baseline_suitability_CV*Response_var1
    future_suitability_CV=future_suitability_CV*Response_var1
    current_bin=current_bin*Response_var1
    future_bin=future_bin*Response_var1
    suitability_delta=suitability_delta*Response_var1
    lost_range=lost_range*Response_var1
    kept_range=kept_range*Response_var1
    gained_range=gained_range*Response_var1    
    masked_text="current_habitat_"
  }else{
    masked_text=""
  }
  
  
  ##align extent/ mask all water
  baseline_masked_suitability=extend_raster(baseline_masked_suitability)
  future_masked_suitability=extend_raster(future_masked_suitability)  
  baseline_suitability=extend_raster(baseline_suitability)
  future_suitability=extend_raster(future_suitability)
  baseline_suitability_CV=extend_raster(baseline_suitability_CV)
  future_suitability_CV=extend_raster(future_suitability_CV)
  current_bin=extend_raster(current_bin)
  future_bin=extend_raster(future_bin)
  suitability_delta=extend_raster(suitability_delta)
  lost_range=extend_raster(lost_range)
  kept_range=extend_raster(kept_range)
  gained_range=extend_raster(gained_range)
    
  ##calculations
  spp_em_masked_current_suitability=spp_em_masked_current_suitability+baseline_masked_suitability
  spp_em_masked_future_suitability=spp_em_masked_future_suitability+future_masked_suitability
  spp_em_current_suitability=spp_em_current_suitability+baseline_suitability
  spp_em_future_suitability=spp_em_future_suitability+future_suitability
  spp_em_masked_current_suitability_CV=spp_em_masked_current_suitability_CV+baseline_suitability_CV
  spp_em_masked_future_suitability_CV=spp_em_masked_future_suitability_CV+future_suitability_CV
  spp_em_current_bin=spp_em_current_bin+current_bin
  spp_em_future_bin=spp_em_future_bin+future_bin
  spp_em_suitability_delta=spp_em_suitability_delta+suitability_delta
  spp_em_lost_range=spp_em_lost_range+lost_range
  spp_em_kept_range=spp_em_kept_range+kept_range
  spp_em_gained_range=spp_em_gained_range+gained_range
    
}

##save rasters
dir.create('output_rasters/spp_ensembles/',showWarnings=F)
Process_raster_data_BadtoGood=function(raster_var,out_nm,min_lim=NULL, max_lim=NULL, mask_data=NULL){
  jpeg_name=paste(out_nm, ".jpg", sep = "")
  out_raster_name=paste(out_nm, ".tif", sep = "")
  jpeg(jpeg_name,
       width = 14, height = 14, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  if (is.null(min_lim)){
    min_lim=minValue(raster_var)
  }
  if (is.null(max_lim)){
    max_lim=maxValue(raster_var)
  }
  library(colorRamps)
  col5 <- colorRampPalette(c('red', 'gray96', 'darkgreen'))  
  plot(raster_var, col=col5(n=99), breaks=seq(min_lim,max_lim,length.out=100) , axes=FALSE, box=FALSE,
       legend=T, legend.width=1, legend.shrink=0.75,
       legend.args=list(text="", side=4, font=2, line=2.5, cex=0.8),
       axis.args=list(at=seq(min_lim,max_lim, (max_lim-min_lim)/10),
                      labels=seq(min_lim,max_lim, (max_lim-min_lim)/10)))
  if (!is.null(mask_data)){
    plot(mask_data,add=T)    
  }
  dev.off()  
  
  writeRaster(raster_var, out_raster_name, format="GTiff", overwrite=TRUE)
}

Process_raster_data_NeutraltoGood=function(raster_var,out_nm,min_lim=NULL, max_lim=NULL, mask_data=NULL){
  jpeg_name=paste(out_nm, ".jpg", sep = "")
  out_raster_name=paste(out_nm, ".tif", sep = "")
  jpeg(jpeg_name,
       width = 14, height = 14, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  if (is.null(min_lim)){
    min_lim=minValue(raster_var)
  }
  if (is.null(max_lim)){
    max_lim=maxValue(raster_var)
  }
  library(colorRamps)
  col5 <- colorRampPalette(c('gray96', 'darkgreen'))#lightgrey  
  plot(raster_var, col=col5(n=99), breaks=seq(min_lim,max_lim,length.out=100) , axes=FALSE, box=FALSE,
       legend=T, legend.width=1, legend.shrink=0.75,
       legend.args=list(text="", side=4, font=2, line=2.5, cex=0.8),
       axis.args=list(at=seq(min_lim,max_lim, (max_lim-min_lim)/10),
                      labels=seq(min_lim,max_lim, (max_lim-min_lim)/10)))
  if (!is.null(mask_data)){
    plot(mask_data,add=T)    
  }
  dev.off()  
  
  writeRaster(raster_var, out_raster_name, format="GTiff", overwrite=TRUE)
}

Process_raster_data_NeutraltoBad=function(raster_var,out_nm,min_lim=NULL, max_lim=NULL, mask_data=NULL){
  jpeg_name=paste(out_nm, ".jpg", sep = "")
  out_raster_name=paste(out_nm, ".tif", sep = "")
  jpeg(jpeg_name,
       width = 14, height = 14, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  if (is.null(min_lim)){
    min_lim=minValue(raster_var)
  }
  if (is.null(max_lim)){
    max_lim=maxValue(raster_var)
  }
  library(colorRamps)
  col5 <- colorRampPalette(c('gray96', 'red'))#lightgrey  
  plot(raster_var, col=col5(n=99), breaks=seq(min_lim,max_lim,length.out=100) , axes=FALSE, box=FALSE,
       legend=T, legend.width=1, legend.shrink=0.75,
       legend.args=list(text="", side=4, font=2, line=2.5, cex=0.8),
       axis.args=list(at=seq(min_lim,max_lim, (max_lim-min_lim)/10),
                      labels=seq(min_lim,max_lim, (max_lim-min_lim)/10)))
  if (!is.null(mask_data)){
    plot(mask_data,add=T)    
  }
  dev.off()  
  
  writeRaster(raster_var, out_raster_name, format="GTiff", overwrite=TRUE)
}

shapedir=paste0(DR_PICCC_data_S,"/climate_data/bioclim_data_Aug2013/original_raw_data/HRMC20130810_20yrs_of_3km_16yrs_of_1km_maui/")
mask_layer=shapefile(paste0(shapedir,"Main_Hawaiian_Islands_simple3.shp"))


Process_raster_data_NeutraltoGood(spp_em_current_suitability, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_current_suitability'), mask_data=mask_layer)
Process_raster_data_NeutraltoGood(spp_em_future_suitability, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_future_suitability'), mask_data=mask_layer)
Process_raster_data_NeutraltoGood(spp_em_masked_current_suitability, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_masked_current_suitability'), mask_data=mask_layer)
Process_raster_data_NeutraltoGood(spp_em_masked_future_suitability, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_masked_future_suitability'), mask_data=mask_layer)
Process_raster_data_NeutraltoBad(spp_em_masked_current_suitability_CV, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_current_suitability_CV'), mask_data=mask_layer)
Process_raster_data_NeutraltoBad(spp_em_masked_future_suitability_CV, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_future_suitability_CV'), mask_data=mask_layer)
Process_raster_data_NeutraltoGood(spp_em_current_bin, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_current_bin'), mask_data=mask_layer)
Process_raster_data_NeutraltoGood(spp_em_future_bin, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_future_bin'), mask_data=mask_layer)
max_val=max(abs(c(cellStats(spp_em_suitability_delta,max),cellStats(spp_em_suitability_delta,min))))
Process_raster_data_BadtoGood(spp_em_suitability_delta, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_suitability_delta'), max_lim=max_val, min_lim=-max_val, mask_data=mask_layer)
Process_raster_data_NeutraltoGood(spp_em_lost_range, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_lost_range'), mask_data=mask_layer)
Process_raster_data_NeutraltoGood(spp_em_kept_range, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_kept_range'), mask_data=mask_layer)
Process_raster_data_NeutraltoGood(spp_em_gained_range, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_gained_range'), mask_data=mask_layer)
