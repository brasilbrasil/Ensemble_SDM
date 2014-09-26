clim_data_dir = clim_data_2000 

#overwrite=0 #if 1, will overwrite past results
if(BPS){
  current_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/BPS/current_veg_mask/" 
}else{
  current_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/current_veg_mask/"
}

####START UNDERHOOD
#setwd(working_dir)
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

Process_raster_data_NeutraltoGood_W_overlay=function(raster_var,out_nm,min_lim=NULL, max_lim=NULL, mask_data=NULL, overlay_data=NULL){
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
  if (!is.null(overlay_data)){
    plot(overlay_data,add=T, border="red", lwd=3)    
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


sp_nm=spp_nm[1]
eval_stat = spp_ensemble_eval_stats[1]
for (eval_stat in spp_ensemble_eval_stats){
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
  spp_em_current_suitability_bin=jnk
  spp_em_future_suitability=jnk
  spp_em_future_suitability_bin=jnk
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
    
    response_raster=raster(paste('output_rasters/main/', sp_nm0,"_response_zones_",eval_stat, "_", spp_ensemble_type, ".tif", sep = ""))  
    baseline_masked_suitability=raster(paste('output_rasters/', sp_nm0,"_", "clipped_suitability_","baseline","_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
    future_masked_suitability=raster(paste('output_rasters/', sp_nm0,"_", "clipped_suitability_","future","_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
    baseline_suitability=raster(paste('output_rasters/', sp_nm0,"_", "suitability_","baseline","_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
    baseline_suitability_bin=baseline_suitability>0
    future_suitability=raster(paste('output_rasters/', sp_nm0,"_", "suitability_","future","_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
    future_suitability_bin=future_suitability>0
    baseline_suitability_CV=raster(paste0('output_rasters/', sp_nm0,"_", "suitability_CV_","baseline","_",eval_stat,"_",spp_ensemble_type, ".tif"))
    future_suitability_CV=raster(paste0('output_rasters/', sp_nm0,"_", "suitability_CV_","future","_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
    current_bin=raster(paste('output_rasters/', sp_nm0,"_", "BIN_","baseline","_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
    future_bin=raster(paste('output_rasters/', sp_nm0,"_", "BIN_","future","_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
    suitability_delta=raster(paste('output_rasters/', sp_nm0,"_", "suitability_change_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
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
      baseline_suitability_bin=baseline_suitability_bin*Response_var1
      future_suitability=future_suitability*Response_var1
      future_suitability_bin=future_suitability_bin*Response_var1
      baseline_suitability_CV=baseline_suitability_CV*Response_var1
      future_suitability_CV=future_suitability_CV*Response_var1
      current_bin=current_bin*Response_var1
      future_bin=future_bin*Response_var1
      suitability_delta=suitability_delta*Response_var1
      lost_range=lost_range*Response_var1
      kept_range=kept_range*Response_var1
      gained_range=gained_range*Response_var1    
      if (BPS){
        masked_text="BPS_habitat_"
      }else{
        masked_text="current_habitat_"      
      }
      
    }else{
      masked_text=""
    }
    
    
    ##align extent/ mask all water
    baseline_masked_suitability=extend_raster(baseline_masked_suitability)
    future_masked_suitability=extend_raster(future_masked_suitability)  
    baseline_suitability=extend_raster(baseline_suitability)
    baseline_suitability_bin=extend_raster(baseline_suitability_bin)
    future_suitability=extend_raster(future_suitability)
    future_suitability_bin=extend_raster(future_suitability_bin)
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
    spp_em_current_suitability_bin=spp_em_current_suitability_bin+baseline_suitability_bin
    spp_em_future_suitability=spp_em_future_suitability+future_suitability
    spp_em_future_suitability_bin=spp_em_future_suitability_bin+future_suitability_bin
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
  
  Process_raster_data_NeutraltoGood(spp_em_current_suitability/spp_em_current_suitability_bin, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_current_suitability_',eval_stat), max_lim=1, min_lim=0, mask_data=mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_future_suitability/spp_em_future_suitability_bin, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_future_suitability_',eval_stat), max_lim=1, min_lim=0, mask_data=mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_masked_current_suitability/spp_em_current_bin, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_masked_current_suitability_',eval_stat), max_lim=1, min_lim=0, mask_data=mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_masked_future_suitability/spp_em_future_bin, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_masked_future_suitability_',eval_stat),max_lim=1, min_lim=0,  mask_data=mask_layer)
  Process_raster_data_NeutraltoBad(spp_em_masked_current_suitability_CV, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_current_suitability_CV_',eval_stat), mask_data=mask_layer)
  Process_raster_data_NeutraltoBad(spp_em_masked_future_suitability_CV, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_future_suitability_CV_',eval_stat), mask_data=mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_current_bin, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_current_bin_',eval_stat), mask_data=mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_future_bin, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_future_bin_',eval_stat), mask_data=mask_layer)
  avg_spp_em_suitability_delta=(spp_em_future_suitability/spp_em_future_suitability_bin)-(spp_em_current_suitability/spp_em_current_suitability_bin)
  max_val=max(abs(c(cellStats(avg_spp_em_suitability_delta,max),cellStats(avg_spp_em_suitability_delta,min))))
  Process_raster_data_BadtoGood(avg_spp_em_suitability_delta, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_suitability_avg_delta_',eval_stat), max_lim=max_val, min_lim=-max_val, mask_data=mask_layer)
  
  max_val=max(abs(c(cellStats(spp_em_suitability_delta,max),cellStats(spp_em_suitability_delta,min))))
  Process_raster_data_BadtoGood(spp_em_suitability_delta, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_suitability_delta_',eval_stat), max_lim=max_val, min_lim=-max_val, mask_data=mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_lost_range, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_lost_range_',eval_stat), mask_data=mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_kept_range, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_kept_range_',eval_stat), mask_data=mask_layer)
  Process_raster_data_NeutraltoGood(spp_em_gained_range, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_gained_range_',eval_stat), mask_data=mask_layer)
  
  #MS figures
  prot_areas=shapefile(paste0("Y:/PICCC_analysis/FB_analysis/habitat_analysis/","protected_areas_20100331_simpleWGS1984.shp"))
  Process_raster_data_NeutraltoGood_W_overlay(spp_em_kept_range, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_kept_ranges_w_prot_areas_',eval_stat), mask_data=mask_layer, overlay_data=prot_areas)
  
  if (BPS){
    #future restoration priority
    FUT_w_BPS_hab=raster(paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_future_bin.tif'))  
    FUT_w_Curr_hab=raster(paste0('output_rasters/spp_ensembles/',"current_habitat_",'spp_em_future_bin.tif'))
    not_cur_hab=FUT_w_Curr_hab==0
    
    restoration_priority=FUT_w_BPS_hab*not_cur_hab
    Process_raster_data_NeutraltoGood(restoration_priority, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_restoration_priority_future_',eval_stat), mask_data=mask_layer)  
    
    #current restoration priority
    CUR_w_BPS_hab=raster(paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_current_bin.tif'))  
    CUR_w_Curr_hab=raster(paste0('output_rasters/spp_ensembles/',"current_habitat_",'spp_em_current_bin.tif'))
    not_cur_hab=CUR_w_Curr_hab==0
    
    restoration_priority=CUR_w_BPS_hab*not_cur_hab
    Process_raster_data_NeutraltoGood(restoration_priority, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_restoration_priority_current_',eval_stat), mask_data=mask_layer)  
    
    #gained restoration priority
    CUR_w_BPS_hab=raster(paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_gained_range.tif'))  
    CUR_w_Curr_hab=raster(paste0('output_rasters/spp_ensembles/',"current_habitat_",'spp_em_gained_range.tif'))
    not_cur_hab=CUR_w_Curr_hab==0
    
    restoration_priority=CUR_w_BPS_hab*not_cur_hab
    Process_raster_data_NeutraltoGood(restoration_priority, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_restoration_priority_gained_',eval_stat), mask_data=mask_layer)  
    
    #kept restoration priority
    CUR_w_BPS_hab=raster(paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_kept_range.tif'))  
    CUR_w_Curr_hab=raster(paste0('output_rasters/spp_ensembles/',"current_habitat_",'spp_em_kept_range.tif'))
    not_cur_hab=CUR_w_Curr_hab==0
    
    restoration_priority=CUR_w_BPS_hab*not_cur_hab
    Process_raster_data_NeutraltoGood(restoration_priority, paste0('output_rasters/spp_ensembles/',masked_text,'spp_em_restoration_priority_kept_',eval_stat), mask_data=mask_layer)  
    
    #MERGE current and future restoration priorities
    #NA/ 0 as no display    
  }
}
