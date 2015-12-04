###USER CONFIGURATION
veg_overlay=T #should always be on
projected_veg_overlay=F #this code is not complete!

if(BPS){
  current_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat_analysis/veg_overlay/BPS/current_veg_mask/" 
  projected_veg_overlay=F
  BPS_str='BPS/'
}else{
  current_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat_analysis/veg_overlay/current_veg_mask/"
  BPS_str=''
}
projected_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat_analysis/veg_overlay/projected_veg_mask/"


####START UNDERHOOD
library(biomod2)
library(stringr)

dir.create('output_rasters/veg_overlays/',showWarnings=FALSE)
if(BPS){
  dir.create('output_rasters/veg_overlays/BPS/',showWarnings=FALSE)  
}
#sp_nm=spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  sp_nm=str_replace_all(sp_nm,"_", ".")
  
  resp_zone_file_name=paste('output_rasters/main/', sp_nm0,"_response_zones_",veg_overlay_eval_stat, "_", spp_ensemble_type, "_", comp_projects[2], ".tif", sep = "")
  resp_zone=raster(resp_zone_file_name)  
  if (veg_overlay){
    jpeg_name=paste('output_rasters/veg_overlays/', BPS_str, sp_nm0,"_response_zones_w_current_veg_distribution.jpg", sep = "")
    if (file.exists(jpeg_name)==F | overwrite==1){ #check to see if the analysis for this species was already done    
      
      mypalette_numbers=c(0, 1, 2, 3, 4, 5, 6, 7)
      mypalette_breaks=c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5)
      mypalette=c("Grey", "lightcoral", "lightgreen", "gold", "darkgrey", "darkred", "darkgreen", "darkorange")
      resp_zone_names0=c("Lost", "Overlap", "Gained")
      
      zone_raster_vals0=unique(resp_zone)
      #legend
      zones_present=zone_raster_vals0[zone_raster_vals0>0]
      resp_zone_names=resp_zone_names0[zones_present]
      jnk0=match(zones_present,mypalette_numbers, nomatch = 0)
      resp_zone_colors=mypalette[jnk0]
      
      #with current veg overlay
      Response_var=paste(current_biome_distribution_dir, sp_nm0,"_current_veg_mask.tif", sep="")
      Response_var=raster(Response_var)
      res_ratio=round(res(resp_zone)[1]/res(Response_var)[1])
      if (res_ratio>1){
        Response_var=aggregate(Response_var,  fact=res_ratio,  fun=max)          
      }
      Response_var=crop(Response_var,resp_zone)
      Response_var1=resample(Response_var,resp_zone,  method="ngb")
      #Response_var1=alignExtent(resp_zone,Response_var)
      temp_resp_zone=(Response_var1*4)+resp_zone
      temp_resp_zone_3Class=(Response_var1*resp_zone)+4
      
      #three lines below are to fix a bug in the plot (it was not considering raster values that had very low freq, hence incorrectly applying the color pallete)
      jnk=freq(temp_resp_zone)
      todel=jnk[jnk[,'count']<=1,'value']
      temp_resp_zone[temp_resp_zone %in% todel]=0
      
      zone_raster_vals=unique(temp_resp_zone)    
      
      jpeg_name=paste('output_rasters/veg_overlays/', BPS_str, sp_nm0,"_response_zones_w_current_veg_distribution.jpg", sep = "")
      
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(temp_resp_zone,  breaks=mypalette_breaks, col=mypalette, legend=F)
      #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
      legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
      title(paste(sp_nm0," distribution shifts with current biome distribution", sep=""))
      dev.off()
      
      #####Response zones with clipped habitat (no light hues)
      #three lines below are to fix a bug in the plot (it was not considering raster values that had very low freq, hence incorrectly applying the color pallete)
      jnk=freq(temp_resp_zone_3Class)
      todel=jnk[jnk[,'count']<=1,'value']
      temp_resp_zone_3Class[temp_resp_zone_3Class %in% todel]=0
      
      zone_raster_vals=unique(temp_resp_zone_3Class)    
      
      jpeg_name=paste('output_rasters/veg_overlays/', BPS_str, sp_nm0,"_response_zones_w_current_veg_distribution_3Class.jpg", sep = "")
      
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(temp_resp_zone_3Class,  breaks=mypalette_breaks, col=mypalette, legend=F)
      #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
      legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
      title(paste(sp_nm0," distribution shifts with current biome distribution", sep=""))
      dev.off()
      
      resp_zone_in_habitat=Response_var1*resp_zone
      out_raster_name0=paste('output_rasters/veg_overlays/', BPS_str, sp_nm0,"_response_zones_w_current_habitat.tif", sep = "")
      writeRaster(resp_zone_in_habitat, out_raster_name0, format="GTiff", overwrite=TRUE)     
    }  
  }
  
  if (projected_veg_overlay){
    jpeg_name=paste('output_rasters/veg_overlays/', BPS_str, sp_nm0,"_response_zones_w_future_clim_env_of_primary_habitat.jpg", sep = "")
    if (file.exists(jpeg_name)==F | overwrite==1){ #check to see if the analysis for this species was already done          
      mypalette_numbers=c(0, 1, 2, 3, 4, 5, 6, 7)
      mypalette_breaks=c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5)
      mypalette=c("Grey", "lightcoral", "lightgreen", "gold", "darkgrey", "darkred", "darkgreen", "darkorange")
      resp_zone_names0=c("Lost", "Overlap", "Gained")
      
      zone_raster_vals0=unique(resp_zone)
      #legend
      zones_present=zone_raster_vals0[zone_raster_vals0>0]
      resp_zone_names=resp_zone_names0[zones_present]
      jnk0=match(zones_present,mypalette_numbers, nomatch = 0)
      resp_zone_colors=mypalette[jnk0]
      
      #with projected veg overlay
      Response_var=paste(projected_biome_distribution_dir, sp_nm0,"_projected_veg_mask.tif", sep="")
      Response_var=raster(Response_var)
      Response_var=Response_var>0 #future envelope
      res_ratio=round(res(resp_zone)[1]/res(Response_var)[1])
      if (res_ratio>1){
        Response_var=aggregate(Response_var,  fact=res_ratio,  fun=max)          
      }
      Response_var=crop(Response_var,resp_zone)
      Response_var1=resample(Response_var,resp_zone,  method="ngb")
      #Response_var1=alignExtent(resp_zone,Response_var)
      temp_resp_zone=(Response_var1*4)+resp_zone
      
      #three lines below are to fix a bug in the plot (it was not considering raster values that had very low freq, hence incorrectly applying the color pallete)
      jnk=freq(temp_resp_zone)
      todel=jnk[jnk[,'count']<=1,'value']
      temp_resp_zone[temp_resp_zone %in% todel]=0
      
      zone_raster_vals=unique(temp_resp_zone)    
      #raster color
      #jnk0=match(zone_raster_vals,mypalette_numbers, nomatch = 0)
      #mypalette1=mypalette[jnk0]
      
      jpeg_name=paste('output_rasters/veg_overlays/', BPS_str, sp_nm0,"_response_zones_w_future_clim_env_of_primary_habitat.jpg", sep = "")
      
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(temp_resp_zone,  breaks=mypalette_breaks, col=mypalette, legend=F)
      #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
      legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
      title(paste(sp_nm0," distribution shifts with future climate envelope of primary habitat", sep=""))
      dev.off()
      
      resp_zone_in_habitat=Response_var1*resp_zone
      out_raster_name0=paste('output_rasters/veg_overlays/', BPS_str, sp_nm0,"_response_zones_w_future_clim_env_of_primary_habitat.tif", sep = "")
      writeRaster(resp_zone_in_habitat, out_raster_name0, format="GTiff", overwrite=TRUE)     
    }  
  }
  
  
}


