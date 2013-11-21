rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))
###USER CONFIGURATION
spp_nm = c('Akikiki','Akekee', 'Hawaii_Amakihi', 'Akiapolauu', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi', 'Kauai_Amakihi','Apapane', 'Iiwi', 'Amakihi', 'Elepaio') 
project_name='finalmodel_P_PA_oldcode'

comp_projects=c('baseline', 'future') #put future second!
ensemble_type="ef.pmw"
eval_stats=c('ROC') 
plot_CV=T
#eval_stats=c('ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS') 
masked=FALSE
veg_overlay=FALSE #this code is not complete!
projected_veg_overlay=FALSE #this code is not complete!
server=0

if (server==1){
  working_dir='Y:/FB_analysis/FB_Base_And_Future/all_DD_merged/'
  clim_data_dir="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/250m/" 
}else{
  working_dir=paste0(resultsDir,project_name,'/')
  clim_data_dir=paste0(bioclimData2013Dir,"all_baseline/500m/")
}
overwrite=1 #if 1, will overwrite past results
current_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/current_veg_mask/"
projected_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/projected_veg_mask/"


####START UNDERHOOD
setwd(working_dir)
library(biomod2)
library(stringr)

save_raster_fx=function(raster_img,out_nm){
  jpeg_name=paste(out_nm, ".jpg", sep = "")
  out_raster_name=paste(out_nm, ".tif", sep = "")
  jpeg(jpeg_name,
       width = 10, height = 8, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  plot(raster_img)
  dev.off()    
  writeRaster(raster_img, out_raster_name, format="GTiff", overwrite=TRUE)        
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


dir.create('output_rasters/main/',showWarnings=F)
sp_nm=spp_nm[1]
eval_stat=eval_stats[1]
for (eval_stat in eval_stats){
  for (sp_nm in spp_nm){
    sp_nm=as.character(sp_nm)  
    cat('\n',sp_nm,'modeling...')
    sp_nm0=sp_nm
    sp_nm=str_replace_all(sp_nm,"_", ".")
    
    out_nm=paste('output_rasters/main/', sp_nm0,"_response_zones_",eval_stat, "_", ensemble_type, sep = "")
    out_raster_name00=paste(out_nm,".tif", sep = "")
    if (file.exists(out_raster_name00)==F | overwrite==1){
      
      ##binary maps
      raster_names=c("EM_suitability1", "EM_suitability2")
      raster_names_bin=c("EM_BIN1", "EM_BIN2")
      i=1
      for (i in c(1,2)){
        proj_nm=comp_projects[i]
        raster_name=raster_names[i]
        raster_name_bin=raster_names_bin[i]
        file_name1=paste(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_",sp_nm,"_TotalConsensus_EMby",eval_stat,".grd", sep = "")
        temp_raster=stack(file_name1)
        band_n=which(names(temp_raster)==paste0(sp_nm,"_TotalConsensus_EMby",eval_stat,"_",ensemble_type))
        assign(raster_name, raster(temp_raster, layer=band_n)/1000)
        
        file_name1_bin=paste(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_",sp_nm,"_TotalConsensus_EMby",eval_stat,"_",eval_stat,"bin.grd", sep = "")
        temp_raster_bin=stack(file_name1_bin)  
        band_n=which(names(temp_raster_bin)==paste0(sp_nm,"_TotalConsensus_EMby",eval_stat,"_",ensemble_type,"_",eval_stat,"bin"))
        assign(raster_name_bin, raster(temp_raster_bin, layer=band_n))
        
        if (plot_CV){
          band_n=which(names(temp_raster)==paste0(sp_nm,"_TotalConsensus_EMby",eval_stat,"_ef.cv"))
          assign(paste0(raster_name,"_CV"), raster(temp_raster, layer=band_n)/1000)
          
          #output suitability rasters for each image
          out_nm=paste('output_rasters/', sp_nm0,"_", "suitability_CV_",proj_nm,"_",eval_stat,"_",ensemble_type, sep = "")
          Process_raster_data_NeutraltoBad(get(paste0(raster_name,"_CV")),out_nm, mask_data=mask_layer)          
        }
        
        #output suitability rasters for each image
        out_nm=paste('output_rasters/', sp_nm0,"_", "suitability_",proj_nm,"_",eval_stat,"_",ensemble_type, sep = "")
        Process_raster_data_NeutraltoGood(get(raster_name),out_nm,min_lim=0, max_lim=1, mask_data=mask_layer)
        
        #output bin rasters for each image    
        out_nm=paste('output_rasters/', sp_nm0,"_", "BIN_",proj_nm,"_",eval_stat,"_",ensemble_type, sep = "")
        save_raster_fx(get(raster_name_bin),out_nm)
        
        
      }
      cat('\n','done with loading baseline and future rasters for ', sp_nm)
      
      #masked suitability
      masked_suitability1=EM_BIN1*EM_suitability1
      out_nm=paste('output_rasters/', sp_nm0,"_", "clipped_suitability_",comp_projects[1],"_",eval_stat,"_",ensemble_type, sep = "")
      #save_raster_fx(masked_suitability1,out_nm)
      Process_raster_data_NeutraltoGood(masked_suitability1,out_nm,min_lim=0, max_lim=1, mask_data=mask_layer)
      
      masked_suitability2=EM_BIN2*EM_suitability2
      out_nm=paste('output_rasters/', sp_nm0,"_", "clipped_suitability_",comp_projects[2],"_",eval_stat,"_",ensemble_type, sep = "")
      #save_raster_fx(masked_suitability2,out_nm)
      Process_raster_data_NeutraltoGood(masked_suitability2,out_nm,min_lim=0, max_lim=1, mask_data=mask_layer)
      suitability_change=EM_suitability2-EM_suitability1
      out_nm=paste('output_rasters/', sp_nm0,"_", "suitability_change_",eval_stat,"_",ensemble_type, sep = "")
      #save_raster_fx(suitability_change,out_nm)
      Process_raster_data_BadtoGood(suitability_change,out_nm,min_lim=-1, max_lim=1, mask_data=mask_layer)
      
      jnk=EM_BIN2*10
      BIN_dif=EM_BIN1+jnk
      m  =  c(9.9,  10.1,  3, 10.9, 11.1, 2)
      rclmat  =  matrix(m,  ncol=3,  byrow=TRUE)
      resp_zone  =  reclassify(BIN_dif,  rclmat)
      
      
      mypalette_numbers=c(0, 1, 2, 3)
      mypalette=c("Grey", "Red", "Green", "Yellow")
      resp_zone_names0=c("Lost", "Overlap", "Gained")
      
      #   jnk=unique(resp_zone)
      #   zones_present=jnk[jnk>0]
      #   zones_present=zones_present[zones_present<=3]
      #   resp_zone_colors=mypalette[zones_present+1]
      #   mypalette_numbers_selected=mypalette[jnk+1] #CHANGED
      #   zones_present=resp_zone_names0[zones_present]
      
      #   jnk=match(mypalette_numbers, jnk, nomatch = 0)
      #   jnk=which(jnk>0)
      #   mypalette1=mypalette[jnk]
      #   jnk=jnk[jnk>0]
      #   jnk=jnk[jnk<=3]
      #   resp_zone_names=resp_zone_names0[jnk]
      
      
      #   mypalette_numbers=c(0, 1, 2, 3)
      #   mypalette=c("Grey", "Red", "Green", "Yellow")
      #   resp_zone_names0=c("Micro refugia", "Tolerate", "Migrate")  
      
      
      
      if (masked){
        current_mask=EM_suitability1>minValue(EM_suitability1)  
        raster_name=paste(sp_nm0,"_analog_climates2100", sep="")
        tif_name=paste(raster_name,".tif", sep="")
        analog_cc = raster( tif_name)
        habitat = raster( paste(clim_data_dir, "veg_areas.grd", sep=""))
        habitat=crop(habitat, analog_cc)
        
        all_mask=habitat+analog_cc*2+current_mask*4 #1 cur, 2 ang, 4 hab, 3 cur/ang, 6 ang/hab, 7 cur/ang/hab
        cum_mask=current_mask*analog_cc*habitat
        
        masked_resp_zone=resp_zone*cum_mask
        cat('\n','created mask for ', sp_nm)
        
        jpeg_name=paste('output_rasters/main/', sp_nm0, "_mask.jpg", sep = "")
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        plot(all_mask)
        dev.off()
        
        ##MASKED bin comparison rasters
        out_nm=paste('output_rasters/main/', sp_nm0,"_response_zones_masked_",eval_stat, "_", ensemble_type, sep = "")
        jpeg_name=paste(out_nm, ".jpg", sep = "")
        out_raster_name=paste(out_nm,".tif", sep = "")
        
        jnk=unique(masked_resp_zone)
        graph_palette=mypalette_numbers
        zones_present=jnk[jnk>0]
        zones_present=zones_present[zones_present<=3]
        resp_zone_colors=mypalette[zones_present+1]
        resp_zone_names=resp_zone_names0[zones_present]
        mypalette_numbers_selected=mypalette[jnk+1] #CHANGED
        
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        plot(masked_resp_zone,  col=mypalette_numbers_selected, legend=F)
        #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
        legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
        dev.off()  
        writeRaster(masked_resp_zone, out_raster_name, format="GTiff", overwrite=TRUE)
        
        future_bin_with_mask=cum_mask*EM_BIN2
        future_suitability_with_mask=cum_mask*EM_suitability2
        
        #output suitability rasters for each image
        #suitability
        out_nm=paste('output_rasters/', sp_nm0,"_suitability_future_masked_",eval_stat,"_",ensemble_type, sep = "")
        jpeg_name=paste(out_nm,".jpg", sep = "")
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        plot(future_suitability_with_mask)
        dev.off()
        
        out_raster_name=paste(out_nm,".tif", sep = "")
        writeRaster(future_suitability_with_mask, out_raster_name, format="GTiff", overwrite=TRUE)
        
        ##binary
        out_nm=paste('output_rasters/', sp_nm0,"_BIN_future_masked_",eval_stat,"_",ensemble_type, sep = "")
        jpeg_name=paste(out_nm, ".jpg", sep = "")
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        plot(future_bin_with_mask)
        dev.off()
        
        out_raster_name=paste(out_nm,".tif", sep = "")
        writeRaster(future_bin_with_mask, out_raster_name, format="GTiff", overwrite=TRUE)
      }
      
      if (veg_overlay){
        zone_raster_vals0=unique(resp_zone)
        #legend
        zones_present=zone_raster_vals0[zone_raster_vals0>0]
        resp_zone_names=resp_zone_names0[zones_present]
        jnk0=match(zones_present,mypalette_numbers, nomatch = 0)
        resp_zone_colors=mypalette[jnk0]
        
        #with current veg overlay
        Response_var=paste(current_biome_distribution_dir, sp_nm0,"_current_veg_mask.tif", sep="")
        Response_var=raster(Response_var)
        temp_resp_zone=(Response_var*4)+resp_zone
        
        #three lines below are to fix a bug in the plot (it was not considering raster values that had very low freq, hence incorrectly applying the color pallete)
        jnk=freq(temp_resp_zone)
        todel=jnk[jnk[,'count']<=1,'value']
        temp_resp_zone[temp_resp_zone %in% todel]=0
        
        zone_raster_vals=unique(temp_resp_zone)    
        #raster color
        jnk0=match(zone_raster_vals,mypalette_numbers, nomatch = 0)
        mypalette1=mypalette[jnk0]
        
        jpeg_name=paste('jpg_outputs/', community,"_response_zones_w_current_", proj_name,".jpg", sep = "")
        
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        plot(temp_resp_zone,  col=mypalette1, legend=F)
        #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
        legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
        title(paste(community_nm, " climate envelope shift", sep=""))
        dev.off()
        
        out_raster_name0=paste(working_dir, "tifs/", community, "_resp_envelopes_with_current_veg", proj_name,".tif", sep="")
        writeRaster(temp_resp_zone, out_raster_name0, format="GTiff", overwrite=TRUE) 
        
      }
      
      ##bin comparison rasters
      out_nm=paste('output_rasters/main/', sp_nm0,"_response_zones_",eval_stat, "_", ensemble_type, sep = "")
      jpeg_name=paste(out_nm, ".jpg", sep = "")
      out_raster_name00=paste(out_nm,".tif", sep = "")
      
      jnk=unique(resp_zone)
      graph_palette=mypalette_numbers
      zones_present=jnk[jnk>0]
      zones_present=zones_present[zones_present<=3]
      resp_zone_colors=mypalette[zones_present+1]
      resp_zone_names=resp_zone_names0[zones_present]
      mypalette_numbers_selected=mypalette[jnk+1] #CHANGED
      
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(resp_zone,  col=mypalette_numbers_selected, legend=F)
      #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
      legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
      dev.off()
      writeRaster(resp_zone, out_raster_name00, format="GTiff", overwrite=TRUE)
    }else{
      cat('\n', sp_nm, " already calculated")
    }
  }
}


