rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))
###USER CONFIGURATION
#source(paste0("Y:/PICCC_analysis/code/","directory_registry.r"))
#local_config_dir=paste0(DR_FB_SDM_results_S,'test_runs_500m/') #'C:/Users/lfortini/'
#spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Anianiau", "Kauai_Amakihi", "Hawaii_Elepaio", "Palila")
spp_nm=c('Akekee', 'Apapane', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Hawaii_Akepa', 'Anianiau', 'Hawaii_Creeper', 'Kauai_Amakihi', 'Hawaii_Elepaio', 'Kauai_Elepaio', 'Amakihi', 'Maui_Alauahio', 'Hawaii_Amakihi', 'Maui_Parrotbill', 'Oahu_Elepaio', 'Palila','Elepaio')
#'Iiwi',  'Puaiohi',  'Oahu_Amakihi',  'Omao', 
project_name='finalmodel_P_A_PA_oldcode'



comp_projects=c('baseline', 'future') #put future second!
ensemble_type="ef.pmw"
#working_dir='Y:/FB analysis/FB SDM/biomod2/'
eval_stat="ROC"
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
overwrite=0 #if 1, will overwrite past results
current_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/current_veg_mask/"
projected_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/projected_veg_mask/"


####START UNDERHOOD
setwd(working_dir)
library(biomod2)
library(stringr)

sp_nm=spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  sp_nm=str_replace_all(sp_nm,"_", ".")
  a=round(1000*runif(1))
  
  out_nm=paste('output_rasters/main/', sp_nm0,"_response_zones_",eval_stat, "_", ensemble_type, sep = "")
  out_raster_name00=paste(out_nm,".tif", sep = "")
  if (file.exists(out_raster_name00)==F | overwrite==1){
    
    raster_names=c("EM_suitability1", "EM_suitability2")
    raster_names_bin=c("EM_BIN1", "EM_BIN2")
    i=1
    for (i in c(1,2)){
      proj_nm=comp_projects[i]
      raster_name=raster_names[i]
      raster_name_bin=raster_names_bin[i]
      file_name1=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat, sep = "")
      file_name2=paste(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_",sp_nm,"_TotalConsensus_EMby",eval_stat,".grd", sep = "")
      if (file.exists(file_name1)){ 
        load(file_name1)
        assign("temp_raster", get(paste(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat, sep = "")))  
        band_names=names(temp_raster)
        band_n=which(band_names==ensemble_type)
        assign(raster_name, raster(temp_raster, layer=band_n))
      }else{
        temp_raster=raster(file_name2)
        assign(raster_name, temp_raster)
        #assign("temp_raster", get(paste(sp_nm,"_AllData_AllRun_EM.",eval_stat, sep = "")))      
      }
      
      file_name1_bin=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat,".bin.",eval_stat, sep = "")
      file_name2_bin=paste(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_",sp_nm,"_TotalConsensus_EMby",eval_stat,"_",eval_stat,"bin.grd", sep = "")
      if (file.exists(file_name1_bin)){       
        load(file_name1_bin)
        assign("temp_raster_bin", get(paste(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat,".bin.",eval_stat, sep = "")))      
        band_names=names(temp_raster_bin)
        band_n=which(band_names==paste(ensemble_type,".bin", sep = ""))
        assign(raster_name_bin, raster(temp_raster_bin, layer=band_n))    
        }else{
        temp_raster_bin=raster(file_name2_bin)  
        assign(raster_name_bin, temp_raster_bin)    
        
      }
            
      #output suitability rasters for each image
      out_nm=paste('output_rasters/', sp_nm0,"_", "suitability_",proj_nm,"_",eval_stat,"_",ensemble_type,".jpg", sep = "")
      jpeg_name=paste(out_nm, ".jpg", sep = "")
      out_raster_name=paste(out_nm, ".jpg", sep = "")
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(get(raster_name))
      dev.off()    
      writeRaster(get(raster_name), out_raster_name, format="GTiff", overwrite=TRUE)
      
      #output bin rasters for each image    
      out_nm=paste('output_rasters/', sp_nm0,"_", "BIN_",proj_nm,"_",eval_stat,"_",ensemble_type,".jpg", sep = "")
      jpeg_name=paste(out_nm, ".jpg", sep = "")
      out_raster_name=paste(out_nm, ".jpg", sep = "")
      
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(get(raster_name_bin))
      dev.off()
      writeRaster(get(raster_name_bin), out_raster_name, format="GTiff", overwrite=TRUE)
      
    }  
    cat('\n','done with loading baseline and future rasters for ', sp_nm)
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


