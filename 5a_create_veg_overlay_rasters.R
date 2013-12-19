rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))

###USER CONFIGURATION
working_dir='Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/BPS/'
clim_data_dir=clim_data_2000 
overwrite=0 #if 1, will overwrite past results
projected_biome=F
BPS=T

####START UNDERHOOD
setwd(working_dir)
library(biomod2)
library(stringr)

biome_map_dir="//10.0.0.5/data2$//PICCC_analysis/community_SRE/"
if (BPS){
  landfire_reclass_resampled=raster(paste(biome_map_dir, "landfire_BPS_reclass_500m.tif", sep=""))  
}else{
  landfire_reclass_resampled=raster(paste(biome_map_dir, "landfire_reclass_wetland_coastal_500m.tif", sep=""))  
}
biome_proj_map_dir="//10.0.0.5/data2$//PICCC_analysis/community_SRE/tifs/"

biome_association_table=(read.csv('bird_biome_association.csv',header=T, stringsAsFactors=F))
dir.create("current_veg_mask/",showWarnings=FALSE)
dir.create("projected_veg_mask/",showWarnings=FALSE)

spp_nm=biome_association_table[,'Species']# "Akekeke", "Akikiki", "Anianiau", "Kauai_Amakihi", "Kauai_Elepaio", "Puaiohi", "Oahu_Amakihi", "Oahu_Elepaio", "Apapane", "Iiwi",)   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"


sp_nm=spp_nm[4]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  out_raster_name=paste("projected_veg_mask/",sp_nm,"_projected_veg_mask", sep="")
  jpeg_name=paste(out_raster_name, ".jpg", sep="")
  if (file.exists(jpeg_name)==F | overwrite==1){ #check to see if the analysis for this species was already done    
    
    sp_biome_associations=biome_association_table[biome_association_table[,'Species']==sp_nm,'biome_association']
    sp_biome_associations=c(unlist(strsplit(sp_biome_associations, split=", ")))
    sp_biome_associations=as.numeric(sp_biome_associations)
    if (any(sp_biome_associations == 14)){
      sp_biome_associations=sp_biome_associations[sp_biome_associations!=14]
      sp_biome_associations=c(sp_biome_associations, 7, 11, 12)
    }
    
    all_current_biome=landfire_reclass_resampled %in% sp_biome_associations
    #plot(all_current_biome)
    
    out_raster_name=paste("current_veg_mask/",sp_nm,"_current_veg_mask", sep="")
    out_tif_name=paste(out_raster_name, ".tif", sep="")
    writeRaster(all_current_biome, out_tif_name, format="GTiff", overwrite=TRUE)
    
    jpeg_name=paste(out_raster_name, ".jpg", sep="")
    jpeg(jpeg_name,
         width = 10, height = 8, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(all_current_biome, legend=F)
    title(paste(sp_nm, " current distribution of primary biomes", sep=""))
    dev.off()  
    
    ##projected biome associations
    if (projected_biome){
      
    sp_biome_association=sp_biome_associations[1]
    sp_biome_association=sp_biome_associations[2]
    for (sp_biome_association in sp_biome_associations){
      sp_biome_association_map=paste(biome_proj_map_dir,sp_biome_association,"_resp_envelopes_bio5_6_12_Q99.tif", sep="")
      sp_biome_association_map=raster(sp_biome_association_map)
      temp_tol=sp_biome_association_map==2
      temp_mig=sp_biome_association_map==3
      if (sp_biome_association==sp_biome_associations[1]){
        tol_assoc_ensemble=temp_tol
        mig_assoc_ensemble=temp_mig
      }else{
        tol_assoc_ensemble=tol_assoc_ensemble+temp_tol
        mig_assoc_ensemble=mig_assoc_ensemble+temp_mig
      }
    }
    tol_assoc_ensemble=(tol_assoc_ensemble>0)*10
    mig_assoc_ensemble=mig_assoc_ensemble>0
    combined_projection_map=tol_assoc_ensemble+mig_assoc_ensemble
    combined_projection_map[combined_projection_map==11]=10
    combined_projection_map[combined_projection_map==10]=2
    combined_projection_map[combined_projection_map==1]=3
    #unique(combined_projection_map)
    
    out_raster_name=paste("projected_veg_mask/",sp_nm,"_projected_veg_mask", sep="")
    out_tif_name=paste(out_raster_name, ".tif", sep="")
    writeRaster(combined_projection_map, out_tif_name, format="GTiff", overwrite=TRUE)
    
    jpeg_name=paste(out_raster_name, ".jpg", sep="")
    jpeg(jpeg_name,
         width = 10, height = 8, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(combined_projection_map, legend=F)
    title(paste(sp_nm, " projected overlap and gained climate envelope for primary biomes", sep=""))
    dev.off()  
    }
    
  }
}


