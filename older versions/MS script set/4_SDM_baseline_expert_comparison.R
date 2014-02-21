rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))
###USER CONFIGURATION
spp_nm = c('Akekee', 'Hawaii_Amakihi', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi','Kauai_Amakihi', 'Iiwi', 'Amakihi', 'Elepaio', 'Apapane')
project_name='finalmodel_P_PA_oldcode_less_PAs'

comp_project='baseline' #put future second!
eval_stats=c('ROC') 
plot_CV=T
#eval_stats=c('ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS') 
masked=FALSE
overwrite=T

working_dir=paste0(resultsDir,project_name,'/')
clim_data_dir=paste0(bioclimData2013Dir,"all_baseline/500m/")
overwrite=0 #if 1, will overwrite past results
current_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/current_veg_mask/"
projected_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/projected_veg_mask/"
veg_areas_loc=paste0(clim_data_dir, "veg_areas.grd")

####START UNDERHOOD
setwd(working_dir)
library(biomod2)
library(stringr)

save_raster_fx=function(raster_img,out_nm, mask_data=NULL, expert_data=NULL){
  jpeg_name=paste(out_nm, ".jpg", sep = "")
  out_raster_name=paste(out_nm, ".tif", sep = "")
  jpeg(jpeg_name,
       width = 10, height = 8, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  plot(raster_img, legend=F)
  if (!is.null(mask_data)){
    plot(mask_data,add=T)    
  }
  if (!is.null(expert_data)){
    plot(expert_data,add=T, border="red", lwd=3)    
  }
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

Process_raster_data_NeutraltoGood=function(raster_var,out_nm,min_lim=NULL, max_lim=NULL, mask_data=NULL, expert_data=NULL){
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
  if (!is.null(expert_data)){
    plot(expert_data,add=T, border="red", lwd=3)    
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

dir.create('output_rasters/',showWarnings=F)
dir.create('output_rasters/expert_comparisons/',showWarnings=F)
sp_nm=spp_nm[1]
eval_stat=eval_stats[1]
for (eval_stat in eval_stats){
  for (sp_nm in spp_nm){
    sp_nm=as.character(sp_nm)  
    cat('\n',sp_nm,'modeling...')
    sp_nm0=sp_nm
    sp_nm=str_replace_all(sp_nm,"_", ".")
    
    out_nm=paste('output_rasters/expert_comparisons/', sp_nm0,"_", "Expert_vs_clipped_suitability_",comp_project,"_",eval_stat, ".tif",sep = "")
    if (file.exists(out_nm)==F | overwrite){
      
      shp_dir="Y:/PICCC_analysis/FB_analysis/FB range expert maps/shapefiles/"
      expert_shp=shapefile(paste0(shp_dir, sp_nm0,"_range.shp"))
      
      ##binary maps
      raster_name="EM_suitability1"
      raster_name_bin="EM_BIN1"
      proj_nm=comp_project
      
      
      
      file_name1=paste(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_",sp_nm,"_ensemble.grd", sep = "")
      if (file.exists(file_name1)){
        ensemble_type="wmean"
        temp_raster=stack(file_name1)
        band_n=which(names(temp_raster)==paste0(sp_nm,"_TotalConsensus_",eval_stat,"_EM",ensemble_type))
        assign(raster_name, raster(temp_raster, layer=band_n)/1000)
        
        file_name1_bin=paste(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_",sp_nm,"_ensemble_",eval_stat,"bin.grd", sep = "")
        temp_raster_bin=stack(file_name1_bin)  
        band_n=which(names(temp_raster)==paste0(sp_nm,"_TotalConsensus_",eval_stat,"_EM",ensemble_type))
        assign(raster_name_bin, raster(temp_raster_bin, layer=band_n))
        
      }else{
        ensemble_type="ef.pmw"
        file_name1=paste(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_",sp_nm,"_TotalConsensus_EMby",eval_stat,".grd", sep = "")
        temp_raster=stack(file_name1)
        band_n=which(names(temp_raster)==paste0(sp_nm,"_TotalConsensus_EMby",eval_stat,"_",ensemble_type))
        assign(raster_name, raster(temp_raster, layer=band_n)/1000)
        
        file_name1_bin=paste(sp_nm,"/proj_", proj_nm, "/proj_", proj_nm, "_",sp_nm,"_TotalConsensus_EMby",eval_stat,"_",eval_stat,"bin.grd", sep = "")
        temp_raster_bin=stack(file_name1_bin)  
        band_n=which(names(temp_raster_bin)==paste0(sp_nm,"_TotalConsensus_EMby",eval_stat,"_",ensemble_type,"_",eval_stat,"bin"))
        assign(raster_name_bin, raster(temp_raster_bin, layer=band_n))      
      }
      
      
      #output bin rasters for each image    
      out_nm=paste('output_rasters/expert_comparisons/', sp_nm0,"_", "Expert_vs_BIN_",proj_nm,"_",eval_stat, sep = "")
      save_raster_fx(get(raster_name_bin),out_nm, mask_data=mask_layer, expert_data=expert_shp)
      
      #masked suitability
      masked_suitability1=EM_BIN1*EM_suitability1
      out_nm=paste('output_rasters/expert_comparisons/', sp_nm0,"_", "Expert_vs_clipped_suitability_",comp_project,"_",eval_stat, sep = "")
      #save_raster_fx(masked_suitability1,out_nm)
      Process_raster_data_NeutraltoGood(masked_suitability1,out_nm,min_lim=0, max_lim=1, mask_data=mask_layer, expert_data=expert_shp)
      
      cat('\n','done with loading baseline and future rasters for ', sp_nm)
    }else{
      cat('\n', sp_nm, " already calculated")
    }
  }
}


