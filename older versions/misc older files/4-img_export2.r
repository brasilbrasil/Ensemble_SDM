rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
local_config_dir='C:/Users/lfortini/' #if specifiying sp to run by file, this is directory of where csv file is located
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Akekee", "Akikiki", "Anianiau", "Oahu_Amakihi", "Puaiohi")   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
comp_projects=c('baseline', 'future') #put future second!
comp_eval_stats=c("ROC", "ROC")
comp_ensemble_type=c("ef.median", "ef.median")
#working_dir='C:/Users/lfortini/Forest bird SDM/biomod2/'
working_dir='Y:/FB SDM/biomod2/'
clim_data_dir0="Y:/SDM_env_data/bioclim_variables/FB_clim_var_data/mixed_data_2000_100m/"
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")
spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))

####START UNDERHOOD
output_dir=paste(working_dir,'output_rasters/', sep="")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
#working_dir='Y:/FB SDM/biomod2/'
all_vars=c("comp_projects", "comp_eval_stats", "comp_ensemble_type")

setwd(working_dir)
library(biomod2)
library(stringr)
#library(RColorBrewer)
##memory.limit(size=24000000)

#sp_nm=spp_nm[1]

for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  #workspace_name=paste(sp_nm,"_FB_EM_proj_", proj_nm, ".RData", sep = "")#set name of file to load workspace data from model run
  #load(workspace_name)
  sp_nm=str_replace_all(sp_nm,"_", ".")
  a=round(1000*runif(1))
  
  raster_names=c("EM_suitability1", "EM_suitability2")
  raster_names_bin=c("EM_BIN1", "EM_BIN2")
  for (i in c(1,2)){
    proj_nm=comp_projects[i]
    eval_stat=comp_eval_stats[i]
    raster_name=raster_names[i]
    raster_name_bin=raster_names_bin[i]
    ensemble_type=comp_ensemble_type[i]
    file_name1=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat, sep = "")
    file_name2=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.",eval_stat, sep = "")
    if (file.exists(file_name1)){ 
      load(file_name1)
      assign("temp_raster", get(paste(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat, sep = "")))           
    }else{
      load(file_name2)
      assign("temp_raster", get(paste(sp_nm,"_AllData_AllRun_EM.",eval_stat, sep = "")))      
    }

    file_name1_bin=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat,".bin.",eval_stat, sep = "")
    file_name2_bin=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_AllRun_EM.",eval_stat,".bin.",eval_stat, sep = "")
    if (file.exists(file_name1_bin)){       
    load(file_name1_bin)
    assign("temp_raster_bin", get(paste(sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat,".bin.",eval_stat, sep = "")))      
    }else{
    load(file_name2_bin)
    assign("temp_raster_bin", get(paste(sp_nm,"_AllData_AllRun_EM.",eval_stat,".bin.",eval_stat, sep = "")))      
    }
    
    band_names=names(temp_raster)
    band_n=which(band_names==ensemble_type)
    assign(raster_name, raster(temp_raster, layer=band_n))
    
    band_names=names(temp_raster_bin)
    band_n=which(band_names==paste(ensemble_type,".bin", sep = ""))
    assign(raster_name_bin, raster(temp_raster_bin, layer=band_n))
    
    #output suitability rasters for each image
    jpeg_name=paste('output_rasters/', a,"_",sp_nm0,"_", "suitability_",proj_nm,"_",eval_stat,"_",ensemble_type,".jpg", sep = "")
    jpeg(jpeg_name,
         width = 10, height = 8, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(get(raster_name))
    dev.off()
    
    out_raster_name=paste('output_rasters/', sp_nm0,"_", "suitability_",proj_nm,"_",eval_stat,"_",ensemble_type,".tif", sep = "")
    writeRaster(get(raster_name), out_raster_name, format="GTiff", overwrite=TRUE)
    
    #output bin rasters for each image
    jpeg_name=paste('output_rasters/', a,"_",sp_nm0,"_BIN_",proj_nm,"_",eval_stat,"_",ensemble_type,".jpg", sep = "")
    jpeg(jpeg_name,
         width = 10, height = 8, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(get(raster_name_bin))
    dev.off()
    
    out_raster_name=paste('output_rasters/', sp_nm0,"_BIN_",proj_nm,"_",eval_stat,"_",ensemble_type,".tif", sep = "")
    writeRaster(get(raster_name_bin), out_raster_name, format="GTiff", overwrite=TRUE)
  }
  
  ##COMPARISONS##
  Suitability_dif=EM_suitability2-EM_suitability1
  jnk=EM_BIN2*10
  BIN_dif=EM_BIN1+jnk
  m  =  c(9.9,  10.1,  3, 10.9, 11.1, 2)
  rclmat  =  matrix(m,  ncol=3,  byrow=TRUE)
  resp_zone  =  reclassify(BIN_dif,  rclmat)
  #mypalette<-brewer.pal(7,"Blues")
  #mypalette=mypalette[c(1, 4, 6, 7)]
  mypalette=c("Grey", "Red", "Green", "Yellow")

  current_mask=EM_suitability1>minValue(EM_suitability1)
  
  sp_index=which(spp_info[,"Species"]==sp_nm)
  raster_res= spp_info[sp_index,"rasterdir"]
  clim_data_dir=paste(clim_data_dir0,raster_res,"/grd/",sep="")
  jnk0=length(env_var_files)
  analog_cc = raster( paste(clim_data_dir, "analog_clim_2100_250m.tif", sep=""))
  habitat = raster( paste(clim_data_dir, "veg_250m.tif", sep=""))
  
  all_mask=current_mask+analog_cc+habitat #1 cur, 2 ang, 4 hab, 3 cur/ang, 6 ang/hab, 7 cur/ang/hab
  cum_mask=current_mask*analog_cc*habitat

  masked_resp_zone=resp_zone*cum_mask
  
  jpeg_name=paste(sp_nm0, "_mask_.jpg", sep = "")
  jpeg(jpeg_name,
       width = 10, height = 8, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  plot(all_mask)
  dev.off()
  
  
  fixed=""
  for (comp_vars in all_vars){
    jnk=get(comp_vars)
    if (jnk[1]!=jnk[2]){
      d_comp_vars=comp_vars
      jnk=get(comp_vars)
      comp_1=jnk[1]
      comp_2=jnk[2]
      }else{
        fixed=paste(fixed,jnk[1], sep="_")
    }
  }
  ##suitability comparison rasters
  jpeg_name=paste('output_rasters/', a,"_",sp_nm0, fixed, "_suitability_with_", comp_1,"_vs_", comp_2, ".jpg", sep = "")
  jpeg(jpeg_name,
       width = 10, height = 8, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  plot(Suitability_dif)
  dev.off()
  
  out_raster_name=paste('output_rasters/', sp_nm0,fixed, "_suitability_with_",comp_1,"_vs_", comp_2, ".tif", sep = "")
  writeRaster(Suitability_dif, out_raster_name, format="GTiff", overwrite=TRUE)
  
  if (comp_1=="baseline" & comp_2=="future"){
    ##bin comparison rasters
    jpeg_name=paste('output_rasters/', sp_nm0,"_response_zones_",eval_stat, "_", ensemble_type, ".jpg", sep = "")
    jpeg(jpeg_name,
         width = 10, height = 8, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(resp_zone,  col=mypalette, legend=F)
    #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
    legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
    dev.off()
    
    out_raster_name=paste('output_rasters/', sp_nm0,"_", "response_zones_",eval_stat, "_", ensemble_type,".tif", sep = "")
    writeRaster(resp_zone, out_raster_name, format="GTiff", overwrite=TRUE)

    ##MASKED bin comparison rasters
    jpeg_name=paste('output_rasters/', sp_nm0,"_masked_response_zones_",eval_stat, "_", ensemble_type, ".jpg", sep = "")
    jpeg(jpeg_name,
         width = 10, height = 8, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(masked_resp_zone,  col=mypalette, legend=F)
    #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
    legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
    dev.off()
    
    out_raster_name=paste('output_rasters/', sp_nm0,"_masked_response_zones_",eval_stat, "_", ensemble_type,".tif", sep = "")
    writeRaster(masked_resp_zone, out_raster_name, format="GTiff", overwrite=TRUE)
    
    }else{  
    ##bin comparison rasters
    jpeg_name=paste('output_rasters/', a,"_",sp_nm0, fixed, "_BIN_with_", comp_1,"_vs_", comp_2,".jpg", sep = "")
    jpeg(jpeg_name,
         width = 10, height = 8, units = "in",
         pointsize = 12, quality = 90, bg = "white", res = 300)
    plot(resp_zone,  col=mypalette, legend=F)
    #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
    #legend("bottomright", "Micro refugia", "Tolerate", "Migrate")
    dev.off()
    
    out_raster_name=paste('output_rasters/', sp_nm0, fixed, "_BIN_with_", comp_1,"_vs_", comp_2,".tif", sep = "")
    writeRaster(resp_zone, out_raster_name, format="GTiff", overwrite=TRUE)
  }
}



