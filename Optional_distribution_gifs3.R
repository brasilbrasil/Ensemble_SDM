rm(list = ls()) #remove all past worksheet variables
library(jpeg)
library(graphics)
library(animation)
###USER CONFIGURATION
local_config_dir='Y:/FB_analysis/FB_SDM/biomod2/' #'C:/Users/lfortini/'
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Iiwi")
comp_projects=c('baseline', 'future') #put future second!
ensemble_type="ef.pmw"
avail_models=c('GBM', 'RF', 'MAXENT')
n_random_pics=20
n_models=100 #number of total runs
#working_dir='Y:/FB_analysis/FB_SDM/biomod2/'
working_dir='Y:/FB_analysis/FB_Base_And_Future/all_DD_merged/'
eval_stat="ROC"
clim_data_dir="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/250m/"
masked=FALSE
overwrite=0 #if 1, will overwrite past results
gif=TRUE #if false, make mpeg #mpeg currently not working, need to instal mpeg2encode

####START UNDERHOOD
setwd(working_dir)
dir.create("output_rasters")
dir.create("output_rasters/gifs")
library(biomod2)
library(stringr)
if (gif){
  ext=".gif"
}else{
  ext=".mpeg"  
}

sp_nm=spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  sp_nm=str_replace_all(sp_nm,"_", ".")
  ii=1
  for (ii in 1:n_random_pics){ 
    a=round((n_models-1)*runif(1))+1
    picked_model=avail_models[sample(1:3, 1)]
    out_nm=paste('output_rasters/gifs/', sp_nm0,"_response_zones_",eval_stat, "_", ensemble_type, sep = "")
    out_raster_name00=paste(out_nm,".tif", sep = "")
    jpeg_name00=paste(out_nm, ii ,"_tempimg.jpg", sep = "")
    if (file.exists(jpeg_name00)==F | overwrite==1){
      
      raster_names=c("EM_suitability1", "EM_suitability2")
      raster_names_bin=c("EM_BIN1", "EM_BIN2")
      i=1
      for (i in c(1,2)){
        proj_nm=comp_projects[i]
        raster_name=raster_names[i]
        raster_name_bin=raster_names_bin[i]
        file_name1=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_RUN",a,"_","AllAlgos_EM.",eval_stat, sep = "")
        load(file_name1)
        assign("temp_raster", get(paste(sp_nm,"_AllData_RUN",a,"_","AllAlgos_EM.",eval_stat, sep = "")))           
        
        file_name1_bin=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_RUN",a,"_","AllAlgos_EM.",eval_stat,".bin.",eval_stat, sep = "")#file_name1_bin=paste(sp_nm,"/proj_", proj_nm, "/",sp_nm,"_AllData_Full_AllAlgos_EM.",eval_stat,".bin.",eval_stat, sep = "")
        load(file_name1_bin)
        assign("temp_raster_bin", get(paste(sp_nm,"_AllData_RUN",a,"_","AllAlgos_EM.",eval_stat,".bin.",eval_stat, sep = "")))      
        
        band_names=names(temp_raster)
        band_n=which(band_names==ensemble_type)
        assign(raster_name, raster(temp_raster, layer=band_n))
        
        band_names=names(temp_raster_bin)
        band_n=which(band_names==paste(ensemble_type,".bin", sep = ""))
        assign(raster_name_bin, raster(temp_raster_bin, layer=band_n))    
        
        #output suitability rasters for each image
        out_nm=paste('output_rasters/gifs/', sp_nm0,"_", "suitability_",proj_nm,"_",eval_stat,"_",ensemble_type, sep = "")
        jpeg_name=paste(out_nm, ii ,"_tempimg.jpg", sep = "")
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        #plot(get(raster_name),main=NULL) #col=rev(rainbow(99, start=0,end=1)),
        #plot(get(raster_name), col=rainbow(99, start=0,end=1), breaks=seq(0,1000,length.out=101),main=NULL)
        plot(get(raster_name), col = rev(terrain.colors( length(seq(0, 1000, by = 50))-1)), axes = FALSE, breaks= seq(0, 1000, by = 50),main=NULL) 
        dev.off()    
        
        #output bin rasters for each image    
        out_nm=paste('output_rasters/gifs/', sp_nm0,"_", "BIN_",proj_nm,"_",eval_stat,"_",ensemble_type, sep = "")
        jpeg_name=paste(out_nm, ii ,"_tempimg.jpg", sep = "")
        
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        #plot(get(raster_name_bin),main=NULL)
        plot(get(raster_name_bin), col = rev(terrain.colors( length(seq(0, 1, by = 0.05))-1)), axes = FALSE, breaks= seq(0, 1, by = 0.05),main=NULL) 
        dev.off()
        
        #output suitability bin rasters for each image    
        out_nm=paste('output_rasters/gifs/', sp_nm0,"_", "threshold_suit_",proj_nm,"_",eval_stat,"_",ensemble_type, sep = "")
        jpeg_name=paste(out_nm, ii ,"_tempimg.jpg", sep = "")
        
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        #plot(get(raster_name_bin),main=NULL)
        jnkk=get(raster_name_bin)*get(raster_name)
        plot(jnkk, col = rev(terrain.colors( length(seq(0, 1000, by = 50))-1)), axes = FALSE, breaks= seq(0, 1000, by = 50),main=NULL) 
        dev.off()
        
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
      
      
      ##bin comparison rasters
      out_nm=paste('output_rasters/gifs/', sp_nm0,"_response_zones_",eval_stat, "_", ensemble_type, sep = "")
      
      jnk=unique(resp_zone)
      graph_palette=mypalette_numbers
      zones_present=jnk[jnk>0]
      zones_present=zones_present[zones_present<=3]
      resp_zone_colors=mypalette[zones_present+1]
      resp_zone_names=resp_zone_names0[zones_present]
      mypalette_numbers_selected=mypalette[jnk+1] #CHANGED
      
      jpeg(jpeg_name00,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(resp_zone,  col=mypalette_numbers_selected, legend=F, main=NULL)
      #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
      legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
      dev.off()
    }else{
      cat('\n', sp_nm, " already calculated")
    }
  }
  
  outgif_dir=paste(working_dir, 'output_rasters/gifs/', sep="")
  proj_nm=comp_projects[2]
  output_file00=paste(outgif_dir,sp_nm0,"_", "threshold_suit_",proj_nm,"_",eval_stat,"_",ensemble_type,ext, sep = "")#"animation.gif"
  if (file.exists(output_file00)==F | overwrite==1){
    ani.options(convert = shQuote('C:/Program Files/ImageMagick-6.8.2-Q16/convert'), outdir=outgif_dir, autobrowse=FALSE)
    
    out_nm=paste('output_rasters/gifs/', sp_nm0,"_response_zones_",eval_stat, "_", ensemble_type, sep = "")
    jpeg_name=paste(out_nm, "*_tempimg.jpg", sep = "")
    output_file=paste(sp_nm0,"_response_zones_",eval_stat, "_", ensemble_type,ext, sep = "")#"animation.gif"
    im.convert(jpeg_name, output = output_file, convert = c("convert", "gm convert"),
               cmd.fun  =  system,  extra.opts  =  "-delay 1x5",  clean  =  FALSE)
    
    for (i in c(1,2)){
      proj_nm=comp_projects[i]
      
      out_nm=paste('output_rasters/gifs/', sp_nm0,"_", "suitability_",proj_nm,"_",eval_stat,"_",ensemble_type, sep = "")
      jpeg_name=paste(out_nm, "*_tempimg.jpg", sep = "")
      output_file=paste(sp_nm0,"_", "suitability_",proj_nm,"_",eval_stat,"_",ensemble_type,ext, sep = "")#"animation.gif"
      im.convert(jpeg_name, output = output_file, convert = c("convert", "gm convert"),
                 cmd.fun  =  system,  extra.opts  =  "-delay 1x5",  clean  =  FALSE)
      
      out_nm=paste('output_rasters/gifs/', sp_nm0,"_", "BIN_",proj_nm,"_",eval_stat,"_",ensemble_type, sep = "")
      jpeg_name=paste(out_nm, "*_tempimg.jpg", sep = "")
      output_file=paste(sp_nm0,"_", "BIN_",proj_nm,"_",eval_stat,"_",ensemble_type,ext, sep = "")#"animation.gif"
      im.convert(jpeg_name, output = output_file, convert = c("convert", "gm convert"),
                 cmd.fun  =  system,  extra.opts  =  "-delay 1x5",  clean  =  FALSE)
      
      out_nm=paste('output_rasters/gifs/', sp_nm0,"_", "threshold_suit_",proj_nm,"_",eval_stat,"_",ensemble_type, sep = "")
      jpeg_name=paste(out_nm, "*_tempimg.jpg", sep = "")
      output_file=paste(sp_nm0,"_", "threshold_suit_",proj_nm,"_",eval_stat,"_",ensemble_type,ext, sep = "")#"animation.gif"
      im.convert(jpeg_name, output = output_file, convert = c("convert", "gm convert"),
                 cmd.fun  =  system,  extra.opts  =  "-delay 1x5",  clean  =  FALSE)
      
    }
  }else{
    cat('\n',output_file00,' previously done...')
  }
}


