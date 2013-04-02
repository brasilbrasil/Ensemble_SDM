rm(list = ls()) #remove all past worksheet variables
library(raster)
library(stringr)
library(biomod2)

###USER CONFIGURATION
overwrite=0
min_n_points=3
working_dir="Y:/plant SRE/" #this is where r will create a mirror folder structure with the 
setwd(working_dir)
spp_data=read.csv('all_spp_data.csv', ,header=T, stringsAsFactors=F)
all_spp=unique(spp_data$sp_code)
all_spp=sort(all_spp)
jnk=which(all_spp==0)
all_spp=all_spp[-jnk]
#unique_species

env_var_names=c("Tmax","Tmin","ppt") ###DEBUG!!

clim_data_dir="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/yearly_clim_data/2000/250m/" #this is the root directory of all the current env rasters
biovars2000 = stack(paste(clim_data_dir, "Yearly_clim_2000_250m_bio.tif", sep=""))
names(biovars2000)<- env_var_names

clim_data_dir="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/yearly_clim_data/2100/250m/" #this is the root directory of all the current env rasters
biovars2100 = stack(paste(clim_data_dir, "Yearly_clim_2100_250m_bio.tif", sep=""))

names(biovars2100)<- env_var_names

#sp=all_spp[1]
#sp=664
for (sp in all_spp){
  jnk=which(spp_data$sp_code==sp)
  sp_data=spp_data[jnk,]
  sp_name=sp_data$sp_name
  sp_name=sp_name[1]
  sp_str=as.character(sp)
  sp_str=str_pad(sp_str, 4, side = "left", pad = "0")
  ptm0 <- proc.time()
  cat('\n',sp_name,'modeling...')
  out_raster_name=paste(working_dir, "tifs/", sp_str, "_resp_envelopes_Tmin_Tmax_ppt.tif", sep="")
  if (file.exists(out_raster_name)==F | overwrite==1){
    sp_data=data.frame(cbind(x=sp_data$x, y=sp_data$y, pres=sp_data$pres))
    n_points=dim(sp_data)[1]
    if (n_points>=min_n_points){    
      Response_var=sp_data
      explanatory=extract(biovars2000, Response_var[,1:2])
      pred <- sre(Response_var, explanatory, biovars2000, Quant = 0.001)
      SRE_raster_present=subset(pred, 3)
      
      tifname=paste("tifs/", sp_str, "_SRE_2000_Tmin_Tmax_ppt.tif")
      writeRaster(SRE_raster_present, tifname, format="GTiff", overwrite=TRUE)
      
      #     plot_var="SRE_raster_present"    
      #       jpeg_name=paste("jpg_outputs/", sp_str, "_SRE_2000_7_vars.jpg", sep = "")
      #       jpeg(jpeg_name,
      #            width = 10, height = 10, units = "in",
      #            pointsize = 12, quality = 90, bg = "white", res = 300)
      #       plot(get(plot_var), legend=F)
      #     title(paste(sp_name, " current envelope"))
      #     dev.off()
      
      #future
      pred <- sre(Response_var, explanatory, biovars2100, Quant = 0.001)
      SRE_raster_future=subset(pred, 3)
      
      tifname=paste("tifs/", sp_str, "_SRE_2100_Tmin_Tmax_ppt.tif")
      writeRaster(SRE_raster_future, tifname, format="GTiff", overwrite=TRUE)
      
      #     plot_var="SRE_raster_future"    
      #     jpeg_name=paste("jpg_outputs/", sp_str, "_SRE_2100_7_vars.jpg", sep = "")
      #     jpeg(jpeg_name,
      #          width = 10, height = 10, units = "in",
      #          pointsize = 12, quality = 90, bg = "white", res = 300)
      #     plot(get(plot_var), legend=F)
      #     title(paste(sp_name, " future envelope"))
      #     dev.off()
      
      jnk=SRE_raster_future*10
      BIN_dif=SRE_raster_present+jnk
      m  =  c(9.9,  10.1,  3, 10.9, 11.1, 2)
      rclmat  =  matrix(m,  ncol=3,  byrow=TRUE)
      resp_zone  =  reclassify(BIN_dif,  rclmat)
      
      mypalette_numbers=c(0, 1, 2, 3)
      mypalette=c("Grey", "Red", "Green", "Yellow")
      resp_zone_names0=c("Micro refugia", "Tolerate", "Migrate")
      
      jnk=unique(resp_zone)
      zones_present=jnk[jnk>0]
      zones_present=zones_present[zones_present<=3]
      resp_zone_colors=mypalette[zones_present+1]
      resp_zone_names=resp_zone_names0[zones_present]
      
      jpeg_name=paste('jpg_outputs/', sp_str,"_response_zones_Tmin_Tmax_ppt.jpg", sep = "")
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300)
      plot(resp_zone,  col=mypalette, legend=F)
      #legend("bottomleft",legend = c("Micro refugia", "Tolerate", "Migrate"), col = mypalette[2:4],pch = 16)
      legend("bottomleft",legend = resp_zone_names, col = resp_zone_colors,pch = 16)
      title(paste(sp_name, " response zones (n=", n_points, ")", sep=""))
      dev.off()
      
      writeRaster(resp_zone, out_raster_name, format="GTiff", overwrite=TRUE) 
      
      ptm1=proc.time() - ptm0
      jnk=as.numeric(ptm1[3])
      jnk=jnk/60
      cat('\n','It took ', jnk, "minutes to model", sp_name)
    }
    
  }else{
    cat('\n', sp_name, " already calculated")
  }
  
} 
##load dynamic downscaling data
#clim_data_dir="Y:/SDM_env_data/bioclim_recalc/DD2100_2000_bioclimdeltas/" #this is the root directory of all the current env rasters
#Deltas2100to2000_DD_bioclim = stack(paste(clim_data_dir, "Delta_bioclimvars_2100_2000.tif", sep=""))



