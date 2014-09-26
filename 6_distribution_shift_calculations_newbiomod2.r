clim_data_dir = clim_data_2000 
current_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/current_veg_mask/"
projected_biome_distribution_dir="Y:/PICCC_analysis/FB_analysis/habitat analysis/veg_overlay/projected_veg_mask/"

####START UNDERHOOD
library(biomod2)
library(stringr)
dir.create("tables/", showWarnings=F)
px_area=model_resolution^2 
vars=c("Eval stat", "Species", "baseline_area", "future_area", "% change", "lost", "kept", "gained", "baseline_suitability", "future_suitability", "% suitability change")
all_stats<- as.data.frame(setNames(replicate(length(vars),numeric(0), simplify = F), vars)) #create empty dataframe to compile results
if (exclude_areas_beyond_primary_habitat){
  mask_text="habitat_mask_"
}else{
  mask_text=""
}

sp_nm=spp_nm[1]
eval_stat = spp_ensemble_eval_stats[1]
for (eval_stat in spp_ensemble_eval_stats){
  for (sp_nm in spp_nm){
    sp_nm=as.character(sp_nm)  
    cat('\n',sp_nm,'modeling...')
    sp_nm0=sp_nm
    sp_nm=str_replace_all(sp_nm,"_", ".")
    a=round(1000*runif(1))
    
    response_raster_nm=paste('output_rasters/main/', sp_nm0,"_response_zones_",eval_stat, "_", spp_ensemble_type, sep = "")
    response_raster_nm=paste(response_raster_nm,".tif", sep = "")
    response_raster=raster(response_raster_nm)
    
    if (exclude_areas_beyond_primary_habitat){
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
    }
    
    jnk=response_raster>0
    jnk=zonal(jnk,response_raster,sum,na.rm=T)
    jnk=matrix(jnk[jnk[,1]>0,],ncol=2)
    jnk[,2]=jnk[,2]*px_area
    zone_area=as.data.frame(matrix(c(1,2,3,0,0,0),nrow=3,ncol=2))
    names(zone_area)=c("zone", "area")
    zone_area[jnk[,1],2]=jnk[,2]
    zone_area=zone_area[,2]
    curr=zone_area[1]+zone_area[2]
    fut=zone_area[2]+zone_area[3]
    results=c(curr,fut,(fut-curr)*100/curr,zone_area)
    
    #masked suitability
    baseline_masked_suitability=raster(paste('output_rasters/', sp_nm0,"_", "clipped_suitability_",comp_projects[1],"_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
    future_masked_suitability=raster(paste('output_rasters/', sp_nm0,"_", "clipped_suitability_",comp_projects[2],"_",eval_stat,"_",spp_ensemble_type, ".tif", sep = ""))
    if (exclude_areas_beyond_primary_habitat){
      baseline_masked_suitability=baseline_masked_suitability*Response_var1
      future_masked_suitability=future_masked_suitability*Response_var1
    }
    baseline_masked_suitability[baseline_masked_suitability==0]=NA
    future_masked_suitability[future_masked_suitability==0]=NA
    jnk1=cellStats(baseline_masked_suitability, stat='mean', na.rm=TRUE)
    jnk2=cellStats(future_masked_suitability, stat='mean', na.rm=TRUE)
    results=matrix(c(results, jnk1, jnk2, 100*(jnk2-jnk1)/jnk1),nrow=1)
    results=as.data.frame(results)
    results=cbind(eval_stat, sp_nm0,results)
    names(results)=vars
    out_csv_name=paste('tables/', sp_nm0,"_metrics_",mask_text,eval_stat, "_", spp_ensemble_type, ".csv",sep = "")
    write.table(results, file = out_csv_name, sep=",", row.names=F)
    all_stats=rbind(all_stats,results)
    
  }
}
if (exclude_areas_beyond_primary_habitat){
  out_csv_name=paste('tables/', "all_spp_distr_shift_metrics_current_habitat.csv",sep = "")
  write.table(all_stats, file = out_csv_name, sep=",", row.names=F)
  
}else{
  out_csv_name=paste('tables/', "all_spp_distr_shift_metrics.csv",sep = "")
  write.table(all_stats, file = out_csv_name, sep=",", row.names=F)
  
}


