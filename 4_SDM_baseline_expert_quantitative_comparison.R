rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))
###USER CONFIGURATION
spp_nm = c('Akekee', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Iiwi', 'Maui_Alauahio', 'Maui_Parrotbill', 'Puaiohi', 'Anianiau', 'Apapane', 'Hawaii_Amakihi', 'Hawaii_Elepaio', 'Kauai_Amakihi', 'Kauai_Elepaio', 'Oahu_Amakihi', 'Oahu_Elepaio', 'Omao', 'Palila')
Reliability=c('High', 'High', 'High', 'High', 'High', 'High', 'High', 'High', 'High', 'High', 'Low', 'Low', 'Low', 'Low', 'Low', 'Low', 'Low', 'Low', 'Low', 'Low')
project_name='finalmodel_P_PA_oldcode_less_PAs'

comp_project='baseline' #put future second!
eval_stats=c('ROC') 
plot_CV=T
#eval_stats=c('ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS') 
masked=FALSE

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

make_confusion_report=function(projected_cover,Response_var0){
  projected_cover_masked=projected_cover #getting rid of 'other' areas
  eval_actual=Response_var0
  actual=c(as.matrix(eval_actual))
  observed=c(as.matrix(projected_cover_masked))
  nonNA_index=!is.na(actual)
  ref=actual[nonNA_index]
  pred=observed[nonNA_index]
  cover_index=ref>=0
  ref=ref[cover_index]
  pred=pred[cover_index]    
  require(caret)
  a=confusionMatrix(pred, ref)
  return(a)
}

all_accuracy = data.frame(matrix(vector(), 0, 2, dimnames=list(c(), c("species", "accuracy"))), stringsAsFactors=F)
dir.create('output_rasters/',showWarnings=F)
dir.create('output_rasters/expert_comparisons/',showWarnings=F)
sp_nm=spp_nm[2]
eval_stat=eval_stats[1]
for (eval_stat in eval_stats){
  for (sp_nm in spp_nm){
    sp_nm=as.character(sp_nm)  
    cat('\n',sp_nm,'modeling...')
    sp_nm0=sp_nm
    sp_nm=str_replace_all(sp_nm,"_", ".")
    
    out_nm=paste('output_rasters/expert_comparisons/', sp_nm0,"_", "Expert_vs_clipped_suitability_",comp_project,"_",eval_stat, ".txt",sep = "")
    expertRaster_dir="Y:/PICCC_analysis/FB_analysis/FB range expert maps/tifs only latlon comparison/aggregate_median/reprojected/setNull/"
    expert_raster=raster(paste0(expertRaster_dir, sp_nm0,"_aggr.tif"))
    
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
    
    #align rasters
    model_raster=EM_BIN1
    expert_raster=resample(expert_raster, model_raster,  method="ngb")
    expert_raster=expert_raster*(model_raster>=0) #get rid of islands not modeled
    #plot(model_raster)
    #plot(expert_raster)      
    
    a=make_confusion_report(model_raster,expert_raster)
    sink(file(out_nm, open="wt"))
    print(a)
    sink(NULL)
    accuracy=a[3][[1]][1]
    jnk=data.frame(species=sp_nm, accuracy=accuracy)
    all_accuracy=rbind(all_accuracy,jnk)
    
    cat('\n','done with ', sp_nm)
  }
}

all_accuracy=cbind(all_accuracy,Reliability)
library(ggplot2)
long=all_accuracy
csv_name=paste0('output_rasters/expert_comparisons/', "Expert_vs_clipped_suitability_aggreement.csv")
jpeg_name=paste0('output_rasters/expert_comparisons/', "Expert_vs_clipped_suitability_",comp_project,"_",eval_stat, ".jpg")
a=qplot(Reliability, accuracy, data=long, geom=c("boxplot", "jitter"), 
        fill=Reliability, main="",
        xlab="Model reliability", ylab="Model expert agreement") #, ylim=c(0,YLim[p])
ggsave(filename=jpeg_name, plot=a)
write.csv(all_accuracy, csv_name, row.names=F)

csv_name=paste0('output_rasters/expert_comparisons/', "Expert_vs_clipped_suitability_aggreement_test.txt")
sink(file(csv_name, open="wt"))
t.test(accuracy~Reliability, data=long)
sink(NULL)
