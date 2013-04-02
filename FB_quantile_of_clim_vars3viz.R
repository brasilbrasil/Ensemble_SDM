
###USER CONFIGURATION
local_config_dir='C:/Users/lfortini/'
#local_config_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/' #'C:/Users/lfortini/'
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Akekee", "Akikiki", "Anianiau")#, "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
server=0
overwrite=1
proj_names=c("Pres", "Pres_Abs")
Quant = 0.025
Quant_val=(1-Quant*2)*100

if (server==1){
  working_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/'
  clim_data_dir0="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/100m/" 
  necessary_run_data='Y:/FB analysis/FB SDM/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/response_curves/'
  clim_data_dir0="C:/Users/lfortini/Data/SDM_env_data/all_grd/all_baseline/100m/"
  necessary_run_data='C:/Users/lfortini/Data/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
}

###START UNDERHOOD
setwd(working_dir)

proj_name=proj_names[2]
for (proj_name in proj_names){
  jnk<-paste("quantiles/Q", Quant_val, "_", proj_name, "_all_spp.csv", sep="")  
  jnk2=read.csv(jnk,header=T, stringsAsFactors=F)
  if (proj_name==proj_names[1]){
    merged_data=jnk2  
  }else{
    merged_data=  rbind(merged_data, jnk2)
  }
}
####SAVE TABLE
FileName=paste("quantiles/Q", Quant_val, "_p_vs_pa_all_spp.csv", sep="")
write.table(merged_data, file = FileName, sep=",", col.names=NA)



