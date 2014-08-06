rm(list = ls()) #remove all past worksheet variables
#for each species modeled, have csv of presence data in working directory for the species named speciesname_Ps.csv formated with 3 cols: x,y,pa where pa = 1
#after running the code for whichever many species, copy results (species output folder and workspace file) to a new directory, along with the maxent.jar file
#use the projection code to project the distribution model on different environmental surfaces (do not forget to change the working directory)

###USER CONFIGURATION
#local_config_dir='C:/Users/lfortini/'
local_config_dir='Y:/FB_analysis/FB_SDM/biomod2/' #'C:/Users/lfortini/'
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_all.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Kauai_Amakihi")#, "Oahu_Amakihi", "Apapane")# "Akekeke", "Akikiki", "Anianiau", "Kauai_Amakihi", "Kauai_Elepaio", "Puaiohi", "Oahu_Amakihi", "Oahu_Elepaio", "Apapane", "Iiwi",)   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
server=1
remove_PA_abs=TRUE
overwrite=0

models_to_run=c('GBM','RF','MAXENT')
if (server==1){
  working_dir='Y:/FB_analysis/FB_SDM/biomod2/response_curves_500/'
  clim_data_dir0="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/100m/" 
  necessary_run_data='Y:/FB_analysis/FB_SDM/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/test/'
  necessary_run_data='C:/Users/lfortini/Data/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
  clim_data_dir0="C:/Users/lfortini/Data/SDM_env_data/all_grd/all_baseline/100m/"
}

csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")


###START UNDERHOOD
setwd(working_dir)
library(stringr)

spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))

#sp_nm=spp_nm[1]
#spp_nm=spp_nm[11:length(spp_nm)]
#spp_nm=spp_nm[1:10]
#spp_nm=rev(spp_nm)
n_abs_removed=c()

sp_nm = spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)
  cat('\n',sp_nm,'modeling...')
  
  # Start the clock!
  ptm0 <- proc.time()
  workspace_name=paste(sp_nm,"_FB_run.RData", sep = "") #set name of file to save all workspace data after model run
  
  #######Loading datasets#######
  mySpeciesOcc=read.csv(paste(csv_dir,sp_nm,'_pres_abs.csv', sep = "")) #FB_data_points4_PAandA
  
  # Select, Count and Remove presence Duplicates
  jnk=dim(mySpeciesOcc)[1]
  dups2<- duplicated(mySpeciesOcc[, c('X','Y')])
  sum(dups2)
  mySpeciesOcc<-mySpeciesOcc[!dups2, ]
  jnk1=dim(mySpeciesOcc)[1]
  mySpeciesOcc=mySpeciesOcc[mySpeciesOcc[,"pa"]==1,] #get rid of absences
  jnk2=dim(mySpeciesOcc)[1]
  head(mySpeciesOcc)
  cat('\n','removed ', jnk-jnk1, "duplicates for", sp_nm)
  cat('\n','removed ', jnk1-jnk2, "absence records for", sp_nm)
  n_P=jnk2
  n_A=jnk1-jnk2
  tmp=c(sp_nm,n_P,n_A)
  if (sp_nm == spp_nm[1]){
    tmp_all=tmp
  }else{
    tmp_all=rbind(tmp_all, tmp)
  }
}
write.table(tmp_all, file = "all_spp_n.csv", sep=",", col.names=NA)
    