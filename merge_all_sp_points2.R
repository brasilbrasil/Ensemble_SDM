rm(list = ls()) #remove all past worksheet variables
###USER CONFIGURATION
spp_nm=c('Iiwi', 'Kauai_Amakihi', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Oahu_Amakihi', 'Oahu_Elepaio', 'Hawaii_Amakihi', 'Puaiohi', 'Omao', 'Palila', 'Anianiau', 'Akekee', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Apapane', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Hawaii_Elepaio')
server=0
if (server==1){
  working_dir='Y:/FB SDM/biomod2/run_1_7_12_15/'
  necessary_run_data='Y:/FB SDM/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/test/'
  necessary_run_data='C:/Users/lfortini/Data/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
}
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")

###START UNDERHOOD
setwd(working_dir)
library(stringr)
dirs=list.dirs(necessary_run_data, full.names = FALSE, recursive = TRUE)
for (dir in dirs){
  layers<-list.files(dir, pattern=NULL, full.names=FALSE, include.dirs = FALSE)
  for (layer in layers){
    layer_full_nm=paste(dir,layer, sep="/")
    if (file.info(layer_full_nm)$isdir==FALSE){
      out_dir_nm=str_replace(dir, necessary_run_data, working_dir)
      dir.create(out_dir_nm, showWarnings = FALSE, recursive = TRUE, mode = "0777")
      out_lyr_nm=str_replace(layer_full_nm, necessary_run_data, working_dir)
      #out_lyr_nm=str_replace(out_lyr_nm, bl_lr, output_filenm)
      if (file.exists(out_lyr_nm)==F){
        cat('\n','found ', layer, 'in ', dir)
        file.copy(layer_full_nm, out_lyr_nm, overwrite = TRUE, recursive = TRUE,
                  copy.mode = TRUE)
        cat('\n','saved as ', out_lyr_nm)
      }
    }
  }
}

spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))

sp_nm=spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)
  cat('\n',sp_nm,'modeling...')  
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
  
  mySpeciesOcc=cbind(matrix(sp_nm,dim(mySpeciesOcc)[1],1),mySpeciesOcc)
  if (sp_nm==spp_nm[1]){
    all_mySpeciesOcc=mySpeciesOcc[,1:4]
  }else{
    all_mySpeciesOcc=rbind(all_mySpeciesOcc, mySpeciesOcc[,1:4])
  }
}

write.table(all_mySpeciesOcc, file = "all_pres_points.csv", sep=",", col.names=NA)
