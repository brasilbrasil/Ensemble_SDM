rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))
#options(error=stop) #this keeps the code from running after errors 
###USER CONFIGURATION
#spp_nm = c('Akekee', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi', 'Kauai_Amakihi', 'Hawaii_Amakihi', 'Apapane', 'Amakihi', 'Elepaio', 'Iiwi')
spp_nm = c('Akekee', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi', 'Kauai_Amakihi', 'Hawaii_Amakihi', 'Apapane', 'Iiwi')

project_name='finalmodel_P_PA_oldcode_less_PAs'
server=1
overwrite=0; paralelize=F
plot_graphs=1
EM_fit=T
EM_ensemble=T
EM_project=T
create_response_curves=F
apply_biomod2_fixes=T #if running large models use this option

if (server==1){
  working_dir=paste0(resultsDir,project_name,'/')
  #fitting_clim_data_dir=paste0(DR_FB_clim_data_2013,"all_baseline/125m/") 
  #necessary_run_data=paste0(resultsDir,'necessary_run_data/') #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.)
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/test/'
  necessary_run_data='C:/Users/lfortini/Data/biomod2/necessary_run_data/' #where all needed files are stored (maxent.jar, species csvs, crop rasters, etc.    
  fitting_clim_data_dir="C:/Users/lfortini/Data/SDM_env_data/all_grd/all_baseline/100m/"
  DR_code_S=paste0(rootDir, "/Dropbox/code/") #HAD TO ADD THIS TO READ THIS FROM OLD DIRECTORY FILE 
  
}

env_var_files=c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif") 
crop_raster_dir=paste(working_dir, 'map_crop/',sep="")
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")

#################################
####CONFIG FOR SPECIFIC STEPS####
#################################
####fit config (script#1)
NbRunEval=5
include_Abs=F #in test phase
PAs_outside_CE=F #if T, will only consider PAs outside climate envelope of all points collected
dens_PAs_outside_CE=1 #if 1 will create PA density that is equal to point density within surveyed areas
PA.nb.rep=40
PA.nb.absences = 10000 #only used if if PAs_outside_CE=F, this will be overridden! (n of PAs will be determined by P/A point density within CE 
candidatePAperPA=100 #only used if if PAs_outside_CE=F, if value ==0, will use PA.nb.absences   
PA.strategy = "random"
equiv_100m=0.0009430131
PA.dist.min = 5*equiv_100m #500 min distance from actual data points 
do.full.models=T
####ensemble config (script#2)

####projection config (script#3)
baseline_or_future=1 #1 for baseline, 4 for future
memory = T #keep.in.memory=memory
dir_for_temp_files<-paste(rootDir,'/temp/', project_name,'/', baseline_or_future, '/', sep='') #dir for temp run data (to avoid memory errors)

if (server==1){
  clim_data_2000wettest="D:/GIS_Data/REnviroLayers/mixed_data_2000_250mwettest/"
  clim_data_2000driest= "D:/GIS_Data/REnviroLayers/mixed_data_2000_250mdriest/"
  clim_data_2100wettest="D:/GIS_Data/REnviroLayers/mixed_data_2100_250mwettest/"
  clim_data_2100driest= "D:/GIS_Data/REnviroLayers/mixed_data_2100_250mdriest/"  
}else{
  clim_data_2000="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/250m/"
  clim_data_2100="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_future/250m/"
}

##########################
####RUNNING SCRIPTS!!!####
##########################
library(stringr)
library(raster)
if (apply_biomod2_fixes){
  maxentWDtmp = paste("maxentWDtmp_", baseline_or_future, sep = "")
  dir.create(dir_for_temp_files, showWarnings=F, recursive=T)
}#dir.create(paste('Y:/temp/', project_name,'/', sep=''), showWarnings=F)
dir.create(working_dir, showWarnings=F)

#this code below will subset species into the right number of instances started with the bat file                        
Sys.sleep(6) #time for script process to show up on tasklist
n_instances=length(list.files(working_dir, pattern="^00instance"))
cpucores=6#as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
if (paralelize){
  if (cpucores>length(spp_nm)){cpucores=length(spp_nm)}
  jnkn=length(spp_nm)
  x=c(1:jnkn)
  chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
  groups=chunk(x,cpucores)
  jnk=groups[n_instances+1][[1]]
  spp_nm=spp_nm[jnk]
  spp_str=""
  for (sp_nm in spp_nm){
    spp_str=paste(spp_str,sp_nm,sep="__")
  }
  time=Sys.time()
  time=str_replace_all(time,":", ".")
  instance_file=paste0("00instance",spp_str,"_",time)
  file.create(paste0(working_dir,instance_file),showWarnings=F)  
}

coverMap=raster("Y:/PICCC_analysis/FB_analysis/habitat analysis/dominant cover/landfire_reclass_wetland_coastal_500m.tif")
hab_assoc=data.frame(cover=c(0:13))
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
  
  data=mySpeciesOcc[,3:4]
  sp_hab=extract(coverMap,data)
  jnk=as.data.frame(table(sp_hab))
  hab_assoc=cbind(hab_assoc,rep(0,dim(hab_assoc)[1]))
  names(hab_assoc)[dim(hab_assoc)[2]]=sp_nm
  hab_assoc[hab_assoc[,1] %in% jnk$sp_hab,dim(hab_assoc)[2]]=jnk$Freq
}

write.table(hab_assoc, file = "hab_assoc.csv", sep=",", col.names=NA)
rel_hab_assoc=data.frame(mapply('/', hab_assoc,colSums(hab_assoc)))
rel_hab_assoc=rel_hab_assoc*100
rel_hab_assoc[,1]=hab_assoc[,1]
write.table(rel_hab_assoc, file = "rel_hab_assoc.csv", sep=",", col.names=NA)

write.table(all_mySpeciesOcc, file = "all_pres_points.csv", sep=",", col.names=NA)
