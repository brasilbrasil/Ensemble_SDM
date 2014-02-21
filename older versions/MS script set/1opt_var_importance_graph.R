rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))
#for each species modeled, have csv of presence data in working directory for the species named speciesname_Ps.csv formated with 3 cols: x,y,pa where pa = 1
#after running the code for whichever many species, copy results (species output folder and workspace file) to a new directory, along with the maxent.jar file
#use the projection code to project the distribution model on different environmental surfaces (do not forget to change the working directory)

###USER CONFIGURATION
local_config_dir='Y:/FB_analysis/FB_SDM/biomod2/' #if specifiying sp to run by file, this is directory of where csv file is located
spp_nm = c('Akekee', 'Hawaii_Amakihi', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi','Kauai_Amakihi', 'Iiwi', 'Amakihi', 'Elepaio', 'Apapane')
#'Iiwi', 'Amakihi', 'Elepaio', 'Apapane', 
project_name='finalmodel_P_PA_oldcode_less_PAs'
working_dir=paste0(resultsDir,project_name,'/')
models_to_run=c('GBM','MAXENT')

###START UNDERHOOD
setwd(working_dir)

varImp0=read.csv('all_VariImp.csv')
colnames=names(varImp0)
model_cols=list()
for (model in models_to_run){
  jnk=grep(paste0(model,"+"), colnames, perl=TRUE, value=FALSE)  
  model_cols[[length(model_cols)+1]]=jnk
}

sp_nm=spp_nm[1]
dir.create("output_rasters/varImp/", showWarnings=F)
for (sp_nm in spp_nm){
  for (model in models_to_run){
    model_n=which(models_to_run==model)
    varImp=varImp0[varImp0[,1]==sp_nm,model_cols[[model_n]]]
    varImp=as.data.frame(t(varImp))
    names(varImp)=varImp0[varImp0[,1]==sp_nm,"rownames.Spp_VariImp."]
    
    library(ggplot2)
    library(reshape2)
    library(plyr)
    long=melt(varImp)
    names(long)=c("Predictor", "Value")
    jpeg_name=paste0("output_rasters/varImp/",sp_nm, "_",model, "_variable_importance_box_plot.jpg")
    a=qplot(Predictor, Value, data=long, geom=c("boxplot", "jitter"), 
            fill=Predictor, main="",
            xlab="", ylab="Variable importance" )
    ggsave(filename=jpeg_name, plot=a)    
    cat("done with ", sp_nm, " ",model, " variable importance box plot", '\n')
  }
}
