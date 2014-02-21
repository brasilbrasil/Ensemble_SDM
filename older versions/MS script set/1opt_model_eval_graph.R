rm(list = ls()) #remove all past worksheet variables
source(paste0("C:/Users/lfortini/","directory_registry.r"))
#for each species modeled, have csv of presence data in working directory for the species named speciesname_Ps.csv formated with 3 cols: x,y,pa where pa = 1
#after running the code for whichever many species, copy results (species output folder and workspace file) to a new directory, along with the maxent.jar file
#use the projection code to project the distribution model on different environmental surfaces (do not forget to change the working directory)

###USER CONFIGURATION
local_config_dir='Y:/FB_analysis/FB_SDM/biomod2/' #if specifiying sp to run by file, this is directory of where csv file is located
spp_nm = c('Akekee', 'Hawaii_Amakihi', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi','Kauai_Amakihi', 'Iiwi', 'Apapane')
#'Iiwi', 'Amakihi', 'Elepaio', 'Apapane', 
project_name='finalmodel_P_PA_oldcode_less_PAs'
working_dir=paste0(resultsDir,project_name,'/')
models_to_run=c('GBM','MAXENT')
eval_stats=c("ROC","KAPPA", "TSS")

###START UNDERHOOD
setwd(working_dir)


dir.create("output_rasters/evalMat/", showWarnings=F)
eval_stat = eval_stats[1]
model = models_to_run[1]
for (eval_stat in eval_stats){
  evalMat0=read.csv(paste0('all_eval_mat_',eval_stat,'.csv'))
  for (model in models_to_run){
    evalMat=evalMat0[evalMat0[,2]==model,3:ncol(evalMat0)]
    evalMat=as.data.frame(t(evalMat))
    names(evalMat)=evalMat0[evalMat0[,2]==model,1]
    
    library(ggplot2)
    library(reshape2)
    library(plyr)
    long=melt(evalMat)
    names(long)=c("Species", "Value")
    long=long[long[,1] %in% spp_nm,]
    long=long[order(long[,1]),]
    long$Species <- factor(long$Species, levels = sort(levels(long$Species)))
    jpeg_name=paste0("output_rasters/evalMat/",eval_stat, "_",model, "_variable_importance_box_plot.jpg")
    a=qplot(Species, Value, data=long, geom=c("boxplot"), 
            fill=Species, main="",
            xlab="", ylab=paste(eval_stat, "model evaluation"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    ggsave(filename=jpeg_name, plot=a)    
    cat("done with ", eval_stat, " ",model, " variable importance box plot", '\n')
  }
}
