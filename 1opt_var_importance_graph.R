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
    jpeg_name=paste0("output_rasters/varImp/",sp_nm, "_",model, "_variable_importance_violin_plot.jpg")
    a=qplot(Predictor, Value, data=long, geom=c("violin"),
            fill=Predictor, main="",
            xlab="", ylab="Variable importance" )
    a=a + geom_violin(scale = "width")
    a=a+theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
    ggsave(filename=jpeg_name, plot=a)    
    cat("done with ", sp_nm, " ",model, " variable importance box plot", '\n')
  }
}
