rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
server=1
models_to_run=c('RF','MAXENT', 'GBM')
overwrite=0
local_config_dir='Y:/FB_analysis/FB_SDM/biomod2/' #'C:/Users/lfortini/'
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run_4.csv', sep = ""),header=F, stringsAsFactors=F))
spp_nm=c("Kauai_Amakihi")# "Akekeke", "Akikiki", "Anianiau", "Kauai_Amakihi", "Kauai_Elepaio", "Puaiohi", "Oahu_Amakihi", "Oahu_Elepaio", "Apapane", "Iiwi",)   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"
if(server==1){
  working_dir='Y:/FB_analysis/FB_SDM/biomod2/response_curves/'
}else{
  working_dir='C:/Users/lfortini/Data/biomod2/test/'
}


####START UNDERHOOD
setwd(working_dir)

library(biomod2)
library(stringr)
#sp_nm="Akepa" #debug
dir.create("output_rasters/", showWarnings = FALSE)
dir.create("output_rasters/response_curves/", showWarnings = FALSE)

#memory.limit(size=24000000)

sp_nm=spp_nm[1]
for (sp_nm in spp_nm){
  sp_nm=as.character(sp_nm)  
  cat('\n',sp_nm,'modeling...')
  sp_nm0=sp_nm
  workspace_name=paste(sp_nm0,"_FB_run.RData", sep = "") #set name of file to load workspace data from model run
  load(workspace_name)
  sp_nm=str_replace_all(sp_nm,"_", ".")
  
  model = models_to_run[1]
  for (model in models_to_run){
    loaded_models <- BIOMOD_LoadModels(myBiomodModelOut, models=model)
    loaded_models=loaded_models[1:length(loaded_models)-1]
    
    # 4.2 plot 2D response plots
    myRespPlot2D <- response.plot2(models  = loaded_models,
                                   Data = getModelsInputData(myBiomodModelOut,'expl.var'), 
                                   show.variables= getModelsInputData(myBiomodModelOut,'expl.var.names'),
                                   do.bivariate = FALSE,
                                   fixed.var.metric = 'mean',
                                   save.file="no", 
                                   name="response_curve", 
                                   ImageSize=480, 
                                   plot=FALSE)
    
    last_bioclim=dimnames(myRespPlot2D)[[3]][length(dimnames(myRespPlot2D)[[3]])]
    last_jpeg_name=paste('output_rasters/response_curves/', sp_nm,"_", "response_curve","_",model,"_",last_bioclim,".jpg", sep = "")
    if (file.exists(last_jpeg_name)==F | overwrite==1){
      #bioclim_cnt=1
      for (bioclim_cnt in 1:dim(myRespPlot2D)[3]){
        ymax_lim=max(myRespPlot2D[,2,bioclim_cnt,])
        ymin_lim=min(myRespPlot2D[,2,bioclim_cnt,])
        xmax_lim=max(myRespPlot2D[,1,bioclim_cnt,])
        xmin_lim=min(myRespPlot2D[,1,bioclim_cnt,])
        var_name=dimnames(myRespPlot2D)[[3]][bioclim_cnt]
        
        jpeg_name=paste('output_rasters/response_curves/', sp_nm,"_", "response_curve","_",model,"_",var_name,".jpg", sep = "")
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300) 
        for (rep in 1:dim(myRespPlot2D)[4]){
          var=myRespPlot2D[,1,bioclim_cnt,rep]
          pred=myRespPlot2D[,2,bioclim_cnt,rep]
          if (rep==1){
            plot(var,pred, type="l", xlim=c(xmin_lim, xmax_lim),ylim=c(ymin_lim,ymax_lim), xlab=var_name, ylab="Response", col="grey")
          }else{
            lines(var,pred, type="l", xlim=c(xmin_lim, xmax_lim),ylim=c(ymin_lim,ymax_lim), col="grey")
          }
        }
        #Add average response line
        var=rowMeans(myRespPlot2D[1:dim(myRespPlot2D)[1],1,bioclim_cnt,])
        pred=rowMeans(myRespPlot2D[1:dim(myRespPlot2D)[1],2,bioclim_cnt,])
        lines(var,pred, type="l", xlim=c(xmin_lim, xmax_lim),ylim=c(ymin_lim,ymax_lim), lwd=3)        
        dev.off()
      } 
      ##all bioclims combined
      jpeg_name=paste('output_rasters/response_curves/', sp_nm,"_", "response_curves","_",model,"all_vars.jpg", sep = "")
      jpeg(jpeg_name,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300) 
      
      par(mfrow=c(2,2))
      for (bioclim_cnt in 1:dim(myRespPlot2D)[3]){
        ymax_lim=max(myRespPlot2D[,2,bioclim_cnt,])
        ymin_lim=min(myRespPlot2D[,2,bioclim_cnt,])
        xmax_lim=max(myRespPlot2D[,1,bioclim_cnt,])
        xmin_lim=min(myRespPlot2D[,1,bioclim_cnt,])
        var_name=dimnames(myRespPlot2D)[[3]][bioclim_cnt]
        
        for (rep in 1:dim(myRespPlot2D)[4]){
          var=myRespPlot2D[,1,bioclim_cnt,rep]
          pred=myRespPlot2D[,2,bioclim_cnt,rep]
          if (rep==1){
            plot(var,pred, type="l", xlim=c(xmin_lim, xmax_lim),ylim=c(ymin_lim,ymax_lim), xlab=var_name, ylab="Response", col="grey")
          }else{
            lines(var,pred, type="l", xlim=c(xmin_lim, xmax_lim),ylim=c(ymin_lim,ymax_lim), col="grey")
          }
        }
        #Add average response line
        var=rowMeans(myRespPlot2D[1:dim(myRespPlot2D)[1],1,bioclim_cnt,])
        pred=rowMeans(myRespPlot2D[1:dim(myRespPlot2D)[1],2,bioclim_cnt,])
        lines(var,pred, type="l", xlim=c(xmin_lim, xmax_lim),ylim=c(ymin_lim,ymax_lim), lwd=3)        
      } 
      dev.off()
    }else{
      cat('\n',model, 'response graph for ',sp_nm,'already done...')
    }
    
  }    
}
