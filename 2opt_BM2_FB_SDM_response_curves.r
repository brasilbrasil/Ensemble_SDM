#rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
####START UNDERHOOD
setwd(working_dir)
last_model = models_to_run[length(models_to_run)]

library(biomod2)
library(stringr)
dir.create("output_rasters/", showWarnings = FALSE)
dir.create("output_rasters/response_curves/", showWarnings = FALSE)
dir.create("output_rasters/response_curves/combo/", showWarnings = FALSE)
#sp_nm = spp_nm[1] #for testing
#spp_nm=rev(spp_nm)
for (sp_nm in spp_nm){
  sp_nm = as.character(sp_nm) #convert any species names to characters 
  cat('\n',sp_nm,'modeling...')
  sp_nm0 = sp_nm
  
  workspace_name = paste0(sp_nm,"_FB_modelfitting.RData") #set name of file to save all workspace data after model run
  if (file.exists(workspace_name)){
    load(workspace_name)    
  }else{
    workspace_name = paste0(sp_nm0,"_FB_run.RData") #set name of file to load workspace data from model run    
    load(workspace_name)    
  }
  
  sp_nm = str_replace_all(sp_nm,"_", ".") #convert any species names with "_" to "."
  
  model = models_to_run[1]
  for (model in models_to_run){
    temp_jpeg_name = paste0('output_rasters/response_curves/combo/', sp_nm,"_", "response_curves","_",model,"_all_vars.jpg")
    if (file.exists(temp_jpeg_name) == F | overwriteData == T){
      
      loaded_models <- BIOMOD_LoadModels(myBiomodModelOut, models = model)
      loaded_models = loaded_models[1:length(loaded_models)-1]
      
      # 4.2 plot 2D response plots
      myRespPlot2D <- response.plot2(models  = loaded_models,
                                     Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                                     show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                                     do.bivariate = FALSE,
                                     fixed.var.metric = 'mean',
                                     save.file="no", 
                                     name="response_curve", 
                                     ImageSize=480, 
                                     plot=FALSE)
      
      ####These next sections need to be changed because "myRespPlot2D" is now a list rather than an array####
      
      bioclim_cnt = 1 #for testing
      for (bioclim_cnt in 1:length(myRespPlot2D)) { #for each bioclimatic variable
        ymax_lim = 0
        ymin_lim = 10000
        for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
          ymax_lim = max(ymax_lim, max(myRespPlot2D[[bioclim_cnt]][rep]))
          ymin_lim = min(ymin_lim, min(myRespPlot2D[[bioclim_cnt]][rep]))
        }
        xmax_lim = max(myRespPlot2D[[bioclim_cnt]][1])
        xmin_lim = min(myRespPlot2D[[bioclim_cnt]][1])
        var_name = names(myRespPlot2D)[bioclim_cnt]
        
        #intiialize the plot as a jpeg
        jpeg_name = paste0('output_rasters/response_curves/', sp_nm,"_", "response_curve", "_", model, "_", var_name, ".jpg")
        jpeg(jpeg_name,
             width = 10, height = 8, units = "in",
             pointsize = 12, quality = 90, bg = "white", res = 300)
        
        var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
        for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
          pred = rowMeans(myRespPlot2D[[bioclim_cnt]][rep])
          if (rep == 2){
            plot(var, pred, type="l", xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim), xlab = var_name, ylab = "Response", col = "grey")
          } else {
            lines(var, pred, type="l", xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim), col = "grey")
          }
        }
        #Add average response line
        var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
        pred = rowMeans(myRespPlot2D[[bioclim_cnt]][2:length(myRespPlot2D[[bioclim_cnt]])])
        lines(var, pred, type="l", xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim), lwd = 3)        
        dev.off()
      }
      
#       OLD CODE - had to be converted to deal with list
#       bioclim_cnt = 1
#       for (bioclim_cnt in 1:dim(myRespPlot2D)[3]){
#         ymax_lim=max(myRespPlot2D[,2,bioclim_cnt,])
#         ymin_lim=min(myRespPlot2D[,2,bioclim_cnt,])
#         xmax_lim=max(myRespPlot2D[,1,bioclim_cnt,])
#         xmin_lim=min(myRespPlot2D[,1,bioclim_cnt,])
#         var_name=dimnames(myRespPlot2D)[[3]][bioclim_cnt]
#         
#         jpeg_name=paste('output_rasters/response_curves/', sp_nm,"_", "response_curve","_",model,"_",var_name,".jpg", sep = "")
#         jpeg(jpeg_name,
#              width = 10, height = 8, units = "in",
#              pointsize = 12, quality = 90, bg = "white", res = 300) 
#         for (rep in 1:dim(myRespPlot2D)[4]){
#           var=myRespPlot2D[,1,bioclim_cnt,rep]
#           pred=myRespPlot2D[,2,bioclim_cnt,rep]
#           if (rep==1){
#             plot(var,pred, type="l", xlim=c(xmin_lim, xmax_lim),ylim=c(ymin_lim,ymax_lim), xlab=var_name, ylab="Response", col="grey")
#           }else{
#             lines(var,pred, type="l", xlim=c(xmin_lim, xmax_lim),ylim=c(ymin_lim,ymax_lim), col="grey")
#           }
#         }
#         #Add average response line
#         var=rowMeans(myRespPlot2D[1:dim(myRespPlot2D)[1],1,bioclim_cnt,])
#         pred=rowMeans(myRespPlot2D[1:dim(myRespPlot2D)[1],2,bioclim_cnt,])
#         lines(var,pred, type="l", xlim=c(xmin_lim, xmax_lim),ylim=c(ymin_lim,ymax_lim), lwd=3)        
#         dev.off()
#       } 
#       
      ###########combo figure########
      ###############################
      jpeg_name_combined = paste0('output_rasters/response_curves/combo/', sp_nm,"_", "response_curves","_",model,"_all_vars.jpg")
      jpeg(jpeg_name_combined,
           width = 10, height = 8, units = "in",
           pointsize = 12, quality = 90, bg = "white", res = 300) 
      
      ymax_lim = 0
      ymin_lim = 10000
      for (bioclim_cnt in 1:length(myRespPlot2D)) { #for each bioclimatic variable
        for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) { #for each of the model reps
          ymax_lim = max(ymax_lim, max(myRespPlot2D[[bioclim_cnt]][rep]))
          ymin_lim = min(ymin_lim, min(myRespPlot2D[[bioclim_cnt]][rep]))
        }
      }
      
      par(mfrow = c(2,2), oma = c(0,0,3,0))
      bioclim_cnt = 1 #for testing
      for (bioclim_cnt in 1:length(myRespPlot2D)) {
        xmax_lim = max(myRespPlot2D[[bioclim_cnt]][1])
        xmin_lim = min(myRespPlot2D[[bioclim_cnt]][1])
        
        var_name = names(myRespPlot2D)[bioclim_cnt]
        
        
        for (rep in 2:length(myRespPlot2D[[bioclim_cnt]])) {
          var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
          pred = rowMeans(myRespPlot2D[[bioclim_cnt]][rep])
          if (rep == 2){
            plot(var, pred, type = "l", xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim), xlab = var_name, ylab = "Response", col = "grey")
          } else {
            lines(var, pred, type = "l", xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim), col = "grey")
          }
        }
        
        #Add average response line
        var = rowMeans(myRespPlot2D[[bioclim_cnt]][1])
        pred = rowMeans(myRespPlot2D[[bioclim_cnt]][2:length(myRespPlot2D[[bioclim_cnt]])])
        lines(var, pred, type="l", xlim = c(xmin_lim, xmax_lim), ylim = c(ymin_lim,ymax_lim), lwd = 3)          
      } 
      
      temp_title = paste0(model, " modeled response for ", sp_nm)
      mtext(temp_title, adj = 0.5, side = 3, outer = TRUE) 
      dev.off()
      cat('\n','response curves for ',sp_nm, model, 'finished...')  
    }else{
      cat('\n','response curves for ',sp_nm, model, 'already done...')  
    }      
  }
}
