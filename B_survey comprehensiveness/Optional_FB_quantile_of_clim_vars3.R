rm(list = ls()) #remove all past worksheet variables

###USER CONFIGURATION
source(paste0("C:/Users/lfortini/","directory_registry.r"))
spp_nm=c('Akekee', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi', 'Kauai_Amakihi', 'Hawaii_Amakihi', 'Apapane', 'Iiwi') #'Amakihi', 'Elepaio', 
server=1
overwrite=1
exclude_abs=F
exclude_pres=F
island_quantile=F

Quant = 0.025

Quant_val=(1-Quant*2)*100

if (!exclude_abs & !exclude_pres & !island_quantile){
  proj_name="Pres_Abs"  
}
if (island_quantile){
  proj_name="Isl"  
}
if (exclude_abs & !exclude_pres & !island_quantile){
  proj_name="Pres"  
}
if (!exclude_abs & exclude_pres & !island_quantile){
  proj_name="Abs"  
}


if (server==1){
#   working_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/'
#   clim_data_dir0="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/100m/" 
  working_dir='Y:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_less_PAs/'
  clim_data_dir0="Y:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg/bioclims_abs/all_baseline/250m/"
}else{
  working_dir='D:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_220runs/'
  clim_data_dir0="D:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg/bioclims_abs/all_baseline/250m/"
}

env_var_files=c("bio1.tif", "bio7.tif", "bio12.tif", "bio15.tif") 
#env_var_fileslabels=c("MAT (c)", "Annual T range (c)", "MAP (mm)", "P seasonality") 
env_var_fileslabels=c(expression("MAT ("*~degree*"C)"), expression("Annual T range ("*~degree*"C)"), "MAP (mm)", "P seasonality") 


crop_raster_dir=paste(working_dir, 'map_crop/',sep="")
csv_dir=paste(working_dir,"single_sp_CSVs/", sep="")

###START UNDERHOOD
setwd(working_dir)
library(biomod2)
library(raster)
library(randomForest)
library(dismo)
library(mda)
library(stringr)

all_quantiles=c()
spp_info=read.csv(paste(csv_dir,'FB_spp_data.csv', sep = ""))
dir.create("quantiles/", showWarnings=FALSE)

var_name=c()
for (env_var_file  in env_var_files){
  a=strsplit(env_var_file,"\\.")
  var_name=c(var_name, a[[1]][1])
}

sp_nm = spp_nm[1]
for (sp_nm in spp_nm){  
  sp_nm=as.character(sp_nm)
  cat('\n',sp_nm,'modeling...')
  sp_str=sp_nm
  
  FileName<-paste("quantiles/Q", Quant_val, "_", proj_name, "_", sp_nm,".csv", sep="")
  if (file.exists(FileName)==F | overwrite==1){
    #######Loading datasets#######
    sp_data=read.csv(paste(csv_dir,sp_nm,'_pres_abs.csv', sep = "")) #FB_data_points4_PAandA
    
    if (!island_quantile){ #this step is skipped if island quantile since species points are not used
      #so processing is not necessary
      # Select, Count and Remove presence Duplicates
      jnk=dim(sp_data)[1]
      dups2<- duplicated(sp_data[, c('X','Y')])
      sum(dups2)
      sp_data<-sp_data[!dups2, ]
      jnk1=dim(sp_data)[1]
      if (exclude_abs){
        sp_data=sp_data[sp_data[,"pa"]==1,] #get rid of absences      
      }
      if (exclude_pres){
        sp_data=sp_data[sp_data[,"pa"]==0,] #get rid of absences            
      }
      jnk2=dim(sp_data)[1]
      head(sp_data)
      cat('\n','removed ', jnk-jnk1, "duplicates for", sp_nm)
      cat('\n','removed ', jnk1-jnk2, "absence records for", sp_nm)      
    }
    
    ##raster_based_env_grid:
    sp_index=which(spp_info[,"Species"]==sp_nm)
    raster_res= spp_info[sp_index,"rasterdir"]
    clim_data_dir=clim_data_dir0 
    jnk0=length(env_var_files)
    crop_raster=raster(paste(crop_raster_dir,raster_res,".grd",sep=""))
    predictors = raster( paste(clim_data_dir, env_var_files[1], sep=""))
    predictors=crop(predictors,  crop_raster)
    for (jj in 2:jnk0){
      temp=raster(paste(clim_data_dir, env_var_files[jj], sep=""))
      if (island_quantile){
        okHab=raster(okrasterDir, "okHabRaster.tif")
        temp=temp*okHab
      }
      temp=crop(temp,  crop_raster)
      predictors = addLayer(predictors, temp)
    }
    names(predictors)<- var_name
    predictors
    
    cat('\n',sp_nm,'modeling...')
    sp_data=data.frame(cbind(x=sp_data$X, y=sp_data$Y, pres=sp_data$pa))
    

    n_points=dim(sp_data)[1]
    Response_var=sp_data
    if (island_quantile){
      layer=1
      for (layer in 1:nlayers(predictors)){
        jnk=c(as.matrix(predictors[[layer]]))
        if (layer==1){
          explanatory=as.matrix(jnk)
        }else{
          explanatory=cbind(explanatory, jnk)
        }
      }
      explanatory=explanatory[!is.na(explanatory[,1]),]
      explanatory=as.data.frame(explanatory)
      names(explanatory)=var_name
    }else{
      explanatory=extract(predictors, Response_var[,1:2])      
    }
    if (proj_name=="Pres_Abs"){
      sp_data2=cbind(sp_data,as.data.frame(explanatory))
      var_name_i = var_name[1]
      for (var_name_i in var_name){        
        data1=sp_data2[,c("pres", var_name_i)]
        data1[data1$pres==0,"pres"]="Absences"
        data1[data1$pres==1,"pres"]="Presences"
        data1$pres=as.factor(data1$pres)
        data1$pres=factor(data1$pres, c("Presences", "Absences"))
        
        env_var_fileslabel=env_var_fileslabels[which(var_name_i == var_name)] 
        
        #Compound histogram
        bw=range(data1[,var_name_i])
        bw=(bw[2]-bw[1])/10
        require(ggplot2)
        library(grid) #for unit argument
        a=qplot(get(var_name_i), data = data1, geom = "histogram" , binwidth=bw, group = pres, fill = pres) + theme_bw( )
        a=a+xlab(env_var_fileslabel)+ylab("Count")+ theme(legend.title=element_blank()) 
        a=a+theme(legend.position=c(.15, .85)) #where x,y is 0,0 (bottom left) to 1,1 (top right)
        a=a+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) #remove grids
        a=a+theme(panel.grid=element_blank(), #remove outward ticks
                  axis.ticks.length = unit(-0.25,'cm'),
                  axis.ticks.margin = unit(0.5,'cm'))        
        #a

        plotname=paste("quantiles/Histogram_", sp_nm,"_",var_name_i,".tiff", sep="")
        ggsave(filename = plotname, plot = a,  
               scale = 1, width = par("din")[1]*0.6, height = par("din")[2]*0.6, 
               units = "in", dpi = 1200, compression = "lzw")#save image c("in", "cm", "mm")

        #proportional presence
        xl=c(floor(min(sp_data2[,var_name_i])), ceiling(max(sp_data2[,var_name_i])))
        bw=range(data1[,var_name_i])
        brks=seq(from = bw[1], to = bw[2], by = (bw[2]-bw[1])/20)
        xl=c(floor(min(sp_data2[,var_name_i])), ceiling(max(sp_data2[,var_name_i])))
        p1 <- hist(sp_data2[sp_data2$pres==0,var_name_i], breaks=brks)
        p2 <- hist(sp_data2[sp_data2$pres==1,var_name_i], breaks=brks)
        
        abs_counts=p1$counts
        pres_counts=p2$counts
        mids=p1$mids
        pres_prop=pres_counts/(pres_counts+abs_counts)
        data2=data.frame(mids, pct_pres=pres_prop)
        
        b=ggplot(data2, aes(x = mids, y=pres_prop)) 
        b=b+geom_line(colour="grey", linetype="dotted", size=1.5)
        b=b+geom_point(colour="black", size=4, shape=21, fill="white")
        #b
        #b=ggplot(data2) + geom_bar(aes(x = mids, y=pres_prop),stat="identity") #http://stackoverflow.com/questions/19389888/ggplot-barplot-does-not-accept-y-value
        #b
        b=b+ ylim(0, 1)
        b=b+theme_bw( )
        b=b+xlab(env_var_fileslabel)+ylab("Proportional presence") 
        b=b+ theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) #remove grids
        b=b+theme(panel.grid=element_blank(), #remove outward ticks
                  axis.ticks.length = unit(-0.15,'cm'),
                  axis.ticks.margin = unit(0.5,'cm'))        
        #b
        
        plotname=paste("quantiles/Proportional_presence_", sp_nm,"_",var_name_i,".tiff", sep="")
        ggsave(filename = plotname, plot = b,  
               scale = 1, width = par("din")[1]*0.4, height = par("din")[2], 
               units = "in", dpi = 1200, compression = "lzw")#save image c("in", "cm", "mm") 
      }
    }
    
    extrem.cond <- t(apply(as.data.frame(explanatory), 2, quantile, probs = c(0 + Quant, 1 - Quant), 
                           na.rm = TRUE))          
    extrem.cond=as.data.frame(extrem.cond)
    extrem.cond=cbind(species=matrix(sp_nm,dim(extrem.cond)[1],1),extrem.cond)
    extrem.cond=cbind(type=matrix(proj_name,dim(extrem.cond)[1],1),extrem.cond)
    extrem.cond=cbind(var=var_name,extrem.cond)
    row.names(extrem.cond) <- NULL 
    write.table(extrem.cond, file = FileName, sep=",", col.names=NA)
    
  }else{
    cat('\n', sp_nm, " already calculated")
  }
  if (sp_nm==spp_nm[1]){
    all_quantiles=extrem.cond
  }else{
    all_quantiles=rbind(all_quantiles,extrem.cond)    
  }
}
names(all_quantiles)[4:5]=c("QL","QU")
FileName<-paste("quantiles/Q", Quant_val, "_", proj_name, "_all_spp.csv", sep="")
write.table((all_quantiles), file = FileName, sep=",", col.names=NA)
