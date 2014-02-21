
###USER CONFIGURATION
local_config_dir='C:/Users/lfortini/'
#local_config_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/' #'C:/Users/lfortini/'
spp_nm=c('Akekee', 'Akiapolauu', 'Akikiki', 'Akohekohe', 'Anianiau', 'Hawaii_Akepa', 'Hawaii_Creeper', 'Oahu_Amakihi','Hawaii_Elepaio', 'Kauai_Elepaio', 'Maui_Alauahio', 'Maui_Parrotbill', 'Omao', 'Oahu_Elepaio', 'Palila', 'Puaiohi', 'Kauai_Amakihi', 'Hawaii_Amakihi', 'Apapane', 'Amakihi', 'Elepaio', 'Iiwi')
server=0
overwrite=1
proj_names=c("Pres", "Pres_Abs", "Isl")
Quant = 0.025
Quant_val=(1-Quant*2)*100

if (server==1){
  working_dir='Y:/FB analysis/FB SDM/biomod2/run_1_7_12_15/'
  clim_data_dir0="Y:/SDM_env_data/bioclim_variables/full extent bioclim data/all_grd/all_baseline/100m/" 
}else{
  working_dir='D:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_220runs/'
  clim_data_dir0="D:/PICCC_data/climate_data/bioclim_data_Aug2013/complete_rasters/allYrs_avg/bioclims_abs/all_baseline/250m/"
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

####GRAPH
library(plyr)
library(reshape2)

merged_data2=merged_data[merged_data$type=="Pres",]
jnk=merged_data[merged_data$type=="Pres_Abs",]
jnk2=merged_data[merged_data$type=="Isl",]
merged_data2=cbind(merged_data2, jnk[,c("QL", "QU")], jnk2[,c("QL", "QU")])
names(merged_data2)=c("X", "var", "type", "species", "P_QL", "P_QU", "PA_QL", "PA_QU", "Is_QL", "Is_QU")
sp_order=c(19, 20, 21, 22, 18, 9, 13, 15, 6, 7, 2, 4, 11, 12, 8, 14, 17, 10, 5, 1, 3, 16)

library(ggplot2)
var="bio1"
for (var in unique(merged_data2[,"var"])){
  df=merged_data2[merged_data2$var==var,]
  df=df[sp_order,]
  df$species <- factor(df$species, levels=df$species)
  q=ggplot(df) +
    geom_crossbar(aes(ymin = Is_QL, ymax = Is_QU, x = species, y = Is_QL), fill = "grey", fatten = 0)+
    geom_crossbar(aes(ymin = PA_QL, ymax = PA_QU, x = species, y = PA_QL), fill = "red", fatten = 0)+
    geom_crossbar(aes(ymin = P_QL, ymax = P_QU, x = species, y = P_QL), fill = "blue", fatten = 0)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  q + ylab(var)
  jpeg_name=paste0("quantiles/Q", Quant_val, "_", var, "_all_spp.jpg")
  ggsave(filename=jpeg_name, plot=q)
}



ggplot(df) +
  geom_crossbar(aes(ymin = min, ymax = max, x = id, y = min),
                fill = "blue", fatten = 0)

df2=df
df2$min2=df2[,2]-1
df2$max2=df2[,3]+1
ggplot(df2) +
  geom_crossbar(aes(ymin = min2, ymax = max2, x = id, y = min), fill = "red", fatten = 0)+
  geom_crossbar(aes(ymin = min, ymax = max, x = id, y = min), fill = "blue", fatten = 0)



####SAVE TABLE
FileName=paste("quantiles/Q", Quant_val, "_p_vs_pa_all_spp.csv", sep="")
write.table(merged_data, file = FileName, sep=",", col.names=NA)



