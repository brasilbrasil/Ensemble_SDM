rm(list = ls()) #remove all past worksheet variables
#for each species modeled, have csv of presence data in working directory for the species named speciesname_Ps.csv formated with 3 cols: x,y,pa where pa = 1
#after running the code for whichever many species, copy results (species output folder and workspace file) to a new directory, along with the maxent.jar file
#use the projection code to project the distribution model on different environmental surfaces (do not forget to change the working directory)

###USER CONFIGURATION
local_config_dir='C:/Users/lfortini/' #if specifiying sp to run by file, this is directory of where csv file is located
spp_nm=(read.csv(paste(local_config_dir,'spp_to_run.csv', sep = ""),header=F, stringsAsFactors=F))
#spp_nm=c("Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi")   #"Akekee", "Akikiki", "Anianiau", "Apapane", "Iiwi", "Kauai_Amakihi", "Kauai_Elepaio", "Oahu_Amakihi", "Oahu_Elepaio", "Puaiohi"

working_dir='Y:/FB SDM/biomod2/'

###START UNDERHOOD
setwd(working_dir)
library(biomod2)

#sp_nm=spp_nm[1]
i=1
for (sp_nm in spp_nm){
  workspace_name=paste(sp_nm,"_FB_run.RData", sep = "") #set name of file to save all workspace data after model run
  load(workspace_name)
  #######Loading datasets#######
  
#   # Outputting the validation metrics for all tests
#   myBiomodModelEval["TSS","Testing.data",,,]
#   Spp_TSS<- data.frame(myBiomodModelEval["TSS","Testing.data",,,])
#   if (i==1){
#     all_eval_mat=Spp_TSS
#   }else{
#     all_eval_mat=rbind(all_eval_mat,Spp_TSS)
#   }
  
  #FileName<-paste(sp_nm, "_TSS.csv")
  #write.table(Spp_TSS, file = FileName, sep=",", col.names=NA)
  
  myBiomodModelEval["ROC","Testing.data",,,]
  Spp_ROC<- data.frame(myBiomodModelEval["ROC","Testing.data",,,])
  all_eval_mat=rbind(all_eval_mat,Spp_ROC)
  #FileName<-paste(sp_nm, "_ROC.csv")
  #write.table(Spp_ROC, file = FileName, sep=",", col.names=NA)
  
#   myBiomodModelEval["KAPPA","Testing.data",,,]
#   Spp_KAP<- data.frame(myBiomodModelEval["KAPPA","Testing.data",,,])
#   all_eval_mat=rbind(all_eval_mat,Spp_KAP)
#   #FileName<-paste(sp_nm, "_KAP.csv")
#   #write.table(Spp_KAP, file = FileName, sep=",", col.names=NA)
  
  ## getting the variable importance ##
  getModelsVarImport(myBiomodModelOut)
  Spp_VariImp<- data.frame(getModelsVarImport(myBiomodModelOut))
  if (i==1){
    all_var_imp_mat=Spp_VariImp
  }else{
    all_var_imp_mat=rbind(all_var_imp_mat,Spp_VariImp)
  }
  
  i=i+1
}
FileName<-paste("all_VariImp.csv")
write.table(all_var_imp_mat, file = FileName, sep=",", col.names=NA)

FileName<-paste("all_eval_mat.csv")
write.table(all_eval_mat, file = FileName, sep=",", col.names=NA)
