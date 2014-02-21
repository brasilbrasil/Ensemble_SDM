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
eval_stats=c("ROC","KAPPA", "TSS")

###START UNDERHOOD
setwd(working_dir)
library(biomod2)

sp_nm=spp_nm[1]
i=1
for (sp_nm in spp_nm){
  workspace_name=paste(sp_nm,"_FB_modelfitting.RData", sep = "") #set name of file to save all workspace data after model run
  load(workspace_name)
  myBiomodModelEval <- getModelsEvaluations(myBiomodModelOut)    
  dimnames(myBiomodModelEval)
  #######Loading datasets#######
  for (eval_stat in eval_stats){
    myBiomodModelEval[eval_stat,"Testing.data",,,]
    Spp_eval<- data.frame(myBiomodModelEval[eval_stat,"Testing.data",,,])
    Spp_eval=cbind(matrix(sp_nm,dim(Spp_eval)[1],1),rownames(Spp_eval),Spp_eval)
    
    if (i==1){
      assign(paste0("all_eval_mat_",eval_stat),Spp_eval)
      #all_eval_mat=Spp_eval
    }else{
      jnk=rbind(get(paste0("all_eval_mat_",eval_stat)),Spp_eval)
      assign(paste0("all_eval_mat_",eval_stat),jnk)      
    }
  }
  
  ## getting the variable importance ##
  getModelsVarImport(myBiomodModelOut)
  Spp_VariImp<- data.frame(getModelsVarImport(myBiomodModelOut))
  Spp_VariImp=cbind(matrix(sp_nm,dim(Spp_VariImp)[1],1),rownames(Spp_VariImp),Spp_VariImp)
  if (i==1){
    all_var_imp_mat=Spp_VariImp
  }else{
    all_var_imp_mat=rbind(all_var_imp_mat,Spp_VariImp)
  }
  
  i=i+1
}
FileName<-paste("all_VariImp.csv")
write.table(all_var_imp_mat, file = FileName, sep=",", row.names = FALSE)

for (eval_stat in eval_stats){
  FileName<-paste0("all_eval_mat_",eval_stat,".csv")
  write.table(get(paste0("all_eval_mat_",eval_stat)), file = FileName, sep=",", row.names = FALSE)
}