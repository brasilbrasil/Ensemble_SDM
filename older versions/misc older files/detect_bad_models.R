myBiomodModelEval[,,"MAXENT", "RUN1", "PA9"]

jnk=myBiomodModelEval[,2,,,]
NAs=which(is.na(jnk))
all_models=list()
for (d in dimnames(jnk)[[4]]){
  for (c in dimnames(jnk)[[3]]){
    for (b in dimnames(jnk)[[2]]){
      for (a in dimnames(jnk)[[1]]){
        jnk_str=paste(sp_nm,d,c,b,a,sep="_")
        jnk_str2=paste(sp_nm,d,c,b,sep="_")
        all_models[length(all_models)+1]=jnk_str2
      }}}}
bad_models_short=all_models[NAs]
bad_models_short=unique(bad_models_short)
jnk_bad=!(all_models %in% bad_models_short)
remaining_models=all_models[jnk_bad]
remaining_models=unique(remaining_models)

grep('_MAXENT', get_built_models(myBiomodModelOut),
     value=TRUE)