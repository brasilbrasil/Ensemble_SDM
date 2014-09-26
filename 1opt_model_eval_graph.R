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
    a=a+theme(legend.position="none")
    ggsave(filename=jpeg_name, plot=a)    
    cat("done with ", eval_stat, " ",model, " variable importance box plot", '\n')
  }
}

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
    jpeg_name=paste0("output_rasters/evalMat/",eval_stat, "_",model, "_variable_importance_violin_plot.jpg")
    a=qplot(Species, Value, data=long, geom=c("violin"), 
            fill=Species, main="",
            xlab="", ylab=paste(eval_stat, "model evaluation"))#+ theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    a=a + geom_violin(scale = "width")
    a=a+guides(fill = guide_legend(keywidth = 1, keyheight = 1.2))
    a=a+theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
    a=a+coord_flip()
    a=a+guides(fill = guide_legend(reverse=TRUE), guide = guide_legend(title = NULL))
    a
    ggsave(filename=jpeg_name, plot=a)    
    cat("done with ", eval_stat, " ",model, " variable importance box plot", '\n')
  }
}
