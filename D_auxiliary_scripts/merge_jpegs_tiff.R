rm(list = ls()) #remove all past worksheet variables
wd="D:/Dropbox/current work/FB SDM/FB MS draft/not in MS/base figs/140807 version/hiRel/arcgis figs/to_merge/"
setwd(wd)
#library(jpeg)
library(tiff)
#library(rasterImage)
files=list.files(wd)
strs=c("Fig 2", "Fig 3", "Fig 4", "Fig 6")

str = strs[1]
for (str in strs){
  jnk=grep(paste0(str,"+"), files, perl=TRUE, value=FALSE)
  selected_files=files[jnk]
  #selected_files=grep(paste0(str,"+"), files, perl=TRUE, value=FALSE)
  selected_files=selected_files[grep("merged",selected_files,invert=TRUE)]
  
  image1 <- readTIFF(selected_files[1], native=T)
  image2 <- readTIFF(selected_files[2], native=T)
  
  jpeg_name=paste(str, "_merged.jpg", sep = "")
  jpeg(jpeg_name,
       width = 14, height = 8, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  par(mai=c(0.1,0.1,0.1,0.1))
  plot(0:1, 0:1, type='n', xaxt='n', yaxt='n', axes=F, ann=F)
  rasterImage(image1, 0, 0, 0.48, 1)
  rasterImage(image2, 0.5, 0, 1, 1)
  dev.off()

  tif_name=paste(str, "_merged.tif", sep = "")
  tiff(tif_name,
       width = 10, height = 6, units = "in",compression = "lzw",
       pointsize = 12, bg = "white", res = 300)
  par(mai=c(0.1,0.1,0.1,0.1))
  plot(0:1, 0:1, type='n', xaxt='n', yaxt='n', axes=F, ann=F)
  rasterImage(image1, 0, 0, 0.48, 1)
  rasterImage(image2, 0.5, 0, 1, 1)
  text(x=0.03, y = 0.85, labels = "(a)", adj = NULL, #x=0.08, y = 0.97
       pos = NULL, offset = 0.5, vfont = NULL,
       cex = 1.8, col = "white", font = 2)
  text(x=0.53, y = 0.85, labels = "(b)", adj = NULL,
       pos = NULL, offset = 0.5, vfont = NULL,
       cex = 1.8, col = "white", font = 2)
  
  dev.off()
  
}