# Defineing the extent of the different islands each Trupanea lives on 
Kauai = c(-159.82,-159.26, 21.84, 22.25)
Oahu = c(-158.32, -157.62,  21.22, 21.73)
Molokai = c(-157.34, -156.69, 21.03, 21.25)
Lanai = c(-157.08, -156.78, 20.70, 20.92)
Maui= c(-156.8, -155.53, 20.46, 21.05)
Hawaii = c(-156.10,-154.74, 18.87, 20.30)
Kahoolawe = c(-156.8, -156.51, 20.46, 20.62)

# Cutting out each island
for (i in 1:length(Islands)){
  e = extent(get(Islands[i]))
  Isras = crop(predstack, e, snap = 'in')
  assign(Islands[i], Isras)
}

# for Maui I need to cut out Kahoolawe, but because of the extent issues one has to first reclass, merge and then reclass the merged Kahoo to NA
rcl <- c(0, 10000, -1)
rcl <- matrix(rcl, ncol=3, byrow=TRUE)
a=reclassify(Kahoolawe, rcl)  # defining the Kahoolawe raster values as -1 so that they can be distinguished once merged with Maui
Maui = merge(Maui, a)

rcl <- c(-1, NA)
rcl <- matrix(rcl, ncol=2, byrow=TRUE)
Maui=reclassify(Maui, rcl) 
# if an error occurs in this section it's usually associated with the spp_nm files

sp_row<- which(sp_islanddata[,"Species"]==sp_nm)
spIsland<- sp_islanddata[sp_row,(2:length(names(sp_islanddata)))]
jnk<-seq(1:length(names(sp_islanddata))-1)*spIsland
jnkx<-jnk[jnk>0]
spIsland2<- names(sp_islanddata)[c(jnkx+1)]

if ((length(spIsland2)>1)==T){
  pred1<-get(spIsland2[1])
  for (isle in 2:length(spIsland2)){
    pred1<-merge(pred1, get(spIsland2[isle]))     
  }
  predstack<-pred1
}else{
  predstack<-get(spIsland2)
}

#   plot(pred1, useRaster=F)
#   predstack<-merge(get(spIsland2[1:length(spIsland2)]))
names(predstack)<-var_name