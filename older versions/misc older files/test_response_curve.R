# 0. Load data & Selecting Data

# species occurances
species_occ <- read.csv(system.file("external/species/species_occ.csv",package="biomod2"))

# we consider only presences of MyocastorCoypus species
myRespName <- 'MyocastorCoypus'
myRespCoord <- species_occ[which(!is.na(species_occ[,myRespName])),c('x','y')]
myResp <- as.numeric(na.omit(species_occ[,myRespName]))

# Environemental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack(system.file("external/climat/current/bio3.grd",package="biomod2"),
                       system.file("external/climat/current/bio4.grd",package="biomod2"),
                       system.file("external/climat/current/bio7.grd",package="biomod2"),
                       system.file("external/climat/current/bio11.grd",package="biomod2"),
                       system.file("external/climat/current/bio12.grd",package="biomod2"))

# 1. Formating Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 200,
                                     PA.strategy = 'random')

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('RF','GLM'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=1, 
                                       DataSplit=80, 
                                       Yweights=NULL, 
                                       VarImport=3, 
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE )


# 4. Plot response curves

# 4.1 Load the models for which we want to extract the predicted response curves
myGLMs <- BIOMOD_LoadModels(myBiomomodModelOut, models='GLM')

# 4.2 plot 2D response plots
myRespPlot2D <- response.plot2(models  = myGLMs,
                               Data = getModelsInputData(myBiomomodModelOut,'expl.var'), 
                               show.variables= getModelsInputData(myBiomomodModelOut,'expl.var.names'),
                               do.bivariate = FALSE,
                               fixed.var.metric = 'mean',
                               save.file="no", 
                               name="response_curve", 
                               ImageSize=480, 
                               plot=TRUE)

# 4.2 plot 3D response plots
## here only for a lone model (i.e "MyocastorCoypus_PA1_RUN1_GLM")
myRespPlot3D <- response.plot2(models  = myGLMs[1],
                               Data = getModelsInputData(myBiomomodModelOut,'expl.var'), 
                               show.variables= getModelsInputData(myBiomomodModelOut,'expl.var.names'),
                               do.bivariate = TRUE,
                               fixed.var.metric = 'mean',
                               save.file="no", 
                               name="response_curve", 
                               ImageSize=480, 
                               plot=TRUE)

### all the values used to produce this plot are stored into the returned object
### you can redo plots by yourself and customised them
dim(myRespPlot2D)
dimnames(myRespPlot2D)

dim(myRespPlot3D)
dimnames(myRespPlot3D)