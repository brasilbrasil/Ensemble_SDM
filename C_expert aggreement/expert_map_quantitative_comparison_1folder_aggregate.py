##this code is not finished trying to do it by shapefile means
#USER INPUT
rootdir=r"Y:/PICCC_analysis/FB_analysis/FB range expert maps/tifs only latlon comparison/" #whichever is data dir, will have to have subfolders: results/, results/la/, la/ (where you place CCE and FCE files)
results_dir=rootdir+"aggregate_median/"
overwrite=True
#END USER INPUT


import os
import arcpy
from random import randrange
from types import *

jnk=randrange(10000)
arcpy.env.overwriteOutput = True
arcpy.env.workspace = rootdir
arcpy.env.compression = "LZW"
arcpy.CreateFileGDB_management("D:/temp/arcgis/", "scratchoutput"+str(jnk)+".gdb")
arcpy.env.scratchWorkspace = "D:/temp/arcgis/scratchoutput"+str(jnk)+".gdb"

if not os.path.exists(results_dir):
    os.mkdir(results_dir)

if arcpy.CheckExtension("Spatial") == "Available":
	arcpy.CheckOutExtension("Spatial")

rasterList = arcpy.ListRasters("*", "tif")

mask_layer=rootdir+"consensus_full_class_12.tif" #ignore pixel values from this mask in aggregation

cell_factor=500/30
f=rasterList[2]
for f in rasterList:
    jnk=f[:-4]+"_aggr.tif"
    out_name="%s%s" %(results_dir,jnk)
    if arcpy.Exists(out_name) and overwrite==False:
		print "raster " + out_name+" already done"
    else:
        #maskedRaster = arcpy.sa.SetNull(mask_layer, f, "Value = 100")
        maskedRaster = arcpy.Raster(f)
        maskedRaster=arcpy.sa.Con(arcpy.sa.IsNull(maskedRaster),0,maskedRaster)
        outRas=arcpy.sa.Aggregate(maskedRaster, cell_factor, "MEDIAN", "TRUNCATE", "DATA")
        arcpy.CopyRaster_management(outRas,out_name,"","","","","","8_BIT_SIGNED")

arcpy.CheckInExtension("Spatial")
