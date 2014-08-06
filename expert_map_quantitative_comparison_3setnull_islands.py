import arcpy
#FOR SOME REASON THIS ONLY WORKS WITHIN ARCPY GUI
#USER INPUT
#this code will find all rasters in a folder and apply a projection from a file
#if you want to apply for only one type of raster, change the "All" value to the extension name desired
prj_file=arcpy.Raster("Y:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_less_PAs/output_rasters/Iiwi_BIN_baseline_ROC_ef.pmw.tif") #had to define datum on arcgis D_WGS_1984
datadir="Y:/PICCC_analysis/FB_analysis/FB range expert maps/tifs only latlon comparison/aggregate_median/reprojected/"
outputdir="Y:/PICCC_analysis/FB_analysis/FB range expert maps/tifs only latlon comparison/aggregate_median/reprojected/setNull/"

#projection="C:/Program Files (x86)/ArcGIS/Desktop10.0/Coordinate Systems/Geographic Coordinate Systems/North America/NAD 1983.prj"
#END USER INPUT
import os
from random import randrange
from types import *

if not os.path.exists(outputdir):
    os.mkdir(outputdir)
jnk=randrange(10000)
arcpy.env.overwriteOutput = True
arcpy.env.workspace = datadir
arcpy.env.compression = "LZW"
arcpy.CreateFileGDB_management("D:/temp/arcgis/", "scratchoutput"+str(jnk)+".gdb")
arcpy.env.scratchWorkspace = "D:/temp/arcgis/scratchoutput"+str(jnk)+".gdb"
if arcpy.CheckExtension("Spatial") == "Available":
	arcpy.CheckOutExtension("Spatial")
#prj_file=datadir+prj_file_name
rasterList = arcpy.ListRasters("*", "TIF")
f=rasterList[0]
for f in rasterList:
    print "starting %s" %(f)
    sp=f[:-9]
    mask_layer0="Y:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_less_PAs/output_rasters/%s_BIN_baseline_ROC_ef.pmw.tif"%(sp)
    mask_layer=arcpy.sa.Con(mask_layer0,0,0,"value >= 0")
    maskName=outputdir+"temp_"+sp+"_mask.tif"
    #maskName=outputdir+"temp_mask.tif"
    arcpy.CopyRaster_management(mask_layer,maskName,"","","","","","8_BIT_SIGNED")
    rastername=os.path.join(datadir, f)
    jnk=arcpy.sa.SetNull(rastername, rastername, "Value = 0")
    spName=outputdir+"temp_"+sp+"_range.tif"
    #spName=outputdir+"temp_range.tif"
    arcpy.CopyRaster_management(jnk,spName,"","","","","","8_BIT_SIGNED")
    arcpy.MosaicToNewRaster_management("%s;%s"%(maskName,spName), outputdir, f, "",\
                                       "8_BIT_UNSIGNED", "", "1", "LAST","")
    ##arcpy.MosaicToNewRaster_management("mask_layer;jnk", outputdir, f, "",\
    ##                arcpy.CheckInExtension("Spatial")