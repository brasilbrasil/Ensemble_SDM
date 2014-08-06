import arcpy
#USER INPUT
#this code will find all rasters in a folder and apply a projection from a file
#if you want to apply for only one type of raster, change the "All" value to the extension name desired
prj_file=arcpy.Raster("Y:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_less_PAs/output_rasters/Iiwi_BIN_baseline_ROC_ef.pmw.tif") #had to define datum on arcgis D_WGS_1984
datadir="Y:/PICCC_analysis/FB_analysis/FB range expert maps/tifs only latlon comparison/aggregate_median/"
outputdir="Y:/PICCC_analysis/FB_analysis/FB range expert maps/tifs only latlon comparison/aggregate_median/reprojected/"

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

#prj_file=datadir+prj_file_name
rasterList = arcpy.ListRasters("*", "TIF")
f=rasterList[0]
for f in rasterList:
    print "starting %s" %(f)
    sp=f[:-9]
    prj_file="Y:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_less_PAs/output_rasters/%s_BIN_baseline_ROC_ef.pmw.tif"%(sp)
    #arcpy.env.extent=prj_file
    rastername=os.path.join(datadir, f)
    outrastername=os.path.join(outputdir, f)
    arcpy.ProjectRaster_management(rastername, outrastername, prj_file, "NEAREST", "#", "#", "#", "#")
