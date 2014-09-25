#must have island map open: open Island extent raster or dem raster
#must have option to show temp output on map DISABLED;
#if multiple frames, click on properties of each frame and fix extent so redrawing does not screw with map layout
#disable background processing (arcgis geoprocessing options) as it does not work well- python will call variables for routines that are not yet complete from previous routine!!
#load common datasets, create legend in mxd (make a copy in case something goes wrong!)
#if template has a legend, after legend is done, rclk legend and turn off all map connection options, keep layers in at least one
#data frame, but unchecked
#for rasters with categorical (unique vals) ref layer must keep non-display values as no color (i.e., do not remove value)
#if categorical raster, make sure raster is not a float raster, create reference layer from same set of files

#USER INPUT
print_map_output=1
overwrite_res=1
df_to_exclude=["inset"]

#START UNDERHOOD
import arcpy, os, string, logging, datetime
from arcpy import env
import csv
import time
arcpy.env.overwriteOutput = True
#arcpy.env.workspace = resultsDir
mxd=arcpy.mapping.MapDocument("CURRENT")
arcpy.env.compression = "LZW"
df = arcpy.mapping.ListDataFrames(mxd)[0]
#refLayer = arcpy.mapping.ListLayers(mxd, "*", df)[0]

def del_layer(layer_name):
    for df in arcpy.mapping.ListDataFrames(mxd):
        for lyr in arcpy.mapping.ListLayers(mxd, "", df):
            if lyr.name.lower() == layer_name:
                arcpy.mapping.RemoveLayer(df, lyr)
    return

def print_map(To_load, ref_layers, names_to_load, overwrite_res, output_path):
    added_layers=[]
    if overwrite_res==0 and arcpy.Exists(output_path):
        print output_path + " calculated already"
    else:
        df = arcpy.mapping.ListDataFrames(mxd)[0]
        i=0
        for df in arcpy.mapping.ListDataFrames(mxd):
            if df.name not in df_to_exclude:
                mxd.activeView = df.name
                for i in range(len(To_load)):
                    fc=To_load[i]
                    sourceLayer0=arcpy.mapping.Layer(ref_layers[i])
                    if fc[-4:]==".tif":
                        result=arcpy.MakeRasterLayer_management(fc, "lyr"+str(i)) #os.path.basename(fc)
                        layer = result.getOutput(0)
                        if df.name==arcpy.mapping.ListDataFrames(mxd)[0].name:
                            added_layers.append("lyr"+str(i))
                    else:
                        layer = arcpy.mapping.Layer(fc)
                        if df.name==arcpy.mapping.ListDataFrames(mxd)[0].name:
                            added_layers.append(os.path.basename(fc)[:-4])
                    #arcpy.ApplySymbologyFromLayer_management(layer, sourceLayer0)
                    arcpy.mapping.UpdateLayer(df, layer, sourceLayer0)
                    arcpy.mapping.AddLayer(df, layer, "TOP")
                arcpy.RefreshActiveView()
                arcpy.RefreshTOC()
        if print_map_output==1:
            arcpy.mapping.ExportToTIFF(map_document=mxd, out_tiff=output_path,resolution=600,tiff_compression="LZW")
            for df in arcpy.mapping.ListDataFrames(mxd):
                #mxd.activeView = df.name
                for lyr in arcpy.mapping.ListLayers(mxd, "", df):
                    if lyr.name in added_layers:
                        arcpy.mapping.RemoveLayer(df, lyr)

##                for i in range(len(To_load)):
##                    fc=To_load[i]
##                    if fc[-4:]==".tif":
##                        arcpy.mapping.RemoveLayer(df, "lyr"+str(i))
##                    else:
##                        arcpy.mapping.RemoveLayer(df, os.path.basename(fc)[:-4])
##                    jnks=arcpy.mapping.ListBrokenDataSources(df)
##                    for jnk in jnks:
##                        arcpy.mapping.RemoveLayer(df, jnk)

            arcpy.RefreshActiveView()
            arcpy.RefreshTOC()
            #remove these map layers!!
    return


#iterating data sets
#shpdir=r"Y:/PICCC_analysis/FB_analysis/habitat_analysis/veg_overlay/current_veg_mask/" #whichever is data dir, will have to have subfolders: results/, results/la/, la/ (where you place CCE and FCE files)
tifdir=r"D:/Dropbox/current work/FB SDM/FB MS draft/not in MS/base figs/140807 version/hiRel/"
refdir=tifdir+"ref_layers_and_mxd/"
resultsDir=tifdir+"arcgis_figs_script/"
import glob
import os
os.chdir(tifdir)
tifs=glob.glob("*.tif")

for tif in tifs:
    #sp="Apapane"
    To_load=[tifdir+tif] #locs
    ref_layers=[refdir+"sp_ensemble_ref_layer.lyr"] #locs
    names_to_load=["response_zones"] #strings for referencing
    output_path=resultsDir+tif[:-4]+"_leg.tif"
    print_map(To_load, ref_layers, names_to_load, overwrite_res, output_path)
