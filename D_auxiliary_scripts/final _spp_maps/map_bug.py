#Have a blank arcmap open
import arcpy, os
from arcpy import env
arcpy.env.overwriteOutput = True
mxd=arcpy.mapping.MapDocument("CURRENT")
df = arcpy.mapping.ListDataFrames(mxd)[0]
#specify data location
datadir=r"Y:/temp/map_bug/" #whichever is data dir, will have to have subfolders: results/, results/la/, la/ (where you place CCE and FCE files)
sp="Iiwi"
To_load=[datadir+"%s_response_zones_ROC_ef.pmw_int.tif"%(sp)] #locs
ref_layers=[datadir+"FB_response_zones4.lyr"] #locs

df = arcpy.mapping.ListDataFrames(mxd)[0]
i=0 #index is here because purpose of original code is to loop operation below
mxd.activeView = df.name
fc=To_load[i]
sourceLayer0=arcpy.mapping.Layer(ref_layers[i])
result=arcpy.MakeRasterLayer_management(fc, "lyr"+str(i)) #os.path.basename(fc)
layer = result.getOutput(0)
arcpy.ApplySymbologyFromLayer_management(layer, sourceLayer0)
arcpy.mapping.AddLayer(df, layer, "TOP")

#this is what does not work
arcpy.ApplySymbologyFromLayer_management("lyr"+str(i), sourceLayer0)
arcpy.mapping.UpdateLayer(df, layer, sourceLayer0)


