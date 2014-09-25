#USER INPUT
#search_terms=["protected_area.csv", "CO_points_w_zones", "tabulate_points_temp", "core_biome_", "fragmentation_", "GBU_", "intersect_", "Migrate_", "refuge", "simplified_", "ungfree_map", "temp_vals", "protected", "COR_CCE", "COR_FCE"]
search_term="pmw"
#rootdir=r"C:/Users/lfortini/toFWS/"
rootdir=r"Y:/PICCC_analysis/FB_analysis/model_results/biomod2finalmodel_P_PA_oldcode_less_PAs/output_rasters/main/"
do_tif_conversion=True

#END USER INPUT
import os
import arcpy
from types import *
arcpy.env.overwriteOutput = True
arcpy.env.workspace = rootdir

for root, dirs, files in os.walk(rootdir):
    for f in files:
        if f[-4:]==".tif":
            if search_term in f:
                if "_int.tif" not in f:
                    print "found " + f + " in " + root
                    f1_lyr=root+f[:-4] +"_int.tif"
                    if do_tif_conversion:
                        arcpy.CopyRaster_management(f,f1_lyr,"DEFAULTS","","0","","","8_BIT_UNSIGNED") #as integet!

