#USER INPUT
#del_terms=["protected_area.csv", "CO_points_w_zones", "tabulate_points_temp", "core_biome_", "fragmentation_", "GBU_", "intersect_", "Migrate_", "refuge", "simplified_", "ungfree_map", "temp_vals", "protected", "COR_CCE", "COR_FCE"]
del_term="veg_mask"
#rootdir=r"C:/Users/lfortini/toFWS/"
rootdir=r"Y:/PICCC_analysis/FB_analysis/habitat_analysis/veg_overlay/current_veg_mask/" #whichever is data dir, will have to have subfolders: results/, results/la/, la/ (where you place CCE and FCE files)
do_tif_conversion=True
do_pol_creation=True

#END USER INPUT
import os
import arcpy
from types import *
arcpy.env.overwriteOutput = True
arcpy.env.workspace = rootdir
#mxd=arcpy.mapping.MapDocument("CURRENT")
#arcpy.env.compression = "LZW"


for root, dirs, files in os.walk(rootdir):
    for f in files:
        if f[-4:]==".tif":
            if del_term in f:
                if "_int.tif" not in f:
                    #os.unlink(os.path.join(root, f))
                    print "found " + f + " in " + root
                    f1_lyr=root+f[:-4] +"_int.tif"
                    fc_lyr=root+f[:-4] +".shp"
                    #fc_lyr=f[:-4] +".shp"
                    print "saving as " + fc_lyr
                    if do_tif_conversion:
                        arcpy.CopyRaster_management(f,f1_lyr,"DEFAULTS","","0","","","8_BIT_UNSIGNED") #as integet!
                    if do_pol_creation:
                        arcpy.RasterToPolygon_conversion(f1_lyr, fc_lyr, "NO_SIMPLIFY", "")

