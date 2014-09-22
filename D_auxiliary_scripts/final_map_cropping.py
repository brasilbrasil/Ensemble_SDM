import os
#import Pillow
#import Image
import glob
from PIL import Image
#http://stackoverflow.com/questions/6997419/how-to-create-a-loop-to-read-several-images-in-a-python-script
wd=r"D:/Dropbox/current work/FB SDM/FB MS draft/not in MS/base figs/140807 version/all/arcgis figs/"
os.chdir(wd)
files=glob.glob("*.tif")
suffix="_cropped"
filename=files[0]
for filename in files:
    if not suffix in filename:
        im=Image.open(filename)
        w, h = im.size
        box=(850, 280, w-280, h-1450) #0,0 at top left # left, top, right, down #2550, 3300
        im_crop=im.crop(box)
        #im_crop.show() to display
        outname=filename[:-4]+suffix+".tif"
        im_crop.save(outname)
