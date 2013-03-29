#!/usr/bin/env python
# build a catalog using sextractor on des image

import sys,glob
sys.path.append('/usr/remote/user/sispi/jiangang/decam-fermi')
#from DECamCCD import *

if len(sys.argv) == 1:
    print 'syntax: '
    print 'getStarCat image.fits'
    print 'The image need to be reduced'
else:
   
   def sex(image, output, sexdir='/usr/remote/user/sispi/jiangang/des-image-hao/imgQualityHao/sexconfig/', check_img=None,config=None):
      '''Construct a sextractor command and run it.'''
      if sexdir[-1] != '/':  sexdir += '/'
      com = ["sex", image, "-c "+sexdir+config,"-CATALOG_NAME "+output,"-CHECKIMAGE_NAME "+check_img]
      com = string.join(com)
      res = os.system(com)
      return res

   img_name = sys.argv[1]
   output=img_name[0:-5]+'_star_catalog.fits'
   ckimg=img_name[0:-5]+'check.fits'
   t=sex(img_name,output,check_img=ckimg,config="initial.sex")

print '----done !---'


