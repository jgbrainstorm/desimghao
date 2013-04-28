#! /usr/bin/env python

import sys,os
import numpy as np

sourceDIR = '/home3/data_local/images/fits/2012B-9999/'

if len(sys.argv) == 1:
    print 'syntax: '
    print 'processdata.py expidFile'
else:
    expidFile=sys.argv[1]
    b = np.genfromtxt(expidFile,dtype='i')
    for expid in b:
        os.system('cp '+sourceDIR+'DECam_00'+str(expid)+'.fits.fz'+' .')
    os.system('funpack *.fz')
    os.system('rm *.fz')
    os.system('desImgQuickReduction.py all')
    os.system('desImgQuality.py all')
    os.system('htmlFig.py 3-26/27-2013 bestimage')

    
