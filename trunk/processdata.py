#! /usr/bin/env python

import sys,os
import numpy as np
from datetime import datetime
day = str(datetime.now())[0:10]

if len(sys.argv) == 1:
    print 'syntax: '
    print 'processdata.py sourceImageDir expidFile'
    print 'example: '
    print 'processdata.py /home3/data_local/images/fits/2012B-9999/ expid.txt'
else:
    sourceDIR = sys.argv[1]
    expidFile=sys.argv[2]
    b = np.genfromtxt(expidFile,dtype='i')
    for expid in b:
        os.system('cp '+sourceDIR+'DECam_00'+str(expid)+'.fits.fz'+' .')
    os.system('funpack *.fz')
    os.system('rm *.fz')
    os.system('desImgQuickReduction.py all')
    os.system('desImgQuality.py all')
    os.system('htmlFig.py '+ day)

    
