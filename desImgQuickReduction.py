#! /usr/bin/env python

"""
This code perform only the overscan subtraction for each image extension, no masterbias, master flat are used. 
J. Hao @ FNAL, 11/28/2012
"""

import numpy as np
import pyfits as pf
import sys,time
import glob as gl


def oscanSub(img):
    """
    Note that this code only subtract the overscan of individual 4k x 2k imaging ccd from DECam. 
    The L and R is only for the fits image. This is the most recent one as of 9/30/2012 J. Hao
    """
    oscanL = img[:,10:50]
    oscanR = img[:,2110:2150]
    mdL=np.median(oscanL,axis=1)
    mdR=np.median(oscanR,axis=1)
    #rowL=np.arange(0,mdL.shape[0])
    #rowR=np.arange(0,mdR.shape[0])
    #(aL,bL,sda,sdb,se)=linefit(rowL,mdL)
    #(aR,bR,sda,sdb,se)=linefit(rowR,mdR)
    #oscanLfit=rowL*bL+aL
    #oscanRfit=rowR*bR+aR
    for i in range(1080):
        img[:,i] = img[:,i] - mdL #oscanLfit
        img[:,1080+i] = img[:,1080+i] - mdR #oscanRfit
    return img

if len(sys.argv) == 1:
    print 'syntax: '
    print '   desImgQuickReduction.py epxid  '
    print 'or:'
    print '   desImgQuickReduction.py all'
    print '   which process all the files in the current directory'
    print 'example:' 
    print '   desImgQuickReduction.py 121414'
    print '   desImgQuickReduction.py all'

else:
    startTime=time.time()
    if sys.argv[1] == 'all':
        filename = gl.glob('*.fits')
        nimg = len(filename)
    else:
        filehead = 'DECam_00'
        nimg=len(sys.argv) - 1
    for i in range(nimg):
        if sys.argv[1] == 'all':
            hdu = pf.open(filename[i],mode='update')
            correctedFilename = filename[i][0:-5]+'_reduced.fits'
        else:
            hdu = pf.open(filehead+sys.argv[1+i]+'.fits',mode='update')
            correctedFilename = filehead+'_'+sys.argv[1+i]+'_reduced.fits'
        hdu[0].header.update('PROCTYPE','Reduced')
        for ext in range(1,71):
            print ext
            hdu[ext].data = oscanSub(hdu[ext].data)
        hdu.writeto(correctedFilename)
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)

    
