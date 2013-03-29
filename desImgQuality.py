#! /usr/bin/env python

'''
This executable module implements the analysis on each image and generate the whisker plot, the R50, whisker, whisker rms

Created on Feb 1, 2013

@author: jghao
'''

import sys, time, glob
sys.path.append('/usr/remote/user/sispi/jiangang/des-image-hao/imgQualityHao')


if __name__ == "__main__":
    from func_def import *

    startTime=time.time()
    if len(sys.argv) == 1:
        print 'syntax: '
        print 'desImgQuality.py expid'
        print 'or'
        print 'desImgQuality.py all'
        print 'Note: The image need to be reduced (bias subtraction, flat fielding'
    elif sys.argv[1] == 'all':
        img_nameList = glob.glob('*_reduced.fits')
        nimg = len(img_nameList)
        for i in range(nimg):
            t=des_img_analysis(img_nameList[i])
    else:   
        expid = sys.argv[1]
        img_name = 'DECam_00'+expid+'_reduced.fits'
        t=des_img_analysis(img_name)
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)
