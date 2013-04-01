#! /usr/bin/env python
import glob as gl
import sys,os

if len(sys.argv) == 1:
    print 'syntax: htmlFig.py date MJD'
    print 'example: htmlFig.py 11-17-2012 56248'
    sys.exit()
date = sys.argv[1]
mjd = sys.argv[2]
htmlName = 'Image_Quality_'+date+'-MJD-'+mjd+'.html'
Fig_coeff = gl.glob('zernike_coeff*.png')
Fig_moments = gl.glob('moments*.png') 

Fig_coeff.sort()
Fig_moments.sort()

nfig = len(Fig_coeff)

htm=open(htmlName,'w')
htm.write('<HTML> \n')

htm.write('<HEAD> \n')
htm.write('<TITLE>Image Quality of '+date+'</TITLE>\n')
htm.write('</HEAD> \n')
htm.write('<BODY> \n')
htm.write('<p>Image Quality of images taken during: '+date+'</p>\n')
htm.write('<p>MJD: '+mjd+'</p>\n')
htm.write('<p>J. Hao @ FNAL</p>\n')


for i in range(nfig):
    htm.write('<p>\n')
    htm.write('<img src="%s" width="1000">\n'%Fig_moments[i])
    htm.write('<img src="%s" width="1000">\n'%Fig_coeff[i])
    htm.write('</p>\n')
htm.write('</BODY> \n')
htm.write('</HTML> \n')
htm.close()
os.system('tar -czf '+htmlName+'.tar.gz *.png *.html')
