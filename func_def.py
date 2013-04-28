'''
This module define the functions used for the image quality measurements. 
This is the final version as of 2/1/2013

Created on Feb 1, 2013

@author: jghao
'''
import sys,os
sys.path.append('/usr/remote/user/sispi/jiangang/des-sv')

try:
    import numpy as np, glob as gl, pyfits as pf, scipy.ndimage as nd, pylab as pl, cPickle as p
    from scipy.misc import factorial as fac
    from scipy.optimize import leastsq
    from DECamCCD_def import *
except ImportError:
    print "Error: missing one of the libraries (numpy, pyfits, scipy, pylab)"
    sys.exit()

scale = 0.264

def robust_mean(x):
    '''
    calculate the mean after outlier removal
    '''
    y = x.flatten()
    if len(y) < 6:
        return -999.
    if len(np.unique(y))==1:
        meany=y[0]
    else:
        n = len(y)
        y.sort()
        ind_qt1 = round((n+1)/4.)
        ind_qt3 = round((n+1)*3/4.)
        IQR = y[ind_qt3]- y[ind_qt1]
        lowFense = y[ind_qt1] - 1.*IQR
        highFense = y[ind_qt3] + 1.*IQR
        if lowFense == highFense:
            meany=lowFense
        else:
            ok = (y>lowFense)*(y<highFense)
            yy=y[ok]
            meany=yy.mean(dtype='double')
    return meany

def image_quality_summary():
    """
    make a summary plot for the r50, whk, whkrms
    """
    f = gl.glob('desIQ_measure*.txt')
    f.sort()
    b=[]
    for ff in f:
        b.append(np.genfromtxt(ff))
    b=np.array(b)
    ind = np.arange(len(f))
    expid = b[:,0].astype('i')
    whk = b[:,1]  # these are all use weighted moments after correcting weight
    phi = b[:,2]  # these are all use weighted moments after correcting weight
    whkrms = b[:,3] # these are all use weighted moments after correcting weight
    r50Sex = b[:,4] # these are from my run of sextractor
    pl.figure(figsize=(15,12))
    pl.subplot(3,1,1)
    pl.plot(ind,whk,'bo')
    pl.hlines(0.2,-0.5,len(expid),color='red')
    pl.ylim(0,0.4)
    pl.xticks(ind,np.repeat('',len(expid)))
    pl.grid()
    pl.ylabel('Whisker')
    pl.title('Weighted moments after correcting for the weights')
    pl.subplot(3,1,2)
    pl.plot(ind,whkrms,'bo')
    pl.hlines(0.2,-0.5,len(expid),color='red')
    pl.grid()
    pl.xticks(ind,np.repeat('',len(expid)))
    pl.ylim(0.,0.4)
    pl.ylabel('Whisker RMS')
    pl.subplot(3,1,3)
    pl.plot(ind,r50Sex,'bo')
    pl.hlines(0.552,-0.5,len(expid),color='red')
    pl.ylim(0,1)
    pl.grid()
    pl.ylabel('R50')
    pl.xticks(ind,expid.astype('S10'),rotation=90)
    pl.savefig('desIQ_summary.png')
    pl.close()
    np.savetxt('desIQ_summary.txt',b,fmt=['%i','%10.5f','%10.5f','%10.5f','%10.5f'])


def robust_std(x):
    '''
    calculate the standard deviation after outlier removal
    '''
    y = x.flatten()
    if len(y) < 6:
        return -999.
    if len(np.unique(y))==1:
        meany=y[0]
    else:
        n = len(y)
        y.sort()
        ind_qt1 = round((n+1)/4.)
        ind_qt3 = round((n+1)*3/4.)
        IQR = y[ind_qt3]- y[ind_qt1]
        lowFense = y[ind_qt1] - 1.*IQR
        highFense = y[ind_qt3] + 1.*IQR
        if lowFense == highFense:
            stdy=lowFense
        else:
            ok = (y>lowFense)*(y<highFense)
            yy=y[ok]
            stdy=yy.std(dtype='double')
    return stdy


def getStamp(data=None,xcoord=None,ycoord=None,Npix = None):
    """
    Input: CCD image in maxtrix, x, y centroid of stars,the stamp npix
    Output: a list of stamp image around the x,y centroid
    """
    Nstar = len(xcoord)
    rowcen = ycoord
    colcen = xcoord
    stampImg=[]
    for i in range(Nstar):
        Img = data[int(rowcen[i]-Npix/2):int(rowcen[i]+Npix/2),int(colcen[i]-Npix/2):int(colcen[i]+Npix/2)]
        if np.sum(Img) > 0:
            stampImg.append(Img)
    return stampImg


def wr(x=None,y=None,xcen=None,ycen=None,sigma=None):
    """
    Returns a circular Gaussian weights with the given parameters
    """
    res=np.exp(-((x-xcen)**2+(y-ycen)**2)/(2.*sigma**2))/(2.*np.pi*sigma**2) 
    return res


def adaptiveCentroid(data=None,sigma=None):
    """
    calculate the centroids using the adaptive approach with a kernel weights function. 
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        IWcol = IWmat.sum(axis=0)
        IWrow = IWmat.sum(axis=1)
        drowmean = np.sum((rowgrid-rowmean)*IWrow)/np.sum(IWrow)
        dcolmean = np.sum((colgrid-colmean)*IWcol)/np.sum(IWcol)
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if drowmean**2+dcolmean**2 <= EP:
            break
    return rowmean,colmean

def complexMoments(data=None,sigma=None):
    """
    This one calcualte the 3 2nd moments and 4 thrid moments with the Gaussian weights.
    col : x direction
    row : y direction
    the centroid is using the adpative centroid.
    sigma is the stand deviation of the measurement kernel in pixel
    The output is in pixel coordinate
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        IWcol = IWmat.sum(axis=0)
        IWrow = IWmat.sum(axis=1)
        IWsum = IWmat.sum()          
        drowmean = np.sum((rowgrid-rowmean)*IWrow)/IWsum
        dcolmean = np.sum((colgrid-colmean)*IWcol)/IWsum
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if drowmean**2+dcolmean**2 <= EP:
            break
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mr = np.sum(rowgrid*IWrow)/IWsum
    Mc = np.sum(colgrid*IWcol)/IWsum
    Mrr = np.sum(rowgrid**2*IWrow)/IWsum
    Mcc = np.sum(colgrid**2*IWcol)/IWsum
    Mrc = np.sum(np.outer(rowgrid,colgrid)*IWmat)/IWsum
    Mrrr = np.sum(rowgrid**3*IWrow)/IWsum
    Mccc = np.sum(colgrid**3*IWcol)/IWsum
    Mrrc = np.sum(np.outer(rowgrid**2,colgrid)*IWmat)/IWsum
    Mrcc = np.sum(np.outer(rowgrid,colgrid**2)*IWmat)/IWsum
    #print Mrrr, Mccc, Mrrc, Mrcc
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    M31 = complex(3*Mc - (Mccc+Mrrc)/sigma**2, 3*Mr - (Mrcc + Mrrr)/sigma**2)
    M33 = complex(Mccc-3*Mrrc, 3.*Mrcc - Mrrr)
    return M20, M22, M31, M33

def complex2ndMoments(data=None,sigma=None):
    """
    This one calcualte the 2nd moments with the Gaussian weights and then subract the weight contribution away
    col : x direction
    row : y direction
    the centroid is using the adpative centroid.
    sigma is the stand deviation of the measurement kernel in pixel
    The output is in pixel coordinate
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))
    maxItr = 50
    EP = 0.0001
    for i in range(maxItr):
        wrmat = wr(ROW,COL,rowmean,colmean,sigma)
        IWmat = data*wrmat
        IWcol = IWmat.sum(axis=0)
        IWrow = IWmat.sum(axis=1)
        IWsum = IWmat.sum()
        drowmean = np.sum((rowgrid-rowmean)*IWrow)/IWsum
        dcolmean = np.sum((colgrid-colmean)*IWcol)/IWsum
        rowmean = rowmean+2.*drowmean
        colmean = colmean+2.*dcolmean
        if IWsum == 0:
            return -999., -999., -999., -999., -999., -999.
        if drowmean**2+dcolmean**2 <= EP:
            break
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mr = np.sum(rowgrid*IWrow)/IWsum
    Mc = np.sum(colgrid*IWcol)/IWsum
    Mrr = np.sum(rowgrid**2*IWrow)/IWsum
    Mcc = np.sum(colgrid**2*IWcol)/IWsum
    Mrc = np.sum(np.outer(rowgrid,colgrid)*IWmat)/IWsum
    Cm = np.matrix([[Mcc,Mrc],[Mrc,Mrr]])
    Cw = np.matrix([[sigma**2,0.],[0.,sigma**2]])
    Cimg = (Cm.I - Cw.I).I
    Mcc_correct = Cimg[0,0]
    Mrr_correct = Cimg[1,1]
    Mrc_correct = Cimg[0,1]
    #M20 = Mrr + Mcc
    #M22 = complex(Mcc - Mrr,2*Mrc)
    return Mcc, Mrr, Mrc, Mcc_correct, Mrr_correct, Mrc_correct

def firstcutStar(b):
    '''
    select stars based on the mag - radius relation of the des firstcut catalog
    input is the sextractor catalog
    '''
    rad = b.FLUX_RADIUS
    mag = b.MAG_AUTO
    flag = b.FLAGS
    ok = (mag>=10.5)*(mag<=12)*(flag ==0)*(rad<5.)
    radmedian = np.median(rad[ok])
    idx = (mag>=10.5)*(mag<=13)*(flag ==0)*(abs(rad-radmedian)<=0.2)
    return b[idx]


def measureIQstamp(stamp=None,bkg=None,sigma=None):
    '''
    measure the weighted moments
    '''
    Nobj = len(stamp)
    Mrr = np.zeros(Nobj)
    Mcc = np.zeros(Nobj)
    Mrc = np.zeros(Nobj)
    Mrr_correct = np.zeros(Nobj)
    Mcc_correct = np.zeros(Nobj)
    Mrc_correct = np.zeros(Nobj)
    for i in range(Nobj):
        if bkg == None:
            Mcc[i],Mrr[i],Mrc[i],Mcc_correct[i], Mrr_correct[i],Mrc_correct[i]=complex2ndMoments(data=stamp[i],sigma=sigma)
        else:
            data = stamp[i]-bkg[i]
            if data.sum > 0.:
                Mcc[i],Mrr[i],Mrc[i],Mcc_correct[i], Mrr_correct[i],Mrc_correct[i]=complex2ndMoments(data=data,sigma=sigma)      
    return robust_mean(Mcc), robust_mean(Mrr), robust_mean(Mrc),robust_mean(Mcc_correct), robust_mean(Mrr_correct), robust_mean(Mrc_correct)


def zernike_rad(m, n, rho):
    """
    Calculate the radial component of Zernike polynomial (m, n) 
    given a grid of radial coordinates rho.
    """
    if (n < 0 or m < 0 or abs(m) > n):
        raise ValueError
    if ((n-m) % 2):
        return rho*0.0
    pre_fac = lambda k: (-1.0)**k * fac(n-k) / ( fac(k) * fac( (n+m)/2.0 - k ) * fac( (n-m)/2.0 - k ) )
    return sum(pre_fac(k) * rho**(n-2.0*k) for k in xrange((n-m)/2+1))

def zernike(m, n, rho, phi):
    """
    Calculate Zernike polynomial (m, n) given a grid of radial
    coordinates rho and azimuthal coordinates phi.
    """
    if (m > 0): return zernike_rad(m, n, rho) * np.cos(m * phi)
    if (m < 0): return zernike_rad(-m, n, rho) * np.sin(-m * phi)
    return zernike_rad(0, n, rho)


def zernikel(j, rho, phi):
    """
    Calculate Zernike polynomial with Noll coordinate j given a grid of radial coordinates rho and azimuthal coordinates phi.
    """
    n = 0
    while (j > n):
        n += 1
        j -= n
    m = -n+2*j
    return zernike(m, n, rho, phi)


def zernikeFit(x, y, z,max_rad=225.,cm=[0,0],max_order=20):
    """
    Fit a set of x, y, z data to a zernike polynomial with the least square fitting. Note that here x, y, z are all 1 dim array. Here the max_rad is by default equal to 225 mm, the size of the decam focal plane.
    It will return the beta and the adjusted R2
    """
    x = x - cm[0]
    y = y - cm[1]
    n = len(x)
    p = max_order
    rho = np.sqrt(x**2+y**2)/max_rad #normalize to unit circle.
    phi = np.arctan2(y,x)
    dataX = []
    ok = rho <= 1.
    for j in range(max_order):
        dataX.append(zernikel(j,rho[ok],phi[ok]))
    dataX=np.array(dataX).T
    beta,SSE,rank,sing = np.linalg.lstsq(dataX,z[ok])# SSE is the residual sum square
    sigma = np.sqrt(SSE/(n-p))
    betaErr = sigma/np.dot(dataX.T,dataX).diagonal()
    SST = np.var(z[ok])*(len(z[ok])-1)# SST is the sum((z_i - mean(z))^2)
    R2 = 1 - SSE/SST
    R2adj = 1-(1-R2)*(len(z[ok])-1)/(len(z[ok])-max_order)# adjusted R2 for quality of fit.             
    fitted = np.dot(dataX,beta) # fitted value
    return beta,betaErr, R2adj,fitted

def gaussian2d(x,y,xc,yc,sigmax,sigmay,rho,A,B):
    """
    2D Gaussian profile with a constant
    """
    res = A*np.exp(-0.5/(1-rho**2)*(x**2/sigmax**2+y**2/sigmay**2-2.*rho*x*y/(sigmax*sigmay)))+B
    return res

def moments(data):
    """
    Returns (height, and width)
    the gaussian parameters of a 2D distribution by calculating its
    moments
    """
    nrow,ncol=data.shape
    Isum = data.sum()
    Icol = data.sum(axis=0) # sum over all rows
    Irow = data.sum(axis=1) # sum over all cols
    colgrid = np.arange(ncol)
    rowgrid = np.arange(nrow)
    rowmean=np.sum(rowgrid*Irow)/Isum
    colmean=np.sum(colgrid*Icol)/Isum
    ROW,COL=np.indices((nrow,ncol))   
    rowgrid = rowgrid - rowmean # centered
    colgrid = colgrid - colmean
    Mrr = np.sum(rowgrid**2*data)/Isum
    Mcc = np.sum(colgrid**2*data)/Isum
    Mrc = np.sum(np.outer(rowgrid,colgrid)*data)/Isum
    width = np.sqrt((Mrr+Mcc)*0.5)
    datasrt = np.sort(data.flatten())
    height = np.median(datasrt[-20:-1])
    return height,width



def AcomplexMoments(img,sigma=1.1/scale):
    """
    calculate M20, M22 using the adaptive moments.
    x is col, y is row
    sigmax -> sigmac, sigmay -> sigmar
    """
    npix = img.shape[0]
    rowCen,colCen = adaptiveCentroid(img,sigma)
    row,col = np.mgrid[0:npix,0:npix]
    row = row - rowCen
    col = col - colCen
    A0,sigmac0 = moments(img)
    sigmar0 = sigmac0
    rho0 = 0.
    B0 = 0.
    p0=np.array([sigmac0,sigmar0,rho0,A0, B0])
    def residualg2d(p,x,y,xc,yc,I):
        sigmax,sigmay,rho,A,B = p
        Ierr = np.sqrt(abs(I))+0.00001 # to avoid those = 0, add a small number 
        res = (gaussian2d(x,y,xc,yc,sigmax,sigmay,rho,A,B) - I)/Ierr
        return res.flatten()
    p = leastsq(residualg2d,p0,args=(col,row,colCen,rowCen,img))[0]
    sigmac,sigmar,rho,A,B = p
    Mcc = sigmac**2
    Mrr = sigmar**2
    Mrc = rho**2*Mcc*Mrr
    M20 = Mrr + Mcc
    M22 = complex(Mcc - Mrr,2*Mrc)
    return M20, M22


def measure_stamp_moments(stamp,bkg=None,sigma=2.,adaptive=False):
    """
    measure the moments of stamp image on one CCD
    return the median moments
    """
    Nobj = len(stamp)
    M20=np.zeros(Nobj)
    M22=np.zeros(Nobj).astype(complex)
    M31=np.zeros(Nobj).astype(complex)
    M33=np.zeros(Nobj).astype(complex)
    for i in range(Nobj):
        if bkg == None:
            M20[i],M22[i],M31[i],M33[i]=complexMoments(data=stamp[i],sigma=sigma)
        else:
            data = stamp[i]-bkg[i]
            if data.sum > 0.:
                if adaptive == False:
                    M20[i],M22[i],M31[i],M33[i]=complexMoments(data=data,sigma=sigma)               
                else:
                    M20[i],M22[i] = AcomplexMoments(img=data,sigma = sigma)
    return [np.median(M20), np.median(M22), np.median(M31), np.median(M33)]


def measure_stamp_coeff(data = None, zernike_max_order=20):
    """
    the convention of data is: x, y, M20, M22, M31, M33
    """
    betaAll=[]
    betaErrAll=[]
    R2adjAll=[]
    beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,2].real,max_order=zernike_max_order)
    betaAll.append(beta)
    betaErrAll.append(betaErr)
    R2adjAll.append(R2_adj)
    for i in range(3,6):
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].real,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
        beta,betaErr,R2_adj = zernikeFit(data[:,0].real,data[:,1].real,data[:,i].imag,max_order=zernike_max_order)
        betaAll.append(beta)
        betaErrAll.append(betaErr)
        R2adjAll.append(R2_adj)
    betaAll = np.array(betaAll)
    betaErrAll = np.array(betaErrAll)
    R2adjAll = np.array(R2adjAll)
    return betaAll,betaErrAll, R2adjAll

def display_2nd_moments(data=None):
    # remove the mean for all moments
    datamean = data.mean(axis = 0)
    data = subMeanAll(data)
    pl.figure(figsize=(11,5.5))
    pl.subplot(1,2,1)
    phi22 = 0.5*np.arctan2(data[:,3].imag,data[:,3].real)
    x = data[:,0].real
    y = data[:,1].real
    #phi22[x<0] = phi22+np.deg2rad(180)
    u = np.abs(data[:,3])*np.cos(phi22)
    v = np.abs(data[:,3])*np.sin(phi22)
    qvr = pl.quiver(x,y,u,v,width = 0.004, color='r',pivot='middle',headwidth=0.,headlength=0.,headaxislength=0.,scale_units='width')
    #qk = pl.quiverkey(qvr, -150,-240,np.max(np.sqrt(u**2+v**2)),str(round(np.max(np.sqrt(u**2+v**2)),3))+' pix^2',coordinates='data',color='blue')
    qk = pl.quiverkey(qvr, -150,-240,0.3,str(0.3)+' pix^2',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('Camera WEST [mm]')
    pl.ylabel('Camera NORTH [mm]')
    pl.title('mean |M22|= '+str(round(abs(datamean[:,3]),5))+' pix^2')
    pl.subplot(1,2,2)
    m20sqr = np.sqrt(data[:,2].real)
    x = data[:,0].real
    y = data[:,1].real
    m20sqr_med = np.median(m20sqr) 
    m20sqr_diff = m20sqr - m20sqr_med
    m20sqr_diff_absmed = np.median(np.abs(m20sqr_diff))
    plotScale = 1./m20sqr_diff_absmed*100
    pos = m20sqr_diff >=0
    neg = m20sqr_diff < 0
    pl.scatter(x[pos],y[pos],s=m20sqr_diff[pos]*plotScale,c='r',alpha=0.5)
    pl.scatter(x[neg],y[neg],s=-m20sqr_diff[neg]*plotScale,c='b',alpha=0.5)
    pl.scatter(-230,-210,s=0.01*plotScale,c='b',alpha=0.5)
    pl.text(-200,-215,'-'+str(0.01)+' pix')
    pl.scatter(-230,-230,s=0.01*plotScale,c='r',alpha=0.5)
    pl.text(-200,-235,str(0.01)+' pix')
    pl.plot(x,y,'y,')
    pl.grid(color='g')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.xlabel('Camera WEST [mm]')
    pl.ylabel('Camera NORTH [mm]')
    pl.title('median '+r'$\sqrt{M20}$: '+str(round(scale*m20sqr_med,3))+' [arcsec]')
    return '---done!--'


def subMeanAll(data=None):
    """
    this subtract the mean of all moments except M20 from the data
    """
    datamean = data.mean(axis = 0)
    data[:,3:] = data[:,3:] - datamean[3:]
    return data

def averageN30(data=None):
    """
    this function use the mean moments of N29,N31 to replace the moments of N30 because of the failure of N30. 
    """
    idxN29 = 59
    idxN30 = 60
    idxN31 = 61
    datanew = data.copy()
    datanew[idxN30,2:] = 0.5*(datanew[idxN29,2:]+datanew[idxN31,2:])
    return datanew

def selectStarFwhm(catname):
    ext = [1,2,3,4]
    fwhm_sex=np.zeros(0)
    mag = np.zeros(0)
    for i in ext:
        cat=pf.getdata(catname,i)
        fwhm_sex=np.concatenate((fwhm_sex,cat.FWHM_IMAGE))
        mag = np.concatenate((mag,cat.MAG_AUTO))
    ok = (mag > -15)*(mag<-13)*(fwhm_sex > 0)*(fwhm_sex < 6.)
    md = np.median(fwhm_sex[ok])
    return md

def whisker4QReduce(X2WIN_IMAGE=None,Y2WIN_IMAGE=None,XYWIN_IMAGE=None):
    """
    This function make the whisker plot based on the sextractor export from QuickReduce
    J.Hao, 12/4/2012
    """
    xpos = np.genfromtxt('/usr/remote/user/sispi/jiangang/desimghao/xpos_ypos_fp.txt').T[0]
    ypos = np.genfromtxt('/usr/remote/user/sispi/jiangang/desimghao/xpos_ypos_fp.txt').T[1]
    temp = np.zeros(62).astype('complex')
    temp.real = X2WIN_IMAGE - Y2WIN_IMAGE
    temp.imag = 2*XYWIN_IMAGE
    data=np.array([xpos, ypos,X2WIN_IMAGE + Y2WIN_IMAGE,temp]).T
    data = averageN30(data)
    data = subMeanAll(data)
    pl.figure(figsize=(11,5.5))
    pl.subplot(1,2,1)
    phi22 = 0.5*np.arctan2(data[:,3].imag,data[:,3].real)
    x = data[:,0].real
    y = data[:,1].real
    u = np.abs(data[:,3])*np.cos(phi22)
    v = np.abs(data[:,3])*np.sin(phi22)
    qvr = pl.quiver(x,y,u,v,width = 0.004, color='r',pivot='middle',headwidth=0.,headlength=0.,headaxislength=0.,scale_units='width')
    qk = pl.quiverkey(qvr, -150,-240,0.3,str(0.3)+' pix^2',coordinates='data',color='blue')
    pl.plot(x,y,'b,')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.grid(color='g')
    pl.xlabel('Camera WEST [mm]')
    pl.ylabel('Camera NORTH [mm]')
    pl.title('M22')
    pl.subplot(1,2,2)
    m20sqr = np.sqrt(data[:,2].real)
    x = data[:,0].real
    y = data[:,1].real
    m20sqr_med = np.median(m20sqr) 
    m20sqr_diff = m20sqr - m20sqr_med
    m20sqr_diff_absmed = np.median(np.abs(m20sqr_diff))
    plotScale = 1./m20sqr_diff_absmed*100
    pos = m20sqr_diff >=0
    neg = m20sqr_diff < 0
    pl.scatter(x[pos],y[pos],s=m20sqr_diff[pos]*plotScale,c='r',alpha=0.5)
    pl.scatter(x[neg],y[neg],s=-m20sqr_diff[neg]*plotScale,c='b',alpha=0.5)
    pl.scatter(-230,-210,s=0.01*plotScale,c='b',alpha=0.5)
    pl.text(-200,-215,'-'+str(0.01)+' pix')
    pl.scatter(-230,-230,s=0.01*plotScale,c='r',alpha=0.5)
    pl.text(-200,-235,str(0.01)+' pix')
    pl.plot(x,y,'y,')
    pl.grid(color='g')
    pl.xlim(-250,250)
    pl.ylim(-250,250)
    pl.xlabel('Camera WEST [mm]')
    pl.ylabel('Camera NORTH [mm]')
    pl.title('median '+r'$\sqrt{M20}$: '+str(round(scale*m20sqr_med,3))+' [arcsec]')
    return '---done!--'


def selectStar(mag,fwhm_sex):
    ok = (mag > -15)*(mag<-13)*(fwhm_sex > 0)*(fwhm_sex < 10.)
    md = np.median(fwhm_sex[ok])
    return md



def hexapodPosition(beta,betaErr,weighted=True):
    """
    the CRAY position to the hexapod position parameters. There is a 15 deg rotation between the two coordinate. However, this is accounted in the sispi. So, the hexapod position in the header is acutally rotated to its designed coordiante, which relate to the CRAY coordinate by the last page of des-docdb #6551
    """
    x,y,z,thetax,thetay = CRAYposLinearModel(beta,betaErr,weighted)
    #xh = x                
    #yh = -y
    yh = x       # change on 11/30/12 after comparing with BCAM and hexapod data on 11/22/12
    xh = -y   
    zh = -z
    thetaxh = -thetay
    thetayh = -thetax
    return np.array([xh,yh,zh,thetaxh,thetayh])

def CRAYposLinearModel(b=None,bErr=None,weighted=True):
    """
    here, the convention for b is: M20 (0-19), M22real(20 - 39),M22imag (40-59) 
    This function return the hexapod parameters based on the zernike coefficients from second moments on the FP. 
    Note that the zero order of the zernike for each moments are removed in the regression to remove the seeing, tracking, guding effects.
    """
    M22realTrefoil2 = b[29] 
    M22imagTrefoil1 = b[48] 
    M22realTrefoil2Err = bErr[29] 
    M22imagTrefoil1Err = bErr[48] 
    if weighted == False:
        M22TrefoilXshift = 0.5*(M22realTrefoil2+M22imagTrefoil1) # for x decenter
    else:
        M22TrefoilXshift = (M22realTrefoil2*M22imagTrefoil1Err**2+M22imagTrefoil1*M22realTrefoil2Err**2)/(M22imagTrefoil1Err**2+M22realTrefoil2Err**2) # for x decenter
    M22realTrefoil1 = b[26] 
    M22imagTrefoil2 = b[49] 
    M22realTrefoil1Err = bErr[26] 
    M22imagTrefoil2Err = bErr[49] 
    if weighted == False:
        M22TrefoilYshift = 0.5*(M22realTrefoil1 - M22imagTrefoil2) # for y decenter
    else:
        M22TrefoilYshift = (M22realTrefoil1*M22imagTrefoil2Err**2 - M22imagTrefoil2*M22realTrefoil1Err**2)/(M22imagTrefoil2Err**2+M22realTrefoil1Err**2) # for y decenter
    M20defocus = b[4]    # for defocus
    M22realComa2 = b[28] 
    M22imagComa1 = b[47]
    M22realComa2Err = bErr[28] 
    M22imagComa1Err = bErr[47]
    if weighted == False:
        M22ComaXtilt = 0.5*(M22realComa2+M22imagComa1) # for x-tilt
    else:
        M22ComaXtilt = (M22realComa2*M22imagComa1Err**2+M22imagComa1*M22realComa2Err**2)/(M22imagComa1Err**2+M22realComa2Err**2) # for x-tilt
    M22realComa1 = b[27] 
    M22imagComa2 = b[48]
    M22realComa1Err = bErr[27] 
    M22imagComa2Err = bErr[48]
    if weighted == False:
        M22ComaYtilt = 0.5*(M22realComa1 - M22imagComa2) # for y-tilt
    else:
        M22ComaYtilt = (M22realComa1*M22imagComa2Err**2 - M22imagComa2*M22realComa1Err**2)/(M22imagComa2Err**2+M22realComa1Err**2) # for y-tilt
    x = -3.0063 * M22TrefoilXshift -0.0053
    y = -2.9318 * M22TrefoilYshift - 0.0005
    z = 0.4046 * M20defocus - 0.0705
    xtilt = 1075.8892 * M22ComaXtilt - 0.4876
    ytilt = -1064.6332 * M22ComaYtilt - 0.1234
    return np.array([x*1000,y*1000,z*1000,xtilt,ytilt])

def dispM202Coeff(betaAll=None,betaErrAll=None,hexinfo=None):
    ind = np.arange(len(betaAll[0]))
    momname = ('M20','M22.Real','M22.imag')
    fmtarr = ['bo-','ro-','go-']
    if betaErrAll == None:
        betaErrAll = np.zeros(len(ind))
    pl.figure(figsize=(17,7))
    for i in range(3):
        pl.subplot(4,1,i+1)
        pl.errorbar(ind[1:],betaAll[i][1:],yerr = betaErrAll[i][1:],fmt=fmtarr[i])
        pl.grid()
        pl.xlim(-1,len(betaAll[i])+1)
        if i==0 and hexinfo != None:
            pl.title('Hexapod deviation from perfect: x'+str(round(hexinfo[0],2))+'  y:'+str(round(hexinfo[1],2))+'  z:'+str(round(hexinfo[2],2))+'  xtilt:'+str(round(hexinfo[3],2))+'  ytilt:'+str(round(hexinfo[4],2)))
        #pl.ylim(min(betaAll[i][1:])-0.01,max(betaAll[i][1:])+0.01)
        pl.xticks(ind,('','','','','','','','','','','','','','','','','','','',''))
        pl.ylim(-0.3,0.3)
        pl.ylabel(momname[i])
    pl.xticks(ind,('Piston','Tip','Tilt','Astignism','Defocus','Astignism','Trefoil','Coma','Coma','Trefoil','Ashtray','Astigm.5th','Spherical','Astigm.5th','Ashtray','16','17','18','19','20'),rotation=90)
    pl.xlabel('Zernike Coefficients')
    return '---done!---'



def des_img_analysis(img_name=None):
    '''
    this function generate the weghted moments whiser plot, measure the R50, whisker, whisker rms
    '''
    catname = img_name[0:-5]+'_star_catalog.fits'
    if not os.path.isfile(catname):
        os.system('getstar.py '+img_name)
    imghdu = pf.open(img_name)
    cathdu = pf.open(catname)
    expid = img_name[6:14]
    dimmfwhm = pf.getheader(img_name,0)['dimmsee']
    kernelSigma = np.sqrt(dimmfwhm**2+0.55**2)/2.35482
    hexposhdr = pf.getheader(img_name,0)['telfocus']
    bcamDX = pf.getheader(img_name,0)['BCAMDX']
    bcamDY = pf.getheader(img_name,0)['BCAMDY']
    bcamAX = pf.getheader(img_name,0)['BCAMAX']
    bcamAY = pf.getheader(img_name,0)['BCAMAY']
    fltr = pf.getheader(img_name,0)['filter'][0]
    data=[]
    ccdpos=[]
    r50Sex = []
    nstarall = 0
    for i in range(1,63):
        print i
        img = imghdu[i].data
        cat = cathdu[i].data
        x = cat.XWIN_IMAGE
        y = cat.YWIN_IMAGE
        rad = cat.FLUX_RADIUS
        mag = cat.MAG_AUTO
        flag = cat.FLAGS
        bkg = cat.BACKGROUND
        fwhm_sex = cat.FWHM_IMAGE
        starFwhm = selectStar(mag,fwhm_sex)
        ok = (np.abs(fwhm_sex - starFwhm) < 0.4)*(x>100)*(x<2050)*(y>100)*(y<4100)*(flag == 0)*(mag<=-11.5)*(mag>-14.5)
        nstar = len(mag[ok])
        nstarall = nstarall + nstar
        print '--- Nstars selected: '+str(nstar)+'---'
        xccd = eval(imghdu[i].header['extname'])[1]
        yccd = eval(imghdu[i].header['extname'])[2]
        if ok.any():
            bkg = bkg[ok]
            x=x[ok]
            y=y[ok]
            stamp = getStamp(data=img,xcoord=x,ycoord=y,Npix=25)
            data.append(measureIQstamp(stamp,bkg,2.))
            ccdpos.append([xccd,yccd])
            r50Sex.append(robust_mean(rad[ok]))
        else:
            data.append([0.,0.,0.,0.,0.,0.])
            ccdpos.append([xccd,yccd])
            r50Sex.append(0.)
    data = np.array(data)
    ccdpos = np.array(ccdpos)
    #----make the whisker plot ----
    whisker4QReduce(X2WIN_IMAGE=data[:,0],Y2WIN_IMAGE=data[:,1],XYWIN_IMAGE=data[:,2])
    pl.figtext(0.4,0.96,'expid: '+str(expid)+' ,'+fltr+'-band')
    pl.savefig('moments_whisker_'+expid+'.png')
    pl.close()
    M20 = data[:,0] + data[:,1]
    M22real = data[:,0] - data[:,1]
    M22imag = 2.*data[:,2]
    #---now use the weight kernel corrected moments to calculate whiskers --
    datamean = np.array([robust_mean(data[:,3]),robust_mean(data[:,4]),robust_mean(data[:,5])])
    datasubmean = data[:,3:6] - datamean
    r50Sex = robust_mean(np.array(r50Sex))*scale
    whk = ((datamean[0]-datamean[1])**2 + (2.*datamean[2])**2)**(0.25)*scale
    phi = np.rad2deg(0.5*np.arctan2(2.*datamean[2],(datamean[0]-datamean[1])))
    whkrms = (robust_mean((datasubmean[:,0] - datasubmean[:,1])**2 + 4.*datasubmean[:,2]**2))**(0.25)*scale
    np.savetxt('desIQ_measure_'+expid+'_'+fltr+'.txt',[int(expid),whk,phi,whkrms,r50Sex],fmt=['%i','%10.5f','%10.5f','%10.5f','%10.5f'],delimiter=',')
    print 'expid -----whisker----whisker Position Angle -----whisker rms ----- r50 ------'
    print expid, whk, phi, whkrms, r50Sex
    #---the hexapod adjustment using M20,M22---
    beta=[]
    betaErr=[]
    betaM20 = zernikeFit(ccdpos[:,0],ccdpos[:,1],M20,max_order=20)
    beta.append(betaM20[0])
    betaErr.append(betaM20[1])
    betaM22real = zernikeFit(ccdpos[:,0],ccdpos[:,1],M22real,max_order=20)
    beta.append(betaM22real[0])
    betaErr.append(betaM22real[1])
    betaM22imag = zernikeFit(ccdpos[:,0],ccdpos[:,1],M22imag,max_order=20)
    beta.append(betaM22imag[0])
    betaErr.append(betaM22imag[1])
    betaforplot = beta
    betaErrforplot = betaErr
    beta=np.array(beta)
    betaErr = np.array(betaErr)
    beta=beta.flatten()
    betaErr = betaErr.flatten()
    hexHao = hexapodPosition(beta,betaErr,weighted=False)
    hexposhdr = np.array(hexposhdr.split(',')).astype(float)[0:5]
    hexBCAM = np.array([bcamDX,bcamDY,-999,bcamAX,bcamAY])
    np.savetxt('hexapod_position_'+expid+'_'+fltr+'.txt',[hexposhdr,hexHao,hexBCAM],fmt='%10.5f')
    #hexHao = hexapod_multilinear(beta)
    dispM202Coeff(betaAll = betaforplot, betaErrAll = betaErrforplot,hexinfo=hexHao)
    pl.figtext(0.4,0.96,'expid: '+str(expid)+' ,'+fltr+'-band')
    pl.savefig('zernike_coeff_'+expid+'_'+fltr+'.png')    
    return '----finished one image ----'

def correctMoments(Mcc=None, Mrr=None, Mrc=None,measureSigma=None,targetSigma=2.):
    '''
    This function correct the weighted moments to remove the weight effect and re-assign it to a new weight
    All based on the Gaussian assumption
    '''
    for i in range(len(Mcc)):
        Cm = np.matrix([[Mcc[i],Mrc[i]],[Mrc[i],Mrr[i]]])
        Cw = np.matrix([[measureSigma[i]**2,0.],[0.,measureSigma[i]**2]])
        Cwtarget = np.matrix([[targetSigma**2,0.],[0.,targetSigma**2]])
        Cimg = (Cm.I - Cw.I).I
        Cmtarget = (Cimg.I + Cwtarget.I).I
        Mcc[i] = Cmtarget[0,0]
        Mrr[i] = Cmtarget[1,1]
        Mrc[i] = Cmtarget[0,1]
    return Mcc, Mrr, Mrc
    
    
def firstcut_moments_analysis(expid=None, catdir=None):
    '''
    this function generate the weghted moments whiser plot, 
    measure the R50, whisker, whisker rms based on the firstcut catalog
    '''
    if expid == None or catdir == None:
            return('syntax: firstcut_moments_analysis(expid,catlog_dir)')
    ff = gl.glob(catdir+'/DECam_00'+expid+'_??_cat.fits')
    ff.sort()
    if len(ff) == 62:
        return -999., -999., -999., -999.,-999. 
    data=[]
    ccdpos=[]
    r50Sex = []
    nstarall = 0
    for i in range(62):
        print i
        if i == 60:
            xccd = np.genfromtxt('xpos_ypos_fp.txt').T[0][i]
            yccd = np.genfromtxt('xpos_ypos_fp.txt').T[1][i]
            data.append([0.,0.,0.])
            ccdpos.append([xccd,yccd])
            r50Sex.append(0.)
            continue
        elif i == 61:
            i = i-1
        cat = firstcutStar(pf.getdata(ff[i],2))
        rad = cat.FLUX_RADIUS
        nstar = len(rad)
        measureSigma = rad*2./2.35482
        Mcc, Mrr, Mrc = correctMoments(cat.X2WIN_IMAGE,cat.Y2WIN_IMAGE,cat.XYWIN_IMAGE,measureSigma,2.)
        nstarall = nstarall + nstar
        print '--- Nstars selected: '+str(nstar)+'---'
        xccd = np.genfromtxt('xpos_ypos_fp.txt').T[0][i]
        yccd = np.genfromtxt('xpos_ypos_fp.txt').T[1][i]
        if rad.any():
            data.append([robust_mean(Mcc),robust_mean(Mrr),robust_mean(Mrc)])
            ccdpos.append([xccd,yccd])
            r50Sex.append(robust_mean(rad))
        else:
            data.append([0.,0.,0.])
            ccdpos.append([xccd,yccd])
            r50Sex.append(0.)
    data = np.array(data)
    ccdpos = np.array(ccdpos)
    #----make the whisker plot ----
    whisker4QReduce(X2WIN_IMAGE=data[:,0],Y2WIN_IMAGE=data[:,1],XYWIN_IMAGE=data[:,2])
    pl.savefig('moments_whisker_'+expid+'.png')
    pl.close()
    M20 = data[:,0] + data[:,1]
    M22real = data[:,0] - data[:,1]
    M22imag = 2.*data[:,2]
    datamean = np.array([robust_mean(data[:,0]),robust_mean(data[:,1]),robust_mean(data[:,2])])
    datasubmean = data - datamean
    r50Sex = robust_mean(np.array(r50Sex))*scale
    whk = ((datamean[0]-datamean[1])**2 + (2.*datamean[2])**2)**(0.25)*scale
    phi = np.rad2deg(0.5*np.arctan2(2.*datamean[2],(datamean[0]-datamean[1])))
    whkrms = (robust_mean((datasubmean[:,0] - datasubmean[:,1])**2 + 4.*datasubmean[:,2]**2))**(0.25)*scale
    p.dump([int(expid),whk,phi,whkrms,r50Sex],open('desIQ_measures_'+expid+'.p','w'))
    print 'expid -----whisker----whisker Position Angle -----whisker rms ----- r50 ------'
    print expid, whk, phi, whkrms, r50Sex
    #---the hexapod adjustment using M20,M22---
    beta=[]
    betaErr=[]
    betaM20 = zernikeFit(ccdpos[:,0],ccdpos[:,1],M20,max_order=20)
    beta.append(betaM20[0])
    betaErr.append(betaM20[1])
    betaM22real = zernikeFit(ccdpos[:,0],ccdpos[:,1],M22real,max_order=20)
    beta.append(betaM22real[0])
    betaErr.append(betaM22real[1])
    betaM22imag = zernikeFit(ccdpos[:,0],ccdpos[:,1],M22imag,max_order=20)
    beta.append(betaM22imag[0])
    betaErr.append(betaM22imag[1])
    betaforplot = beta
    betaErrforplot = betaErr
    beta=np.array(beta)
    betaErr = np.array(betaErr)
    beta=beta.flatten()
    betaErr = betaErr.flatten()
    hexHao = hexapodPosition(beta,betaErr,weighted=False)
    #hexHao = hexapod_multilinear(beta)
    dispM202Coeff(betaAll = betaforplot, betaErrAll = betaErrforplot,hexinfo=hexHao)
    pl.savefig('zernike_coeff_'+expid+'.png')    
    return '----finished one image ----'

    
