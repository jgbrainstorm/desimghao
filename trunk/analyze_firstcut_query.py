# this code calcuate the quantities from the firstcut catlogs from Robert's query. 
from func_def import *
import glob as gl

f =gl.glob('/home/jghao/research/data/firstcutcat/iqcat_robert_feb_2013/*.csv') 
f.sort()

expid = []
n=len(f)

for i in range(n):
    b=np.genfromtxt(f[i],delimiter=',')
    for ext in range(1,63):
        idx = b[:,3] = ext
        flux_radius = b[idx,8]
        ixx = b[idx,10]
        iyy = b[idx,11]
        ixy = b[idx,12]
        sigixx = b[idx,].ERRX2_IMAGE
        sigiyy = b.ERRY2_IMAGE
        sigixy = b.ERRXY_IMAGE
        wl = ( (ixx-iyy)**2 + (2.*ixy)**2 )**0.25
        sigma_wl = 0.5 * ( (ixx-iyy)**2 * sigixx**2 + 16. * ixy**2 * sigixy**2 )**0.5 / wl**3
        ok = (sigma_wl <= sigma_maxval)
        ixx = ixx[ok]
        ixy = ixy[ok]
        iyy = iyy[ok]
        size = np.log(ixx+iyy)
        wl = ( (ixx-iyy)**2 + (2.*ixy)**2 )**0.25 
        good = remOutlierIdx(size)*remOutlierIdx(wl)
        Mcc = ixx[good].mean()
        Mrr = iyy[good].mean()
        Mrc = ixy[good].mean()
        fluxrad = robust_mean(b.FLUX_RADIUS)
            fwhmworld =  robust_mean(b.FWHM_WORLD)
            data.append([Mcc,Mrr,Mrc,fluxrad,fwhmworld])
    e