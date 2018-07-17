import numpy as np
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
from scipy.spatial import cKDTree

rad = 0.5 # matching radii, arcsec
err = 0.3 # mag error, Vega
crd = 0.5 # crowding effect, Vega

def DoAll(infile='14786_N6846-SN2017.gst.fits',outfile='targets1'):
    ref = ascii.read('n6946_CAT')
    tol = rad/3600
    filt  = ['F606W','F814W']
    ra2,dec2,mm1,ee1,dd1,mm2,ee2,dd2 = ref['col1'],ref['col2'],ref['col3'],ref['col4'],ref['col5'],ref['col6'],ref['col7'],ref['col8']
    t = (ee1<err)&(ee2<err)&(dd1<crd)&(dd2<crd)
    ra2,dec2,mm1,ee1,mm2,ee2=ra2[t],dec2[t],mm1[t],ee1[t],mm2[t],ee2[t]
    X,Y,RA,DEC,vega_mags,mag_errors,detector = read_phot_fits_table(infile,filt)
    t = (((vega_mags[0]-vega_mags[1])>2)&(vega_mags[0]<29)&(vega_mags[1]<24.5))
    ra1,dec1,mag1,mag2,err1,err2 = RA[t],DEC[t],vega_mags[0][t],vega_mags[1][t],mag_errors[0][t],mag_errors[1][t]
    in1 = matchLists(tol,ra1,dec1,ra2,dec2)
    ra1,dec1,mag1,mag2,err1,err2 = ra1[in1!=-1],dec1[in1!=-1],mag1[in1!=-1],mag2[in1!=-1],err1[in1!=-1],err2[in1!=-1]
    mm1,mm2,ee1,ee2 = mm1[in1[in1!=-1]],mm2[in1[in1!=-1]],ee1[in1[in1!=-1]],ee2[in1[in1!=-1]]
    return write_target_list(ra1,dec1,mag1,err1,mag2,err2,mm1,ee1,mm2,ee2,outfile)

def read_phot_fits_table(filename,filt):
    photTable = fits.open(filename)
    detector = photTable[0].header['CAMERA']
    data = photTable[1].data; del photTable
    vega_mags  = [data[filt[0]+'_VEGA'], data[filt[1]+'_VEGA']]
    mag_errors = [data[filt[0]+'_ERR'], data[filt[1]+'_ERR']]
    X,Y,RA,DEC = data['X'],data['Y'],data['RA'],data['DEC']
    return X,Y,RA,DEC,vega_mags,mag_errors,detector

''' Quick match; returns index of 2nd list coresponding to position in 1st '''
def matchLists(tol,x1,y1,x2,y2):
    d1 = np.empty((x1.size, 2))
    d2 = np.empty((x2.size, 2))
    d1[:,0],d1[:,1] = x1,y1
    d2[:,0],d2[:,1] = x2,y2
    t = cKDTree(d2)
    tmp, in1 = t.query(d1, distance_upper_bound=tol)
    in1[in1==x2.size] = -1
    return in1

def write_target_list(ra,dec,m1,err1,m2,err2,m3,err3,m4,err4,fileroot):
    outfile = fileroot+'.txt'
    tab = [ra,dec,m1,err1,m2,err2,m3,err3,m4,err4]
    nms = ('RA','DEC','F606W','ER606','F814W','ER814','IRAC1','ERR36','IRAC2','ERR45')
    fmt = {'RA':'%12.7f','DEC':'%12.7f','F606W':'%8.3f','ER606':'%8.3f',
    'F814W':'%8.3f','ER814':'%8.3f','IRAC1':'%8.3f','ERR36':'%8.3f','IRAC2':'%8.3f','ERR45':'%8.3f'}
    t = Table(tab, names=nms)
    ascii.write(t, outfile, format='fixed_width', delimiter='', formats=fmt, overwrite=True)
    return print('Wrote out: ', outfile)

###

DoAll('14786_N6846-SN2017.gst.fits','targets1')

DoAll('14786_N6946-111ne-26900.gst.fits','targets2')

DoAll('20180119_Vanisher.gst.fits','targets3')

