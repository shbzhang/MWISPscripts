' mwisp module: cubemoment '

__version__ = '1.0'
__author__ = 'Shaobo Zhang'

import os, math, time
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.coordinates as coord

def stdrange(ran, naxis):
#standardize 2-element range
    lower = int(ran.min().round())
    upper = int(ran.max().round())
    if (lower >= naxis) | (upper < 0):
        return None
    if lower < 0:
        lower = 0
    if upper >= naxis:
        upper = naxis-1
    return [lower, upper]

def c2v(header, channel):
    channel = np.array(channel)
    nc = header['NAXIS3']
    v0 = header['CRVAL3']
    c0 = header['CRPIX3']
    dv = header['CDELT3']
    velocity = (channel-c0+1)*dv+v0
    return velocity

def v2c(header, velocity):
    velocity = np.array(velocity)
    nc = header['NAXIS3']
    v0 = header['CRVAL3']
    c0 = header['CRPIX3']
    dv = header['CDELT3']
    channel = (velocity-v0)/dv+c0-1
    return channel

def cubemoment(cubefile=None, crange=None, direction='v', threshold=-np.inf, zeroth_only=False,
    coveragefile=None, goodlooking=False, rmsfile=None, outname=None):
    '''
    Calculate moment for a datacube

    Parameters
    ----------
    cubefile : str
        datacube filename for calculation
    crange : list of two float
        range over which the moments should be computed.
        For 'L' or 'B' direction, the unit should be arcdeg; for 'V', the unit should be km/s
    direction : str, optional
        default calculation is done along the "Velocity" direction. 
        Other directions as "L" or "B" could be provided to calculate along the gl or gb.
    threshold : float, optional
        the minimum intensity under which a pixel should no longer be considered.
    zeroth_only : bool, optional
        set to True if you only need the 0th moment.
    coveragefile : str, optional
        the filename of coverage map for the datacube.
        default is cubefile(no suffix)+'_coverage.fits'
    goodlooking : bool, optional
        set to True to active the goodlooking mode.
        in the "goodlooking" mode, only pixels with at least 3 consecutive channels above 3 sigma are kept.
        this mode will make the 1st and 2nd moments more reliable.
    rmsfile: str, optional
        the filename of sigma map for the datacube, only use in goodlooking mode.
        default is cubefile(no suffix)+'_rms.fits'.
    
    Versions
    --------
    May,08,2012,v1.0
    Dec,25,2013,v1.1
        goodlooking mode: only keep those pixels with at least 3 "consecutive" channels above 3 sigma
    Jun,24,2014,v1.2
        fix minor error when apply to all FITS file
    Nov,06,2018,v1.3
        rewrite in python3

    Examples
    --------
    For general use:
    >>> cubemoment('cube.fits', [-10,10])
    >>> cubemoment('cube.fits', [-10,10], threshold=1.0, zeroth_only=True, coveragefile='coverage.fits')
    Integrate along b direction to create an l-v map:
    >>> cubemoment('cube.fits', [-2,-1], direction='b')
    Goodlooking mode:
    >>> cubemoment('cube.fits', [-10,10], goodlooking=True, rmsfile='rms.fits', coveragefile='coverage.fits')
    For mosaiced PMO survey data produce by mosaic.pro
    >>> cubemoment, 'mosaic_U.fits', [-10,10], /goodlooking
    '''

    #standardize and test input
    if (cubefile == None) or (crange == None):
        print("Syntax - mwisp.cubemoment(cubefile, crange, direction='v', threshold=-inf, zeroth_only=False,")
        print("             coveragefile=None, goodlooking=False, rmsfile=None, outname=None)")
        return

    time_start = time.time()
    print('Parameters\n----------')

    #cubefile
    if os.path.exists(cubefile):
        print(str('%-12s = %s' % ('Datacube', cubefile)))
    else:
        print('Error - datacube file does not exist.')
        return

    #crange
    if len(crange) >= 2:
        crange = crange[0:2]
        print(str('%-12s = %f %f' % ('Range', crange[0], crange[1])))
    else:
        print('Error - invalid range.')
        return

    #direction
    direction = direction.upper()
    if direction in 'LBV':
        print(str('%-12s = %s' % ('Direction', direction)))
    else:
        print("Error - direction must be 'L', 'B', or 'V'.")
        return

    #threshold & zeroth_only
    print(str('%-12s = %s' % ('Threshold', threshold)))
    print(str('%-12s = %s' % ('Zeroth_only', zeroth_only)))

    #coveragefile
    if coveragefile == None:
    #if coverage file is not provided, try default coverage file name
        pathname, ext = os.path.splitext(cubefile)
        defaultcoveragefile = pathname+'_coverage'+ext
        if os.path.exists(defaultcoveragefile): #existance test
            cubhdr = fits.getheader(cubefile)
            covhdr = fits.getheader(defaultcoveragefile)
            if (cubhdr['NAXIS1'] == covhdr['NAXIS1']) and (cubhdr['NAXIS2'] == covhdr['NAXIS2']): #dimension test
                coveragefile = defaultcoveragefile
        print(str('%-12s = %s' % ('Coverage', coveragefile)))
    else:
        print(str('%-12s = %s' % ('Coverage', coveragefile)))
        if not os.path.exists(coveragefile): #existance test
            print('Error - coveragefile does not exist.')
            return
        else:
            cubhdr = fits.getheader(cubefile)
            covhdr = fits.getheader(coveragefile)
            if (cubhdr['NAXIS1'] != covhdr['NAXIS1']) or (cubhdr['NAXIS2'] != covhdr['NAXIS2']): #dimension test
                print('Error - dimension of coveragefile and cubefile does not match.')
                return

    #goodlooking & rmsfile
    if goodlooking:
        #test rms file only in goodlooking mode
        print(str('%-12s = %s' % ('Goodlooking', goodlooking)))
        if rmsfile == None:
        #if rms file is not provided, try default rms file name
            pathname, ext = os.path.splitext(cubefile)
            defaultrmsfile = pathname+'_rms'+ext
            if os.path.exists(defaultrmsfile): #existance test
                cubhdr = fits.getheader(cubefile)
                rmshdr = fits.getheader(defaultrmsfile)
                if (cubhdr['NAXIS1'] == rmshdr['NAXIS1']) and (cubhdr['NAXIS2'] == rmshdr['NAXIS2']): #dimension test
                    rmsfile = defaultrmsfile
            print(str('%-12s = %s' % ('RMS', rmsfile)))
            if rmsfile == None:
                print('Error - a rmsfile must be given in goodlooking mode.')
                return
        else:
            print(str('%-12s = %s' % ('RMS', rmsfile)))
            if not os.path.exists(rmsfile): #existance test
                print('Error - rmsfile does not exist.')
                return
            else:
                cubhdr = fits.getheader(cubefile)
                covhdr = fits.getheader(rmsfile)
                if (cubhdr['NAXIS1'] != covhdr['NAXIS1']) or (cubhdr['NAXIS2'] != covhdr['NAXIS2']): #dimension test
                    print('Error - dimension of rmsfile and cubefile does not match.')
                    return

    #outname
    if outname == None:
    #if outname is not provided, use the basename of input cubefile name
        outname = os.path.basename(cubefile)
        outname, ext = os.path.splitext(outname)
    print(str('%-12s = ./%s_*.fits' % ('Output', outname)))

    print('\nCalculate Moment\n----------------')
    print('Loading fits file......(this may take a while!)')
    cubhdu = fits.open(cubefile)[0]
    cubhdu.data = np.squeeze(cubhdu.data)
    w = WCS(cubhdu.header)

    if coveragefile != None:
        covhdu = fits.open(coveragefile)[0]
        for c in range(cubhdu.data.shape[0]):
            cubhdu.data[c,:,:][covhdu.data == 0] = np.nan


    if direction == 'L':
        print('Creating L-V map.')
        #range crange --> cposition -wcs-> px --std--> prange --wcs--> realposition --> realrange
        cposition = coord.SkyCoord(crange, np.repeat(cubhdu.header['CRVAL2'],2), unit='deg', frame='galactic')
        px,py = cposition.to_pixel(w, origin=0)
        prange = stdrange(px, cubhdu.header['NAXIS1'])
        if prange == None:
            print('Error: Range is out of the given datacube.')
            return
        realposition = coord.SkyCoord(0,0,unit='deg').from_pixel(prange, py, w)
        realrange = realposition.l.value
        #subcube and moment
        subcube = cubhdu.data[:,:,prange[0]:prange[1]+1]
        subcube[subcube < threshold] = 0
        nan = ~(np.isfinite(subcube).any(axis=2))
        if goodlooking:
            rmshdu = fits.open(rmsfile)[0]
            #only keep those pixels with at least 3 "consecutive" channels above 3 sigma
            mskcube = np.ndarray(subcube.shape, dtype=np.bool)
            for i in range(subcube.shape[0]):
                mskcube[i,:,:] = subcube[i,:,:] > rmshdu.data[:,prange[0]:prange[1]+1]*3
            mskcube = mskcube & np.roll(mskcube, 1, axis=0) & np.roll(mskcube, 2, axis=0)
            mskcube = mskcube | np.roll(mskcube, -1, axis=0) | np.roll(mskcube, -2, axis=0)
            subcube *= mskcube
        mom0 = np.nansum(subcube, axis=2) * np.abs(cubhdu.header['CDELT1'])
        mom0[nan] = np.nan
        #output
        cubhdu.header['BUNIT'] = 'K * arcdeg'
        for kwd in ['NAXIS','CTYPE','CRVAL','CRPIX','CDELT']:
            cubhdu.header[kwd+'1'] = cubhdu.header[kwd+'2']
            cubhdu.header[kwd+'2'] = cubhdu.header[kwd+'3']
            cubhdu.header.remove(kwd+'3')
            cubhdu.header.remove(kwd+'4', ignore_missing = True)
        cubhdu.header['HISTORY'] = str('Integrate L between %f and %f' % tuple(realrange))
        cubhdu.data = mom0
        cubhdu.writeto(outname+'_bvmap.fits', overwrite = True)

    elif direction == 'B':
        print('Creating B-V map.')
        #range
        cposition = coord.SkyCoord(np.repeat(cubhdu.header['CRVAL1'],2), crange, unit='deg', frame='galactic')
        px,py = cposition.to_pixel(w, origin=0)
        prange = stdrange(py, cubhdu.header['NAXIS2'])
        if prange == None:
            print('Error: Range is out of the given datacube.')
            return
        realposition = coord.SkyCoord(0,0,unit='deg').from_pixel(px, prange, w)
        realrange = realposition.b.value
        #subcube and moment
        subcube = cubhdu.data[:,prange[0]:prange[1]+1,:]
        subcube[subcube < threshold] = 0
        nan = ~(np.isfinite(subcube).any(axis=1))
        if goodlooking:
            rmshdu = fits.open(rmsfile)[0]
            #only keep those pixels with at least 3 "consecutive" channels above 3 sigma
            mskcube = np.ndarray(subcube.shape, dtype=np.bool)
            for i in range(subcube.shape[0]):
                mskcube[i,:,:] = subcube[i,:,:] > rmshdu.data[prange[0]:prange[1]+1,:]*3
            mskcube = mskcube & np.roll(mskcube, 1, axis=0) & np.roll(mskcube, 2, axis=0)
            mskcube = mskcube | np.roll(mskcube, -1, axis=0) | np.roll(mskcube, -2, axis=0)
            subcube *= mskcube
        mom0 = np.nansum(subcube, axis=1) * np.abs(cubhdu.header['CDELT2'])
        mom0[nan] = np.nan
        #output
        cubhdu.header['BUNIT'] = 'K * arcdeg'
        for kwd in ['NAXIS','CTYPE','CRVAL','CRPIX','CDELT']:
            cubhdu.header[kwd+'2'] = cubhdu.header[kwd+'3']
            cubhdu.header.remove(kwd+'3')
            cubhdu.header.remove(kwd+'4', ignore_missing = True)
        cubhdu.header['HISTORY'] = str('Integrate B between %f and %f' % tuple(realrange))
        cubhdu.data = mom0
        cubhdu.writeto(outname+'_lvmap.fits', overwrite = True)

    elif direction == 'V':
        print('Creating 0th moment map.')
        #range
        #from converter import c2v, v2c
        pc = v2c(cubhdu.header, np.array(crange)*1e3)
        prange = stdrange(pc, cubhdu.header['NAXIS3'])
        if prange == None:
            print('Error: Range is out of the given datacube.')
            return
        realrange = c2v(cubhdu.header, prange)/1e3
        #subchan/velo for moment 1/2
        subchan = np.arange(prange[0], prange[1]+1)
        subvelo = c2v(cubhdu.header, subchan)
        #subcube and moment
        subcube = cubhdu.data[prange[0]:prange[1]+1,:,:]
        subcube[subcube < threshold] = 0
        nan = ~(np.isfinite(subcube).any(axis=0))
        if goodlooking:
            rmshdu = fits.open(rmsfile)[0]
            #only keep those pixels with at least 3 "consecutive" channels above 3 sigma
            mskcube = np.ndarray(subcube.shape, dtype=np.bool)
            for i in range(subcube.shape[0]):
                mskcube[i,:,:] = subcube[i,:,:] > rmshdu.data*3
            mskcube = mskcube & np.roll(mskcube, 1, axis=0) & np.roll(mskcube, 2, axis=0)
            mskcube = mskcube | np.roll(mskcube, -1, axis=0) | np.roll(mskcube, -2, axis=0)
            subcube *= mskcube
        sumi = np.nansum(subcube, axis=0) 
        mom0 = sumi * np.abs(cubhdu.header['CDELT3'])
        mom0[nan] = np.nan
        #output
        cubhdu.header['BUNIT'] = 'K * km / s'
        for kwd in ['NAXIS','CTYPE','CRVAL','CRPIX','CDELT']:
            cubhdu.header.remove(kwd+'3')
            cubhdu.header.remove(kwd+'4', ignore_missing = True)
        cubhdu.header['HISTORY'] = str('Integrate V between %f and %f' % tuple(realrange))
        cubhdu.data = mom0 / 1e3   #K m/s -> K km/s
        cubhdu.writeto(outname+'_m0.fits', overwrite = True)
        if not zeroth_only:
            print('Creating 1st moment map.')
            #subcube and moment
            datacopy = subcube.copy()
            for i in range(len(subvelo)):
                datacopy[i,:,:] *= subvelo[i]
            sumi[sumi == 0] = np.nan
            mom1 = np.nansum(datacopy, axis=0) / sumi
            mom1[nan] = np.nan
            #output
            cubhdu.header['BUNIT'] = 'km / s'
            cubhdu.data = mom1 / 1e3
            cubhdu.writeto(outname+'_m1.fits', overwrite = True)

            print('Creating 2nd moment map.')
            #subcube and moment
            for i in range(len(subvelo)):
                subcube[i,:,:] *= (subvelo[i] - mom1)**2
            mom2 = np.sqrt(np.nansum(subcube, axis=0) / sumi * 8*np.log(2)) #to FWHM
            mom2[nan] = np.nan
            #output
            cubhdu.data = mom2 / 1e3
            cubhdu.writeto(outname+'_m2.fits', overwrite = True)
    time_end = time.time()
    print('Done in %f seconds.' % (time_end-time_start))


if __name__ == '__main__':
    cube='/Users/shaobo/Work/mwips/L935/mosaic_U_clipv.fits'
    #cover='/Users/shaobo/Work/mwips/L935/mosaic_coverage.fits'
    rms='/Users/shaobo/Work/mwips/L935/mosaic_U_rms.fits'
    #cubemoment(cube,[-15,15],direction='v',outname='test',zeroth_only=True)
    cubemoment(cube,[-15,15],direction='v',coveragefile=None,goodlooking=True,rmsfile=rms,zeroth_only=False,outname='test',threshold=-100)
