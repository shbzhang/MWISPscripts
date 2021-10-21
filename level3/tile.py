' mwisp module: tile '

__version__ = '1.0'
__author__ = 'Shaobo Zhang'

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

#modify projection frame to other galactic longitude
def modifyframe(hdr, gl_ref=180):
	if hdr['CTYPE1'] == 'GLON-CAR':
		hdr['CRVAL1'] = gl_ref
		if hdr['CRPIX1']<0:
			hdr['CRPIX1'] += 43200
		hdr['CRPIX1'] -= gl_ref/0.0083333333333
	return hdr

#use linear convert for non-celetial axes
def linear_WCS(hdr, axis):
	return {key:hdr[key+axis] for key in ['NAXIS','CTYPE','CRVAL','CDELT','CRPIX','CROTA']}

def linear_pixel_to_world(p, wcs, base=0):
	return (np.array(p)-wcs['CRPIX']+1-base)*wcs['CDELT']+wcs['CRVAL']

def linear_world_to_pixel(w, wcs, base=0):
	return (np.array(w)-wcs['CRVAL'])/wcs['CDELT']+wcs['CRPIX']-1+base

def tile(files, output='tile.fits', gl_ref=120):
	'''
	Tile images.
	Could be used to mosaic integrated intensity maps/lv-maps/cubes from separately mosaicked datacube.

	Parameters
	----------
	files : list
		a list of image name in fits format
	output : str, optional
		output file name. Default name is 'tile.fits'.
	Versions
	--------
	Nov,09,2011,v1.0
	Examples
	--------
	>>> from glob import glob
	>>> #tile lv maps
	>>> tile(glob('*_U_lvmap.fits'))
	>>> #tile integrated maps
	>>> tile(glob('*_U_m0.fits'))
	'''
	nfile = len(files)
	#Syntax prompt
	if nfile < 2:
		print("Syntax - tile(files, outputs='tile.fits')")
		return
	#open a fits as basic template
	thdr = fits.getheader(files[0])
	thdr = modifyframe(thdr, gl_ref)
	try:
		twcs = WCS(thdr, naxis=2)
	except:
		twcs = [linear_WCS(thdr,'1'),linear_WCS(thdr,'2')]
	#find xy range for each file
	xran = np.ndarray((nfile,2))
	yran = np.ndarray((nfile,2))
	for i,afile in enumerate(files):
		hdr = fits.getheader(afile)
		hdr = modifyframe(hdr, gl_ref)
		if type(twcs) is not list:
			wcs = WCS(hdr, naxis=2)
			ad = wcs.pixel_to_world([0, hdr['NAXIS1']-1], [0, hdr['NAXIS2']-1])
			xy = twcs.world_to_pixel(ad)
		else:
			a = linear_pixel_to_world([0, hdr['NAXIS1']-1], linear_WCS(hdr,'1'))
			d = linear_pixel_to_world([0, hdr['NAXIS2']-1], linear_WCS(hdr,'2'))
			xy = [linear_world_to_pixel(a, twcs[0]), linear_world_to_pixel(d, twcs[1])]
		xran[i] = xy[0]
		yran[i] = xy[1]
	xran=np.round(xran).astype(int)
	yran=np.round(yran).astype(int)
	print(xran)

	#header
	#thdr['NAXIS'] = 2
	thdr['NAXIS1'] = xran.max()-xran.min()+1
	thdr['NAXIS2'] = yran.max()-yran.min()+1
	thdr['CRPIX1'] -= xran.min()
	thdr['CRPIX2'] -= yran.min()
	xran -= xran.min()
	yran -= yran.min()
	mshp = []
	for i in range(thdr['NAXIS'],0,-1):
		mshp.append(thdr['NAXIS%1i' % i])
	mdat = np.zeros(mshp,dtype=np.float32)-1e3
	#print(mdat.shape)
	for i,afile in enumerate(files):
		dat = fits.open(afile)[0].data
		msk = np.isnan(dat)
		dat[msk] = 0
		mdat[...,yran[i,0]:yran[i,1]+1,xran[i,0]:xran[i,1]+1] = mdat[...,yran[i,0]:yran[i,1]+1,xran[i,0]:xran[i,1]+1]*msk + dat*(1-msk)
		#print(afile,xran[i])
	mdat[mdat == -1e3] = np.nan

	mhdu = fits.PrimaryHDU(mdat, thdr)
	mhdu.writeto(output, overwrite=True)


if __name__ == '__main__':
	import glob
	files = glob.glob('*_U_lvmap.fits')
	tile(files)
