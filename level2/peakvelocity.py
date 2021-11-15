' mwisp module: peakvelocity '

__version__ = '1.0'
__author__ = 'Shaobo Zhang'

import os, re, warnings, time
import numpy as np
import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube
from astropy.wcs import WCS

def peakvelocity(cubefile=None, vrange=[-np.inf, np.inf], unit=u.km/u.s):
	'''
	Calculate temperature and velocity of peak for datacube
	Parameters
	----------
	files : str
		the datacube filename for calculation.
	vrange : list, optional
		contains the range over which calculation should be done.
	unit : astropy.units, optional
		unit of the velocity. Default is km/s
	--------
	Nov,10,2021,v1.0
		rewrite in python3
		add 'unit' keyword.
	Examples
	--------
	>>> peakvelocity('xxx.fits', [-20,20])
	'''
	if cubefile is None:
		print, "Syntax - peakvelocity(fitsfile, vrange=[v1, v2], unit='km/s')"
		return

	time_start = time.time()
	print('Parameters\n----------')

	#cubefile
	if os.path.exists(cubefile):
		print('%-12s = %s' % ('Datacube', cubefile))
	else:
		print('Error - datacube file does not exist.')
		return

	#vrange
	if len(vrange) >= 2:
		vrange = [min(vrange)*u.Unit(unit), max(vrange)*u.Unit(unit)]
		print('%-12s = %s to %s' % ('Range', *vrange))
	else:
		print('Error - invalid range.')
		return

	print('\nCalculate peak and velocity\n----------------')
	warnings.filterwarnings('ignore', r'Could not parse unit')
	
	#read cube and slab velocity
	#hdu = SpectralCube.read(cubefile).with_spectral_unit(unit)
	hdu = fits.open(cubefile)[0]
	header = hdu.header

	hdu = SpectralCube(data=np.squeeze(hdu.data),wcs=WCS(hdu.header,naxis=3)).with_spectral_unit(unit)
	velocity = hdu.spectral_axis
	if vrange[0]<velocity.min(): vrange[0] = velocity.min()
	if vrange[1]>velocity.max(): vrange[1] = velocity.max()
	subhdu = hdu.spectral_slab(*vrange)

	#calculate
	subcube = np.squeeze(subhdu._data)
	tpeakdata = np.nanmax(subcube, axis=0)
	subcube[~np.isfinite(subcube)] = -np.inf
	peakidx = np.argmax(subcube, axis=0)
	vpeakdata = subhdu.spectral_axis.value[peakidx].astype(np.float32)

	for kwd in ['NAXIS','CTYPE','CRVAL','CRPIX','CDELT']:
		header.remove(kwd+'3')
		header.remove(kwd+'4', ignore_missing = True)

	#output
	basename = os.path.basename(cubefile)
	tpeakname = '_Tpeak'.join(os.path.splitext(basename))
	header['BUNIT'] = 'K'
	tpeak = fits.PrimaryHDU(data=tpeakdata, header=header)
	tpeak.writeto(tpeakname, overwrite=True)

	vpeakname = '_Vpeak'.join(os.path.splitext(basename))
	header['BUNIT'] = str(unit)
	vpeak = fits.PrimaryHDU(data=vpeakdata, header=header)
	vpeak.writeto(vpeakname, overwrite=True)

	time_end = time.time()
	print('Done in %f seconds.' % (time_end-time_start))


if __name__ == '__main__':
	peakvelocity('/Users/sz268601/Work/DeepOutflow/procedure/prediction/028.825+1.275+000_U.fits',[10,20])
