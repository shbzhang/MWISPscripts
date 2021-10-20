'''
;by ShaoboZhang (shbzhang@pmo.ac.cn)
;CUBERMS
;Calculate rms for a datacube
;Usage: cuberms, fitsfile, range, [window=[w1, w2, ..., wn]]
;Input:
;  cubefile: a string scalar that contains datacube filename for calculation
;  range: 2-element vector. Range over which the rms should be computed.
;Optional input keyword:
;  window: ranges avoided by baseline fitting.
;  silent: if set, then the display of the process description on the terminal will be suppressed 
;Example:
;  cuberms, 'cube.fits', [-200,200]
;  cuberms, 'cube.fits', [-100,100], window=[-50,-30,-10,10]
;History:
;Dec,23,2013,v1.0
;Oct,17,2014,v1.1
;  consider the situation that velocity range in cube is less then the given range.
'''
import os, re, warnings
import numpy as np
import astropy.units as u
from astropy.io import fits
from spectral_cube import SpectralCube

def readinfo(infofile, cellname):
	try:
		f = open(infofile,'r')
		for line in f:
			if cellname in line:
				break
		for line in f:
			searchmode = re.search('vrange.+?(-?\d+\.?\d*).+?(-?\d+\.?\d*)', line, re.I)
			if searchmode is not None:
				mode = [float(searchmode.group(1)),float(searchmode.group(2))]
			searchwindow = re.search('signal\s+range.+?(-?\d+\.?\d*).+?(-?\d+\.?\d*)', line, re.I)
			if searchwindow is not None:
				window = [float(searchwindow.group(1)),float(searchwindow.group(2))]
			if re.match('\s*\n',line):
				break
		f.close()
		return mode, window
	except:
		return [], []


def cuberms(*files, window=[], unit=u.km/u.s, silent=False):
	if len(files)<2:
		print("Syntax - cube(cubefile, infofile, silent=False)")
		print("Syntax - cube(cubefile, mode, window=[], silent=False)")
		return

	#assign files
	cubefile = files[0]
	cellname = os.path.splitext(os.path.basename(files[0]))[0]
	if type(files[1]) is str:
		#an info file
		mode, window = readinfo(files[1], cellname)
	else:
		#an range list
		mode = files[1]
	if len(mode)!=2:
		print("Error - invalid mode range")
		return

	#open cube
	if not silent:
		print("RMS calculation for '%s'' begins:" % cubefile)
		print('  Loading fits file......(this may take a while!)')
	try:
		cube = fits.open(cubefile)[0]
	except:
		print("Error - unable to read datacube '%s'" % cubefile)
		return

	#obtain spectral flag
	velocity = (np.arange(cube.header['NAXIS3'])-cube.header['CRPIX3']+1)*cube.header['CDELT3']+cube.header['CRVAL3']
	velocity = (velocity*u.m/u.s).to(unit).value
	#velocity = cube.spectral_axis.to(unit).value
	flag = (velocity >= min(mode)) & (velocity <= max(mode))
	for i in range(len(window)//2):
		subwin=window[i*2:i*2+2]
		flag[(velocity > min(subwin)) & (velocity < max(subwin))] = False
	if flag.sum()<3:
		print('Error - No enough range for baseline fitting.')
		return

	#calculate rms and output
	rms = np.squeeze(cube.data)[flag]
	nan = cube.data.max() if cube.header['BITPIX']>0 else -1000
	rms[rms == nan] = np.nan
	warnings.filterwarnings('ignore', r'Mean of empty slice')
	rms = np.nanmean((rms**2),axis=0)**0.5
	hdu = fits.PrimaryHDU(data=rms, header=cube.header)
	hdu.writeto(cellname+'_rms.fits', overwrite=True)


if __name__ == '__main__':
	#cuberms('/share/data/mwisp/G010+00/0190+040U.fits',[-10,10],window=[-5,-3,3,6,7])
	cuberms('test/0150+015U.fits','test/0150+015_info.txt')
