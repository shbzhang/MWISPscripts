' mwisp module: mosaic '

__version__ = '1.9'
__author__ = 'Shaobo Zhang'

import os, math, time, glob
import numpy as np
from astropy.io import fits
from cuberms import cuberms

output_dtype=np.float32

def mosaic(*crange, sb='U', path=None, output='mosaic', silent=False, display=False, fillrms=False):
	'''
	Mosaic fits file of DLH survey
	Parameters
	----------
	crange : six floats indicate gl1,gl2,gb1,gb2,velo1,velo2
		Galactic coordinate range and velocity range (in km/s).
	sb : str
		ideband, for 12CO, sb='U'; for 13CO, sb='L'; for C18O, sb='L2'.
	path : str list, optional
		list of the directories containing datacube, rms, and info files, default path is './'.
	output : str, optional
		prefix of output file names, default name is 'mosaic'.
	silent : bool, optional
		whether to suppress message in the terminal. Default is False.
	display : bool, optional
		whether to display the progress on a window. Default is True.
	fillrms : bool, optional
		whether to fill the cell whose rms file is missing with calculations according to mode and window recorded in info file
	Notes
	-----
	All the input fits file (data and rms) must be the same as the grid define by template.	
	Versions
	--------
	Nov,09,2011,v1.0
	Nov,28,2011,v1.1
		add the fitspath and rmspath keywords.
	Jan,16,2012,v1.2
		fix a error when mosaic region with gl > 180
	Apr,12,2012,v1.3
		change the output file name.
		output a mask fits file for the mosaic.
	May,09,2012,v1.4
		complete the header of the result fits file.
		the mask file is now called 'coverage' file to fit its function.
		add the keywords "display" for users who do not support window.
	Dec,25,2013,v1.5
		accept rms file created by procedure cuberms.
	Oct,13,2015,v1.6
		remove keyword "fitspath" and "rmspath".
		new keyword "path" now accepts string array, procedure will
	    	automatically search datacube and rms file in these directories.
	Dec,09,2016,v1.7
		add history to indicate the completeness of the mosaic.
		add keyword 'output' for the names of the output files.
		add keyword 'silent' to suppress any message in terminal.
		fix spelling mistakes.
	May,17,2018,v1.8
		fix a minor error that certain input range does not match the output range.
	Nov,02,2018,v1.9
		rewrite in python3
	Oct,20,2020,v1.10
		revise for MWISP-II
		add keyword 'fillrms' to mosaic cubes whose rms files are missing.
	Examples
	--------
	>>> pa = ['./', '/share/data/mwisp/G020+00', '/share/data/mwisp/G030+00']
	>>> mwisp.mosaic(29, 32.5, -1, 1.2, -30, 30, sb='L', path=pa, output='Test', silent=True)
	'''

	#Syntax prompt
	if len(crange) != 6:
		print("Syntax - mosaic(l1, l2, b1, b2, v1, v2, sb='U', path='./', output='mosaic', display=False, silent=False, fillrms=False)")
		return

	#Parameter Checking: CRANGE
	l1, l2, b1, b2, v1, v2 = crange
	#modify gl in the range of -60 ~ 300 deg, since there is no observation at 300 deg for MWISP
	l1, l2 = l1 % 360, l2 % 360
	if l1 > 300:
		l1 -= 360
	if l2 > 300:
		l2 -= 360
	#sort range
	if l1 > l2:
		l2, l1 = l1, l2
	if b1 > b2:
		b2, b1 = b1, b2
	if v1 > v2:
		v2, v1 = v1, v2
	#return if out of MWISP range
	### Modify the range if you really need to mosaic data of your own ###
	if l1 < -11 or l2 > 251 or b1 < -13 or b2 > 13:
		print("Galactic coordinate is beyond the coverage of DLH survey!")
		print("l = 0 ~ 250 deg, b = -10 ~ 10 deg")
		return

	#Parameter Checking: SB
	sb = sb.upper()
	if sb not in ['U','L','L2']:
		print("Error - unknown sideband. Please specify 'U', 'L' or 'L2'.")
		return

	#Parameter Checking: PATH
	paths=[]
	if path is None:
		#no path specified, use default for server 119.78.210.193
		paths = glob.glob('/share/data/mwisp/R*') + glob.glob('/share/data/mwisp/G*')
		paths.append('/share/data/mwisp/infofiles')
	else:
		#list all possible paths
		if type(path) is str:
			path=[path]
		for apath in path:
			if type(apath) is str:
				paths += glob.glob(apath)

	#set display
	if display:
		import matplotlib.pyplot as plt
		plt.ion()

	#shape of output FITS
	nx = int(math.floor(l2*120) - math.ceil(l1*120) + 1)
	if nx==0:
		nx=1
	cx = int(math.floor(l2*120) + 1)
	ny = int(math.floor(b2*120) - math.ceil(b1*120) + 1)
	cy = int(-math.ceil(b1*120) + 1)	#1-base for FITS file

	#related FITS file range with value*10)
	fitsl1 = int(math.floor(l1*2)*5)
	fitsl2 = int(math.ceil(l2*2)*5)
	fitsb1 = int(math.floor(b1*2)*5)
	fitsb2 = int(math.ceil(b2*2)*5)
	fitsnx, fitsny = (fitsl2-fitsl1)//5+1, (fitsb2-fitsb1)//5+1
	fitsnum = fitsnx * fitsny

	time_start = time.time()
	if not silent:
		print('Mosaic begins:')
		print('GL from %f to %f' % (l1 % 360, l2 % 360))
		print('GB from %f to %f' % (b1, b2))
		print('V  from %f to %f' % (v1, v2))

	count, num = 0, 0
	for lx10 in range(fitsl1, fitsl2+1, 5):
		for bx10 in range(fitsb1, fitsb2+1, 5):

			count += 1
			cellname = _getcellname(lx10, bx10)

			#search for cube file
			observed = False
			for apath in paths:
				cubefile = os.path.join(apath, cellname+sb+'.fits')
				if os.path.exists(cubefile):
					observed = True
					break
			if not observed:
				if not silent:
					print('[%d/%d]%s%s - FITS file not found.' % (count, fitsnum, cellname, sb))
				continue

			#search for rms file
			haverms = False
			for apath in paths:
				rmsfile = os.path.join(apath, cellname+sb+'_rms.fits')
				if os.path.exists(rmsfile):
					haverms = True
					break
			if not haverms:
				if fillrms:
					#search for info file
					haveinfo = False
					for apath in paths:
						infofile = os.path.join(apath, cellname+'_info.txt')
						if os.path.exists(infofile):
							haveinfo = True
							break
					if not haveinfo:
						if not silent:
							print('[%d/%d]%s%s - RMS or INFO file not found.' % (count, fitsnum, cellname, sb))
						continue
					else:
						cuberms(cubefile, infofile, silent=True)
						rmsfile = cellname+sb+'_rms.fits'
						#if not silent:
						#	print('[%d/%d]RMS file for %s%s is generated according to the INFO file.' % (count, fitsnum, cellname, sb))
				else:
					if not silent:
						print('[%d/%d]%s%s - RMS file not found.' % (count, fitsnum, cellname, sb))
					continue

			#read cube file
			if not silent:
				print("[%d/%d]%s%s - Mosaic file '%s'" % (count, fitsnum, cellname, sb, cubefile))
			hdu = fits.open(cubefile)[0]
			if hdu.header['BITPIX'] > 0:
				nan = (hdu.data == np.nanmax(hdu.data))
			else:
				nan = ((hdu.data == -1000) | ~np.isfinite(hdu.data))
				#nan = ((hdu.data != -1000) & np.isfinite(hdu.data))
			hdu.data[nan] = 0
			#hdu.data *= nan

			#read rms file
			rmshdu = fits.open(rmsfile)[0]
			if rmshdu.header['BITPIX'] > 0:
				nan = rmshdu.data == np.nanmax(rmshdu.data)
			else:
				nan = (rmshdu.data == -1000) | ~np.isfinite(rmshdu.data)
			rmshdu.data[nan] = np.inf
			wei = np.squeeze(1/rmshdu.data**2)

			#clip data velocity
			flag, hdu = _clipv(hdu, v1, v2)
			if not flag:
				print('Error - Velocity is out of the channel range!')
				return

			num += 1
			#initiate a mosaic cube and header
			if num == 1:
				nv = hdu.header['NAXIS3']
				mdat = np.zeros([1, nv, ny, nx], dtype=output_dtype)
				mwei = np.zeros([ny, nx], dtype=output_dtype)
				mcov = np.zeros([ny, nx], dtype=np.uint8)
				mhdr = hdu.header
				mvran = _channel2velocity(mhdr, np.array([0,mhdr['NAXIS3']-1]))

			#check velocity inconsistency
			cstnc = mhdr['NAXIS3'] == hdu.header['NAXIS3']	#consistent channel number
			cstv1 = abs(_channel2velocity(hdu.header,0) - mvran[0]) < mhdr['CDELT3']/1e3	#consistent velocity range
			cstv2 = abs(_channel2velocity(hdu.header,hdu.header['NAXIS3']-1) - mvran[1]) < mhdr['CDELT3']/1e3
			if not (cstnc and cstv1 and cstv2):
				print('Error - Inconsistent velocity dimension!')
				return

			#mosaic dat to mdat
			for i in range(nv):
				hdu.data[0,i,:,:] *= wei
			fcx = cx-lx10*12-1
			fcy = cy+bx10*12-1
			fnx = hdu.header['NAXIS1']
			fny = hdu.header['NAXIS2']
			x1 = fcx-(fnx-1)//2
			x2 = fcx+(fnx-1)//2	#subscript on the mosaic array, within x1:x2
			y1 = fcy-(fny-1)//2
			y2 = fcy+(fny-1)//2
			if x2 < 0 or x1 > (nx-1) or y2 < 0 or y1 > (ny-1):	#fits is out of the mosaic range
				continue
			xcut1 = (x1 < 0)*(-x1)
			xcut2 = fnx-1-(x2 > nx-1)*(x2-nx+1)
			ycut1 = (y1 < 0)*(-y1)
			ycut2 = fny-1-(y2 > ny-1)*(y2-ny+1)
			hdu.data = hdu.data[:,:,ycut1:ycut2+1,xcut1:xcut2+1]
			wei = wei[ycut1:ycut2+1,xcut1:xcut2+1]
			mdat[:,:,max([y1,0]):min([y2,ny-1])+1,max([x1,0]):min([x2,nx-1])+1] += hdu.data
			mwei[max([y1,0]):min([y2,ny-1])+1,max([x1,0]):min([x2,nx-1])+1] += wei
			fnx = 61
			fny = 61
			x1 = fcx-(fnx-1)//2
			x2 = fcx+(fnx-1)//2
			y1 = fcy-(fny-1)//2
			y2 = fcy+(fny-1)//2
			if x2 < 0 or x1 > (nx-1) or y2 < 0 or y1 > (ny-1):
				continue
			xcut1 = (x1 < 0)*(-x1)
			xcut2 = fnx-1-(x2 > nx-1)*(x2-nx+1)
			ycut1 = (y1 < 0)*(-y1)
			ycut2 = fny-1-(y2 > ny-1)*(y2-ny+1)
			mcov[max([y1,0]):min([y2,ny-1])+1,max([x1,0]):min([x2,nx-1])+1] = 1
			if display:
				plt.imshow(mwei)
				plt.pause(0.01)

	if num == 0:
		print('Error - No data in the selected region!')
		return
	mwei[mwei == 0] = np.nan
	for i in range(nv):
		mdat[:,[i],:,:] /= mwei
	mhdr['NAXIS1'] = nx
	mhdr['NAXIS2'] = ny
	mhdr['CTYPE1'] = 'GLON-CAR'
	mhdr['CTYPE2'] = 'GLAT-CAR'
	mhdr['CRVAL1'] = 0.
	mhdr['CRVAL2'] = 0.
	if cx > 21601:
		mhdr['CRPIX1'] = cx-43200
	else:
		mhdr['CRPIX1'] = cx
	mhdr['CRPIX2'] = cy
	mhdr['CDELT1'] = -30/3600	#add for MWISP-II
	mhdr['CDELT2'] = 30/3600
	mhdr['BUNIT'] = 'K (T_MB)'
	mhdr.remove('OBJECT')
	mhdr.remove('GLAT')
	mhdr.remove('GLON')
	mhdr['HISTORY'] = 'MOSAIC: '+time.strftime("%b %d %Y %H:%M:%S", time.localtime())
	mhdr['HISTORY'] = str('L from %f to %f' % (l1, l2))
	mhdr['HISTORY'] = str('B from %f to %f' % (b1, b2))
	mhdr['HISTORY'] = str('V from %f to %f' % (v1, v2))
	mhdr['HISTORY'] = str('Mosaic Completeness: [%d/%d]' % (num, fitsnum))

	#output
	mhdu = fits.PrimaryHDU(mdat, mhdr)
	mhdu.writeto(output+'_'+sb+'.fits', overwrite = True)
	mhdu = fits.PrimaryHDU(mwei, mhdr)
	mhdu.writeto(output+'_'+sb+'_weight.fits', overwrite = True)
	mhdu = fits.PrimaryHDU(mcov, mhdr)
	mhdu.writeto(output+'_'+sb+'_coverage.fits', overwrite = True)
	mrms = 1/np.sqrt(mwei)
	#mrms[~np.isfinite(mrms)] = np.nan
	mhdu = fits.PrimaryHDU(mrms, mhdr)
	mhdu.writeto(output+'_'+sb+'_rms.fits', overwrite = True)
	time_end = time.time()
	if not silent:
		print('Successfully mosaic %d cubes.' % (num))
		print('Done in %f seconds.' % (time_end-time_start))

def _getcellname(gl, gb):
	gl = gl % 3600
	return str('%04d%+04d' % (gl, gb))

def _clipv(hdu, v1, v2):
	#clip a velocity range (v in km/s) from data according to its fits header
	nc = hdu.header['NAXIS3']
	cv = hdu.header['CRVAL3']
	cp = hdu.header['CRPIX3']
	cd = hdu.header['CDELT3']
	if cd < 0:
		hdu.data = hdu.data[:, ::-1, :, :]
		cd = -cd
		cp = nc-cp+1
	c1 = int(math.floor((v1*1e3 - cv) / cd + cp -1))	#convert faster than v2c
	c2 = int(math.ceil((v2*1e3 - cv) / cd + cp -1))
	if c1 > nc-1 or c2 < 0:
		return False, None
	if c1 < 0:
		c1 = 0
	if c2 > nc-1:
		c2 = nc-1
	hdu.data = hdu.data[:,c1:c2+1,:,:]
	hdu.header['NAXIS3'] = c2-c1+1
	hdu.header['CRPIX3'] = cp-c1
	hdu.header['CDELT3'] = cd
	return True, hdu

#c is 0-base; v in km/s
def _channel2velocity(hdr, c):
	v = (c - hdr['CRPIX3'] + 1) * hdr['CDELT3'] + hdr['CRVAL3']
	return v/1e3

def _velocity2channel(hdr, v):
	c = (v * 1e3 - hdr['CRVAL3']) / hdr['CDELT3'] + hdr['CRPIX3'] -1
	return c
'''
def test(a,*b,c=0,**d):
	#argument to argument
	#other argument to *argument
	#keyword to keyword
	#other keyword to **keyword
	print(a)
	print(b)
	print(c)
	print(d)
'''
if __name__ == '__main__':
	mosaic(14.75,15.75,0.75,1.75,-10,10, sb='U',path='test',fillrms=True)

