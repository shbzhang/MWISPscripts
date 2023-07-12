' mwisp module: mosaic '

__version__ = '2.0'
__author__ = 'Shaobo Zhang'

import os, time, glob, logging
import numpy as np
from astropy.io import fits
from math import floor, ceil

###the valid l/b range of mosaic, slightly larger than the MWISP project.
valid_l1 = -11
valid_l2 = 251
valid_b1 = -12
valid_b2 = 12
###output data type
output_dtype = np.float32

def mosaic(*crange, sb='U', path=None, prefix='', suffix='.fits', \
	undone=None, output='mosaic', weightcube=False, silent=False, display=False, ):
	'''
	Mosaic fits file of DLH survey

	Parameters
	----------
	crange : six floats includes gl1, gl2, gb1, gb2, velo1, velo2.
		which indicate galactic coordinate range (in deg) and velocity range (in km/s).
		If a "undone" is provide, these ranges are ignored.
	sb : str
		Sideband, for 12CO, sb='U'; for 13CO, sb='L'; for C18O, sb='L2'.
	path : str or list of strs, optional
		A directories or a List of the directories containing datacube and rms files.
		If there are multiply data in different directories, the path on the top of the list will be used.
		Default path is the directories of data on SASMA.
	prefix : str or tuple of two strs, optional
		Prefix of datacube and rms files.
		if a string, datacube and rms will both use it.
		if a tuple of two strings, datacube will use the first one, while rms will use the second.
	suffix : str or tuple of two strs, optional
		Suffix of datacube and rms files. Similar to "prefix".
	undone : str, optional
		Prefix of the set of mosaicked files, which may miss several cells in the last run.
		You can set undone to the prefix of a mosaiked datacube to patch cells on it without re-mosaic the whole one.
		"crange" will be ignored, and mosaic will continue to use ranges in the undone datacube.
		e.g. >>> mosaic(undone='mosaic', output='newmosaic')
	output : str, optional
		Prefix of output file names, default name is 'mosaic'.
	weightcube: bool, optional
		Whether to use a 3-dimensional weight cube. Default is False.
		Set to True if part of cube is masked as nan. Such as cleaned cubes in which velocities are chopped.
	silent : bool, optional
		Whether to suppress message in the terminal. Default is False.
	display : bool, optional
		Whether to display the progress on a window. Default is False.
		display=True is not recommended when running on remote servers.


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
	May,13,2021,v1.10
		add keyword "undone"
	Jan,12,2023,v2.0
		log file is added to the output set to help keyword "undone" to work better.
		correct the CDELT1 error when MWISP-2 data are involving.
		add keyword "prefix", "suffix" to mosaic cleaned datacubes.


	Examples
	--------
	>>> pa = ['./', '/share/data/mwisp/G020+00', '/share/data/mwisp/G030+00']
	>>> mwisp.mosaic(29, 32.5, -1, 1.2, -30, 30, sb='L', path=pa, output='Test', silent=True)
	'''

	###Syntax prompt
	if len(crange) != 6 and undone is None:
		print("Syntax - mosaic(l1, l2, b1, b2, v1, v2, sb='U', path='./', prefix='', suffix='.fits', \n\
			undone=None, output='mosaic', weightcube=False, display=False, silent=False)")
		return

	###log file config
	logging.basicConfig(filename=_getMosaicLogName(output, sb), filemode='w', \
		level=logging.INFO, format='%(asctime)s %(message)s')
	def _prompt(message):
		if not silent: print(message)
		logging.info(message)

	#prepare paths
	paths=[]
	if path is None:
		#no path specified, use default
		paths = glob.glob('/share/data/mwisp/R*') + glob.glob('/share/data/mwisp/G*')
	else:
		#list all possible paths
		if type(path) is str:
			path=[path]
		for apath in path:
			if type(apath) is str:
				paths += glob.glob(apath)

	#prepare prefix
	if type(prefix) is str: prefix = (prefix, prefix)
	if type(suffix) is str: suffix = (suffix, suffix)

	#set display
	if display:
		import matplotlib.pyplot as plt
		plt.ion()

	if undone is None:
		#standardize crange
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
		if l1 < -11 or l2 > 251 or b1 < -12 or b2 > 12:
			silent = False
			_prompt("Error - Galactic coordinate is out of the range of DLH survey!")
			_prompt("        l = 350 ~ 0 ~ 250 deg, b = -5 ~ 5 deg")
			return
		###get shape of output FITS
		mscNX = int(floor(l2*120) - ceil(l1*120) + 1)
		if mscNX==0: mscNX = 1
		mscCX = int(floor(l2*120) + 1)
		mscNY = int(floor(b2*120) - ceil(b1*120) + 1)
		mscCY = int(-ceil(b1*120) + 1)	#1-base for FITS file
		mscCells = []
	else:
		if os.path.abspath(undone) == os.path.abspath(output):
			_prompt("Error - output should be different from undone.")
			return
		try:
			mscHdr = fits.getheader(_getMosaicCubeName(undone, sb))
			mscNX = mscHdr['NAXIS1']
			mscNY = mscHdr['NAXIS2']
			mscNC = mscHdr['NAXIS3']
			mscCX = mscHdr['CRPIX1']
			mscCY = mscHdr['CRPIX2']
			l1 = _carx2l(mscHdr, mscNX-1)
			l2 = _carx2l(mscHdr, 0)
			b1 = _cary2b(mscHdr, 0)
			b2 = _cary2b(mscHdr, mscNY-1)
			v1 = _c2v(mscHdr, 0)
			v2 = _c2v(mscHdr, mscNC-1)

			mscHdr = fits.getheader(_getMosaicWeightName(undone, sb))
			if mscHdr['NAXIS']>2: weightcube=True

			mscLog = open(_getMosaicLogName(undone, sb), 'r')
			mscCells = mscLog.readlines()[-1]
			mscCells = mscCells.split('Done Cells:')[-1]
			mscCells = mscCells.split()
		except:
			silent = False
			_prompt('Error - Unable to load undone files!')
			_prompt('        Make sure datacube/weight/coverage/log files are all ready.')
			return


	#related FITS file range with value*10)
	cellL1x10 = int(floor(l1*2)*5)
	cellL2x10 = int(ceil(l2*2)*5)
	cellsLx10 = range(cellL1x10, cellL2x10+1, 5)
	cellB1x10 = int(floor(b1*2)*5)
	cellB2x10 = int(ceil(b2*2)*5)
	cellsBx10 = range(cellB1x10, cellB2x10+1, 5)
	cellNum = len(cellsLx10) * len(cellsBx10)

	time_start = time.time()
	_prompt('Mosaic begins:')
	_prompt('GL from %f to %f' % (l1 % 360, l2 % 360))
	_prompt('GB from %f to %f' % (b1, b2))
	_prompt('V  from %f to %f' % (v1, v2))

	index, valid = 0, 0
	for cellLx10 in cellsLx10:
		for cellBx10 in cellsBx10:

			index += 1
			cellName = _getCellName(cellLx10, cellBx10)

			if undone is not None:
				if cellName in mscCells:
					_prompt('[%d/%d]%s%s: mosaic is already in the undone.' % (index, cellNum, cellName, sb))
					continue

			###Prepare files
			#search for cube file
			observed = False
			for apath in paths:
				cellCubeFile = _getCellCubeName(apath, prefix, cellName, sb, suffix)
				if os.path.exists(cellCubeFile):
					observed = True
					break
			if not observed:
				_prompt('[%d/%d]%s%s: datacube file not found.' % (index, cellNum, cellName, sb))
				continue

			#search for rms file
			haveRms = False
			for apath in paths:
				cellRmsFile = _getCellRmsName(apath, prefix, cellName, sb, suffix)
				if os.path.exists(cellRmsFile):
					haveRms = True
					break
			if not haveRms:
				_prompt('[%d/%d]%s%s: RMS file not found.' % (index, cellNum, cellName, sb))
				continue


			###Load data
			#read cube file
			_prompt("[%d/%d]%s%s: mosaic using: '%s'" % (index, cellNum, cellName, sb, cellCubeFile))
			try:
				cellCubeHDU = fits.open(cellCubeFile)[0]
				cellHdr = cellCubeHDU.header
				cellDat = cellCubeHDU.data

				while cellDat.ndim<4: cellDat = cellDat[np.newaxis]

				###Clip velocity
				nc = cellHdr['NAXIS3']
				cv = cellHdr['CRVAL3']
				cp = cellHdr['CRPIX3']
				cd = cellHdr['CDELT3']
				#flip if necessary
				if cd < 0:
					cellDat = np.flip(cellDat, axis=-3)
					cd = -cd
					cp = nc-cp+1
				#regulate channel range
				c1 = int(floor((v1*1e3 - cv) / cd + cp -1))	#convert faster than v2c
				c2 = int(ceil((v2*1e3 - cv) / cd + cp -1))
				if c1 > nc-1 or c2 < 0: #out of range
					_prompt('Error - Velocity from %f to %s is out of the range!' % (v1, v2))
					return
				if c1 < 0: c1 = 0
				if c2 > nc-1: c2 = nc-1
				#cellDat = cellDat[:, c1:c2+1, :, :]
				cellDat = np.take(cellDat, range(c1,c2+1), axis=-3)
				cellHdr['NAXIS3'] = c2-c1+1
				cellHdr['CRPIX3'] = cp-c1
				cellHdr['CDELT3'] = cd

				###deal with NAN
				if cellHdr['BITPIX'] > 0:
					nanDat = cellDat == np.nanmax(cellDat)
				else:
					nanDat = (cellDat == -1000) | ~np.isfinite(cellDat)
				cellDat[nanDat] = 0
			except:
				_prompt('Error - Unable to load datacube file: %s' % (index, cellNum, cellName, sb, cellCubeFile))
				continue


			#read rms file
			try:
				cellRmsHDU = fits.open(cellRmsFile)[0]
				cellRms = np.squeeze(cellRmsHDU.data)

				###deal with NAN
				if cellRmsHDU.header['BITPIX'] > 0:
					nanRms = cellRms == np.nanmax(cellRms)
				else:
					nanRms = (cellRms == -1000) | ~np.isfinite(cellRms)
				cellRms[nanRms] = np.inf
				cellWei = 1/cellRms**2

				if weightcube:
					cellWei = np.zeros_like(cellDat)+cellWei
					cellWei[nanDat] = 0
			except:
				_prompt('Error - Unable to load rms file: %s' % (index, cellNum, cellRmsFile))
				continue



			###Check SHAPE !!!



			#cube and weight are ready
			valid += 1

			#initiate a mosaic cube and header
			if valid == 1:
				if undone is None:
					###init an empty one
					mscNC = cellHdr['NAXIS3']
					mscDat = np.zeros((cellDat.shape[0], mscNC, mscNY, mscNX), dtype=output_dtype)
					if weightcube:
						mscWei = mscDat.copy()
					else:
						mscWei = np.zeros((mscNY, mscNX), dtype=output_dtype)
					mscCov = np.zeros((mscNY, mscNX), dtype=np.uint8)
					mscHdr = cellHdr.copy()
				else:
					###load a mosaicked undone
					mscHDU = fits.open(_getMosaicCubeName(undone, sb))[0]
					mscHdr = mscHDU.header
					mscNC = mscHdr['NAXIS3']
					mscDat = mscHDU.data
					mscDat[np.isnan(mscDat)] = 0
					mscWei = fits.open(_getMosaicWeightName(undone, sb))[0].data
					mscWei[np.isnan(mscWei)] = 0
					mscDat *= mscWei
					mscCov = fits.open(_getMosaicCoverageName(undone, sb))[0].data
				mscV1 = _c2v(mscHdr, 0)
				mscV2 = _c2v(mscHdr, mscHdr['NAXIS3']-1)
			else:
				#check velocity consistency
				sameNC = mscHdr['NAXIS3'] == cellHdr['NAXIS3']	#consistent channel number
				sameV1 = abs(_c2v(cellHdr,0) - mscV1) < mscHdr['CDELT3']/1e3	#consistent velocity range
				sameV2 = abs(_c2v(cellHdr,cellHdr['NAXIS3']-1) - mscV2) < mscHdr['CDELT3']/1e3
				if not (sameNC and sameV1 and sameV2):
					_prompt('Error - Inconsistent velocity dimension!')
					continue


			#mosaic dat to mscDat
			cellDat *= cellWei

			###clip cube/weight
			cellCX = mscCX-cellLx10*12-1
			cellCY = mscCY+cellBx10*12-1
			cellNX = cellHdr['NAXIS1']
			cellNY = cellHdr['NAXIS2']
			cellX1 = cellCX-(cellNX-1)//2
			cellX2 = cellCX+(cellNX-1)//2	#subscript on the mosaic array, within cellX1:cellX2
			cellY1 = cellCY-(cellNY-1)//2
			cellY2 = cellCY+(cellNY-1)//2
			if cellX2 < 0 or cellX1 > (mscNX-1) or cellY2 < 0 or cellY1 > (mscNY-1): continue	#fits is out of the mosaic range
			#cut cell if on boundary of msc
			cutX1 = -cellX1 if cellX1<0 else 0
			cutX2 = cellNX-cellX2+mscNX-2 if cellX2>mscNX-1 else cellNX-1
			cutY1 = -cellY1 if cellY1<0 else 0
			cutY2 = cellNY-cellY2+mscNY-2 if cellY2>mscNY-1 else cellNY-1
			#position of cell on msc
			posX1 = max(cellX1, 0)
			posX2 = min(cellX2, mscNX-1)
			posY1 = max(cellY1, 0)
			posY2 = min(cellY2, mscNY-1)
			cellDat = cellDat[..., cutY1:cutY2+1, cutX1:cutX2+1]
			mscDat[..., posY1:posY2+1, posX1:posX2+1] += cellDat
			if weightcube:
				cellWei = cellWei[..., cutY1:cutY2+1, cutX1:cutX2+1]
				mscWei[..., posY1:posY2+1, posX1:posX2+1] += cellWei
			else:
				cellWei = cellWei[cutY1:cutY2+1, cutX1:cutX2+1]
				mscWei[posY1:posY2+1, posX1:posX2+1] += cellWei
			
			cellNX = 61
			cellNY = 61
			cellX1 = cellCX-(cellNX-1)//2
			cellX2 = cellCX+(cellNX-1)//2
			cellY1 = cellCY-(cellNY-1)//2
			cellY2 = cellCY+(cellNY-1)//2
			if cellX2 < 0 or cellX1 > (mscNX-1) or cellY2 < 0 or cellY1 > (mscNY-1): continue
			posX1 = max(cellX1, 0)
			posX2 = min(cellX2, mscNX-1)
			posY1 = max(cellY1, 0)
			posY2 = min(cellY2, mscNY-1)
			mscCov[posY1:posY2+1, posX1:posX2+1] = 1

			if display:
				plt.imshow(mscWei)
				plt.pause(0.01)

			mscCells.append(cellName)

	if valid == 0:
		_prompt('Error - No file in the selected region!')
		return

	_prompt('Writing to files...')
	mscWei[mscWei == 0] = np.nan
	mscDat /= mscWei
	if weightcube: mscWei = np.nanmax(mscWei, axis=-3)

	mscHdr['NAXIS1'] = mscNX
	mscHdr['NAXIS2'] = mscNY
	mscHdr['CTYPE1'] = 'GLON-CAR'
	mscHdr['CTYPE2'] = 'GLAT-CAR'
	mscHdr['CRVAL1'] = 0.
	mscHdr['CRVAL2'] = 0.
	mscHdr['CDELT1'] = -30/3600
	mscHdr['CDELT2'] = 30/3600
	if 'CTYPE4' not in mscHdr:
		mscHdr['CTYPE4'] = '            '
		mscHdr['CRVAL4'] =  1
		mscHdr['CDELT4'] =  1
		mscHdr['CRPIX4'] =  1
		mscHdr['CROTA4'] =  0
	if mscCX > 21601:
		mscHdr['CRPIX1'] = mscCX-43200
	else:
		mscHdr['CRPIX1'] = mscCX
	mscHdr['CRPIX2'] = mscCY
	mscHdr['BUNIT'] = 'K (T_MB)'

	if undone is None:
		mscHdr.remove('OBJECT')
		mscHdr.remove('GLAT')
		mscHdr.remove('GLON')
		mscHdr['HISTORY'] = 'MOSAIC: '+time.strftime("%b %d %Y %H:%M:%S", time.localtime())
	else:
		mscHdr['HISTORY'] = 'MOSAIC PATCH: '+time.strftime("%b %d %Y %H:%M:%S", time.localtime())
	mscHdr['HISTORY'] = str('L from %f to %f' % (l1, l2))
	mscHdr['HISTORY'] = str('B from %f to %f' % (b1, b2))
	mscHdr['HISTORY'] = str('V from %f to %f' % (v1, v2))
	mscHdr['HISTORY'] = str('Mosaic Completeness: [%d/%d]' % (len(mscCells), cellNum))

	#output
	mscHDU = fits.PrimaryHDU(mscDat, mscHdr)
	mscHDU.writeto(_getMosaicCubeName(output, sb), overwrite = True)
	mscHDU = fits.PrimaryHDU(mscWei, mscHdr)
	mscHDU.writeto(_getMosaicWeightName(output, sb), overwrite = True)
	mscHDU = fits.PrimaryHDU(mscCov, mscHdr)
	mscHDU.writeto(_getMosaicCoverageName(output, sb), overwrite = True)
	mscRms = 1/np.sqrt(mscWei)
	#mscRms[~np.isfinite(mscRms)] = np.nan
	mscHDU = fits.PrimaryHDU(mscRms, mscHdr)
	mscHDU.writeto(_getMosaicRmsName(output, sb), overwrite = True)
	time_end = time.time()
	_prompt('Mosaic successfully.')
	_prompt('Done in %f seconds.' % (time_end-time_start))
	_prompt('Done Cells: '+' '.join(mscCells))


def _getCellName(gl, gb):
	gl = gl % 3600
	return '%04d%+04d' % (gl, gb)

def _getCellCubeName(path, prefix, cellName, sb, suffix):
	return os.path.join(path, '%s%s%s%s' % (prefix[0], cellName, sb, suffix[0]))
def _getCellRmsName(path, prefix, cellName, sb, suffix):
	return os.path.join(path, '%s%s%s_rms%s' % (prefix[1], cellName, sb, suffix[1]))

def _getMosaicLogName(output, sb):
	return '%s_%s.log' % (output, sb)
def _getMosaicCubeName(output, sb):
	return '%s_%s.fits' % (output, sb)
def _getMosaicWeightName(output, sb):
	return '%s_%s_weight.fits' % (output, sb)
def _getMosaicCoverageName(output, sb):
	return '%s_%s_coverage.fits' % (output, sb)
def _getMosaicRmsName(output, sb):
	return '%s_%s_rms.fits' % (output, sb)

def _carx2l(hdr, x):
	return (x - hdr['CRPIX1'] + 1) * hdr['CDELT1'] + hdr['CRVAL1']
def _cary2b(hdr, y):
	return (y - hdr['CRPIX2'] + 1) * hdr['CDELT2'] + hdr['CRVAL2']

#c is 0-base; v in km/s
def _c2v(hdr, c):
	v = (c - hdr['CRPIX3'] + 1) * hdr['CDELT3'] + hdr['CRVAL3']
	return v/1e3
def _v2c(hdr, v):
	c = (v * 1e3 - hdr['CRVAL3']) / hdr['CDELT3'] + hdr['CRPIX3'] -1
	return c

if __name__ == '__main__':
	#mosaic()
	#mosaic(95.5, 97, 4.5, 6, -10, 10, path='/Users/shaobo/Downloads/', output='mosaic', weightcube=True)
	#mosaic(95.5, 97, 4.5, 6, -120, 120, path='/Users/shaobo/Downloads/', output='clip1', weightcube=True, prefix=('CLIP',''))
	#mosaic(path='/Users/shaobo/Downloads/', output='clip3', undone='clip1', weightcube=True, prefix=('CLIP',''))

	#mosaic(214,217,-1,2,-300,300, path='/share/public/qzyanShare/mwispIfinal/longBaseClean/all_cubes/', output='test1', prefix='NC')
	mosaic(214,217,-1,2,-300,300, path='/share/public/qzyanShare/mwispIfinal/longBaseClean/CO13/all_cubes', output='test1', prefix='NC', sb='L')

