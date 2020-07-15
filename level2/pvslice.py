' mwisp module: pvslice '

__version__ = '2.1'
__author__ = 'Shaobo Zhang'

import os, math, time
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.coordinates as coord
import matplotlib.pyplot as plt
from scipy.interpolate import interpn

def gauss2d(x, y):
	#a demo gaussian kernel function
	return np.exp(-(x**2+y**2)/2/0.5778**2)

def hashtagmultiply(a,b):
	return np.dot(a[:,np.newaxis],b[np.newaxis,:])

def convolve(img, y, x, kfunction):
	#x,y must have the same size
	ishp = img.shape
	o = np.ndarray(x.size, dtype=img.dtype)

	gx, gy = np.meshgrid(np.arange(ishp[1]),np.arange(ishp[0]))
	for i in range(x.size):
		k = kfunction(gx-x.ravel()[i], gy-y.ravel()[i])
		k = k/np.nansum(k)
		o[i] = np.nansum(img * k)
	return o.reshape(x.shape)

def pvslice(fitsfile=None, path_coord=None, width=1, step=0.5, spline=False, pathgal=True, kernel=None):
	'''
	Extract position-velocity map from a data cube

	Notes
	-----
	Only accept fits in celestial and Galactic system with CDELT1=CDELT2

	Parameters
	----------
	fitsfile : str
		the file name of datacube.
	path_coord : str or array-like
		if path_coord is a string, it should be the file name of the path, in "x y" column format,
	    or it should be a list/tuple/ndarray in [[x1,x2,x3,...],[y1,y2,y3,...]] format.
	width : int, optional
		the width (in pixel unit) of the belt to be averaged, default is 1.
	step : float, optional
		resample the input slice path with interval of step (in pixel unit), default is 0.5.
	spline : bool, optional
		set to True to do spline interpolation to the input path, the path will become smooth spline. default is False.
	pathgal : bool, optional
		set to True if the given path_coord is in galactic system, default is True.
	kernel: function object, optional
		kernel=f where f=f(x,y), normalization is not needed.
		if not provided, the pixels on pvslice image are bilinear interpolation of the cube.
		when a 2-dimensional function is given, the pixels will be the convolution result of cube and kernel
			at the begining of this file, there is an example of gaussian kernel function.
			WARNING - using kernel keyword will be very time-consuming.

	Outputs
	-------
	"pvslice.fits": pv map
	"pvslice.path": the resampled path
	"pvslice.belt": the outline of the belt
	
	Versions
	--------
	Apr,26,2012,v1.0
	May,16,2012,v1.1
		move the reference pixel of position to the center of the axis
	Dec,01,2015,v1.2
		encounter CTYPE error while calling EXTAST in astrolib,
	    	change CTYPE from 'VELOCITY' to 'VEL', from 'POSITION' to 'POS'
	Jan,04,2016,v2.0
		rewrite the procedure to accept arbitary path with WIDTH
		correct the path resampling error
	Nov,02,2018,v2.1
		rewrite in python3

	Examples
	--------
	>>> pvslice('XXX.fits', [[l1,l2,l3],[b1,b2,b3]], width=5, step=0.2)
	>>> pvslice('XXX.fits', 'path.cat', pathgal=False, spline=True)
	>>> from mwisp.mosaic import gauss2d
	>>> pvslice('XXX.fits', 'path.cat', kernel=gauss2d)

	Todo list
	---------
	give position with emission more weight when averaging
	multi-kernel
	'''

	if (fitsfile == None) | (path_coord == None):
		print('Syntax - mwisp.pvslice(fitsfile, catalog, width=1, step=0.5, pathgal=True, spline=False, kernel=None)')
		print('Syntax - mwisp.pvslice(fitsfile, [a, d], width=1, step=0.5, pathgal=True, spline=False, kernel=None)')
		return

	if not os.path.exists(fitsfile):
		print('Error - fits file does not exist.')
		return
	
	if step <= 0:
		print('Error - step must be positive.')
		return

	if width <= 0:
		print('Error - width must be positive.')
		return

	#read fits
	hdu = fits.open(fitsfile)[0]
	hdu.data = np.squeeze(hdu.data)
	w = WCS(hdu.header)

	#read catalog
	if type(path_coord) == type(''):
		if not os.path.exists(path_coord):
			print('Error - catalog file does not exits!')
			return
		catalog = open(path_coord,'r')
		path_coord = []
		for line in catalog:
			path_coord.append(line.split())
		catalog.close()
	else:
		path_coord = np.array(path_coord).T

	#convert coordinate
	path_coord = coord.SkyCoord(path_coord, unit='deg', frame='galactic' if pathgal else 'icrs')
	x, y = path_coord.to_pixel(w, origin=0)

	#convert polyline to spline
	if spline & x.size >=4:
		from scipy.interpolate import splprep, splev
		tck, u = splprep([x, y], s=0)
		x, y = splev(np.linspace(0,1,len(x)*20), tck)
	elif x.size < 4:
		print('Warning - path nod < 4, skip SPLINE interpolation.')

	#resample path with interval of step
	nodnum = len(x)
	px = x[0]	#give px, py the starting point
	py = y[0]	#then px, py will be a point on the path
	nod = 1	#index of ending nod on the segment
	path = [[px, py]]
	while nod <= nodnum-1:
		rest = math.sqrt((x[nod]-px)**2+(y[nod]-py)**2)
		if rest >= step:	#the rest of current segment is enough for taking another step
			px += (x[nod]-px)/rest*step
			py += (y[nod]-py)/rest*step
			path.append([px,py])
		else:	#not enough, find from the next segment or the next next one...except the last nod
			if nod == nodnum-1:	#the last nod, extend a little
				last = math.sqrt((x[nod]-px)**2+(y[nod]-py)**2)
				if last > step/1e3:
					px+=(x[nod]-px)/last*step
					py+=(y[nod]-py)/last*step
					path.append([px, py])
				break
			while nod < nodnum-1:
				ptop2 = math.sqrt((x[nod+1]-px)**2 + (y[nod+1]-py)**2)
				if ptop2 < step:
					nod += 1
				else:
					p1top2 = math.sqrt((x[nod+1]-x[nod])**2 + (y[nod+1]-y[nod])**2)
					dx21 = x[nod+1] - x[nod]
					dy21 = y[nod+1] - y[nod]
					crossx = x[nod]*y[nod+1]*dy21 - x[nod+1]*y[nod]*dy21 + dx21**2*px + dx21*dy21*py
					crossx /= dx21**2 + dy21**2
					crossy = x[nod+1]*y[nod]*dx21 - x[nod]*y[nod+1]*dx21 + dx21*dy21*px + dy21**2*py
					crossy /= dx21**2 + dy21**2
					crosstonext = math.sqrt((step**2 -(px-crossx)**2 -(py-crossy)**2)) #>0?
					px = crossx + dx21/p1top2*crosstonext
					py = crossy + dy21/p1top2*crosstonext
					path.append([px,py])
					nod += 1
					break
	path = np.array(path)
	pathx = path[:,0]
	pathy = path[:,1]

	#export the path
	pathout = coord.SkyCoord(0,0,unit='deg').from_pixel(pathx,pathy,w)
	pathout = pathout.galactic if pathgal else pathout.icrs
	f=open('pvslice.path','w')
	for i in range(len(pathout)):
		f.write(pathout[i].to_string()+'\n')
	f.close()
	
	#expand path to belt
	stepnum = len(pathx)
	width = int(width)
	if width > 1:
		dy = np.roll(pathy, 1) - np.roll(pathy, -1)
		dy[0] = pathy[0]-pathy[1]
		dy[-1] = pathy[-2]-pathy[-1]
		dx = np.roll(pathx, 1) - np.roll(pathx, -1)
		dx[0] = pathx[0]-pathx[1]
		dx[-1] = pathx[-2]-pathx[-1]
		pa = np.arctan(dy/dx)+np.pi/2
		pa[dx < 0] += np.pi
		dx = np.cos(pa)
		dy = np.sin(pa)
		pathx = np.repeat(pathx[:,np.newaxis], width, axis=1) + hashtagmultiply(dx, np.arange(width)-(width-1)/2.)
		pathy = np.repeat(pathy[:,np.newaxis], width, axis=1) + hashtagmultiply(dy, np.arange(width)-(width-1)/2.)
	else:
		width=1

	#export path boundary of the belt
	pathout = coord.SkyCoord(0,0,unit='deg').from_pixel(pathx,pathy,w)
	pathout = pathout.galactic if pathgal else pathout.icrs
	f=open('pvslice.belt','w')
	for i in range(len(pathout)):
		f.write(pathout[i,0].to_string()+'\n')
	for i in range(len(pathout)-1,-1,-1):
		f.write(pathout[i,-1].to_string()+'\n')
	f.write(pathout[0,0].to_string())
	f.close()

	#extract value from cube
	channum = hdu.header['NAXIS3']
	slice = np.ndarray([stepnum, width, channum], dtype=hdu.data.dtype)
	if type(kernel) != type(convolve):
		for i in range(channum):
			slice[:,:,i] = interpn((range(hdu.header['NAXIS2']),range(hdu.header['NAXIS1'])), hdu.data[i,:,:],
				np.array([pathy.ravel(),pathx.ravel()]).T, bounds_error=False, fill_value=np.nan).reshape(pathx.shape)
	else:
		for i in range(channum):
			print('convolving: %d/%d' % (i+1,channum))
			slice[:,:,i] = convolve(hdu.data[i,:,:], pathy, pathx, kernel)
	slice = np.nanmean(slice, axis=1)

	#output slice and resampled path
	pvhdu = fits.PrimaryHDU(slice)
	pvhdu.header['CTYPE1'] = 'V'
	pvhdu.header['CRPIX1'] = hdu.header['CRPIX3']
	pvhdu.header['CRVAL1'] = hdu.header['CRVAL3']/1e3	#convert m/s to km/s
	pvhdu.header['CDELT1'] = hdu.header['CDELT3']/1e3
	pvhdu.header['CTYPE2'] = 'P'
	pvhdu.header['CRPIX2'] = 1
	pvhdu.header['CRVAL2'] = 0
	pvhdu.header['CDELT2'] = step*abs(hdu.header['CDELT1'])
	pvhdu.header['HISTORY'] = 'PV file: '+fitsfile
	pvhdu.header['HISTORY'] = 'PV path:'
	for i in range(len(path_coord)):
		pvhdu.header['HISTORY'] = path_coord[i].to_string()
	pvhdu.header['HISTORY'] = str('PV width: %d pixels' %  width)
	pvhdu.header['HISTORY'] = 'Position in Degree'
	pvhdu.header['HISTORY'] = 'Velocity in km/s'
	pvhdu.writeto('pvslice.fits', overwrite = True)

#pvslice('/Users/shaobo/Work/mwips/L935/mosaic_U_clipv.fits',[[84,85,85],[-1.2,-0.9,-1.2]],width=5,step=0.5,pathgal=True,spline=True,kernel=gauss2d)