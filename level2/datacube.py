' mwisp module: datacube '

__version__ = '1.1'
__author__ = 'Shaobo Zhang'

'''
A class for MWISP datacube reduction
subclass of spectral_cube.SpectralCube
add new methods to mask, pvslice, resample (more efficient), and plot the datacube.

Usage
>>> from datacube import DataCube
>>> cube = DataCube.openMWISP('cubename.fits')
'''

###use spectral
import copy, warnings
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

import radio_beam
from astropy.io import fits
from astropy.units.quantity import Quantity
from spectral_cube import SpectralCube, OneDSpectrum
from spectral_cube.lower_dimensional_structures import Projection
from astropy.convolution import Gaussian1DKernel, Box1DKernel
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from reproject import reproject_interp
import astropy.units as u
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def open(file):
	return DataCube.open(file)

def openMWISP(file):
	return DataCube.openMWISP(file)

class DataCube(SpectralCube):
	'''
	Creating with:
	>>> DataCube(data, wcs, header=, rms=, **kws)
	***check more keywords in SpectralCube

	Reading with:
	>>> cube = DataCube.open(filename)
	or
	>>> cube = DataCube.open(HDU)
	or to open a MWISP cube
	>>> cube = DataCube.openMWISP(filename)
	***check more keywords in SpectralCube.read


	New property:
	rms


	New method:
	open
	openMWISP
	copy
	with_rms
	x_mask
	y_mask
	z_mask
	v_mask
	pvslice
	average
	channel
	velocity
	with_window
	basline
	get_rms
	rebinvelocity
	resample
	channelmap
	gridmap
	peakvelocity
	huevelocity


	General manupulations with inherited method:

	Subcube (cf. SpectralCube.subcube)
	use slicing:
	>>> newcube = cube[:, 10:25, 14:32]
	or use subcube:
	>>> newcube = cube.subcube(xlo=84*u.deg, xhi=86*u.deg, ylo=-2*u.deg, yhi=1*u.deg, zlo=-5*u.km/u.s, zhi=5*u.km/u.s) 

	Spatial Reprojection (cf. SpectralCube.reproject)
	>>> newcube = cube.reproject(referenceheader, order='bilinear', use_memmap=False, filled=True)

	Spatial Smooth (cf. SpectralCube.convolve_to)
	if cube has no beam info:
	>>> beam = radio_beam.Beam(major=52*u.arcsec, minor=52*u.arcsec, pa=0*u.deg)
	>>> cube = cube.with_beam(beam)
	then smoothing:
	>>> newbeam = radio_beam.Beam(major=4*u.arcmin, minor=4*u.arcmin, pa=0*u.deg)
	>>> newcube = cube.with_beam(beam).convolve_to(newbeam)
	
	Cube Moment (cf. SpectralCube.moment)
	>>> m0 = cube.moment(order=0)
	or a goodlooking moment map
	>>> cube = cube.with_rms(rms)
	>>> mask = cube.v_mask(userms=True, threshold=3, consecutive=3)
	>>> m0 = cube.with_mask(mask).moment(order=0)

	LVmap/BVmap (cf. SpectralCube.moment)
	>>> lv = cube.moment(order=0, axis=2)
	>>> bv = cube.moment(order=0, axis=1)

	RGB 3 Components
	rgb3c = cube.rebinvelocity(-7,7,3,largecube=True)
	ax1 = plt.subplot(projection=rgb3c.wcs, slices=('x', 'y', 3))
	ax1.imshow(rgb3c._data[::-1].transpose(1,2,0)/np.nanmax(rgb3c)*2)
	'''

	def __init__(self, *args, rms = None, **kws):
		super().__init__(*args, **kws)
		self._rms = DataCube._load_rms(rms, self.header)


	def __getitem__(self, index):
		new = super().__getitem__(index)
		if self._rms is not None:
			new._rms = self._rms[index[1:]]
		return new


	def _view(sc):
		###view a SpectralCube as DataCube
		return DataCube(sc._data, sc._wcs, mask=sc._mask, meta=sc._meta, fill_value=sc._fill_value, header=sc._header, \
			allow_huge_operations=sc.allow_huge_operations, beam=sc._beam, wcs_tolerance=sc._wcs_tolerance)


	def read(*args, **kws):
		'''
		Read an ordinary FITS file or HDU.
		'''
		sc = SpectralCube.read(*args, **kws)
		return DataCube._view(sc)


	def open(*args, **kws):
		return DataCube.read(*args, **kws)


	def openMWISP(file):
		'''
		Open a MWISP FITS file or HDU
		'''
		hdu = fits.open(file)[0]
		hdu.header['BUNIT'] = 'K'
		return DataCube.read(hdu)


	def copy(self):
		'''
		deepcopy a DataCube
		'''
		return copy.deepcopy(self)


	def _load_rms(obj, hdr):
		###load array from a filename/HDU/array
		if obj is None: return None
		elif isinstance(obj, str): arr = fits.open(obj)[0].data
		elif isinstance(obj, fits.PrimaryHDU): arr = obj.data
		elif isinstance(obj, np.ndarray): arr = obj
		else: raise TypeError('rms must be a FITS filename (str), a HDU (fits.PrimaryHDU), or an array (numpy.ndarray).')
		arr = np.squeeze(arr)
		###check if the arr is broadcastable to data
		wcs = WCS(hdr, naxis=2)
		if arr.ndim == 2 and arr.shape[-1]==hdr['NAXIS1'] and arr.shape[-2]==hdr['NAXIS2']:
			return Projection(arr, wcs=wcs)
		else: raise ValueError('rms shape is not broadcastable to data shape: %s vs %s.' % (str(arr.shape), str(wcs.array_shape)))


	def with_rms(self, rms):
		'''
		return a new DataCube with RMS

		Parameters
		----------
		rms: array, str, or PrimaryHDU
			use as rms of the datacube.
			if a str, load rms from the fits file;
			if a HDU, load data as rms;

		Returns
		-------
		DataCube with rms

		Versions
		--------
		Jul,07,2023,v2.2

		Examples
		--------
		>>> #load rms from a fits file
		>>> cube.with_rms('rms.fits')
		>>> #load array as rms to the datacube
		>>> cube.with_rms(np.ones(100,150))
		'''
		new = DataCube._view(self)
		new._rms = DataCube._load_rms(rms, self.header)
		return new


	@property
	def rms(self):
		return self._rms
	

	#  __  __                 _      _                 
	# |  \/  |   __ _   ___  | | __ (_)  _ __     __ _ 
	# | |\/| |  / _` | / __| | |/ / | | | '_ \   / _` |
	# | |  | | | (_| | \__ \ |   <  | | | | | | | (_| |
	# |_|  |_|  \__,_| |___/ |_|\_\ |_| |_| |_|  \__, |
	#                                            |___/ 
	'''
	X_MASK

	Y_MASK

	Z_MASK

	V_MASK
	'''

	def x_mask(self, xfunc):
		'''
		return a mask according to f(y,z).
		x values above f(y,z) are marked as True

		Parameters
		----------
		xfunc: Quantity or function
			x above which are marked as True

		Returns
		-------
		mask array

		Versions
		--------
		Jul,07,2023,v2.2

		Examples
		--------
		>>> llo = lambda b,v: b+86*u.deg
		>>> mask = cube.x_mask(llo)
		>>> maskedCube = cube.with_mask(mask)
		'''
		zg, yg, xg = self.world[:]
		if callable(xfunc): return xg > xfunc(yg, zg)
		else: return xg > xfunc


	def y_mask(self, yfunc):
		'''
		return a mask according to f(x,z).
		y values above f(x,z) are marked as True

		Parameters
		----------
		yfunc: Quantity or function
			y above which are marked as True

		Returns
		-------
		mask array

		Versions
		--------
		Jul,07,2023,v2.2

		Examples
		--------
		>>> blo = lambda l,v: l-86*u.deg
		>>> mask = cube.y_mask(blo)
		>>> maskedCube = cube.with_mask(mask)
		'''
		zg, yg, xg = self.world[:]
		if callable(yfunc): return yg > yfunc(xg, zg)
		else: return yg > yfunc


	def z_mask(self, zfunc):
		'''
		return a mask according to f(x,y).
		z values above f(x,y) are marked as True

		Parameters
		----------
		zfunc: Quantity or function
			z above which are marked as True

		Returns
		-------
		mask array

		Versions
		--------
		Jul,07,2023,v2.2

		Examples
		--------
		>>> #use different vrange for moment
		>>> vflo = lambda l,b: ((l/u.deg-84.5)**2*1.5-12.5)*u.km/u.s
		>>> vfhi = 10*u.km/u.s
		>>> mask = cube.z_mask(vflo) & ~cube.z_mask(vfhi)
		>>> maskedCube = cube.with_mask(mask)
		'''
		zg, yg, xg = self.world[:]
		if callable(zfunc): return zg > zfunc(xg, yg)
		else: return zg > zfunc


	def v_mask(self, userms=True, threshold=3, consecutive=3):
		'''
		return a mask including "consecutive" chanels above "threshold".

		Parameters
		----------
		userms: bool, optional
			if set to True, use 'rms * threshold' as threshold
		threshold: float, optional
			values above which are marked as True
		consecutive: int, optional
			minimum number of consecutive channel

		Returns
		-------
		mask array

		Versions
		--------
		Jul,07,2023,v2.2

		Examples
		--------
		>>> cube = cube.with_rms(rms)
		>>> mask = cube.v_mask(userms = True, threshold = 5, consecutive=5)
		>>> cube = cube.with_mask(mask)
		'''
		if userms:
			if self.rms is None: print('DataCube.v_mask - RMS not specified, only threshold will be used.')
			else: threshold = self.rms * threshold
		mask = self._data > threshold
		for i in range(consecutive-1):
			mask = mask & np.roll(mask, +1, axis=0)
		mask[:consecutive] = False
		for i in range(consecutive-1):
			mask = mask | np.roll(mask, -1, axis=0)
		return mask


	def grow_mask(mask, step=1, axis=0):
		'''
		return a mask grown from input mask

		Parameters
		----------
		step: int, optional
			step to grow the boundary
		axis: int or tuple, optional
			axis or axes to grow

		Return
		------
		the grown mask

		Versions
		--------
		Jul,07,2023,v2.2

		Examples
		--------
		>>> largeMask = DataCube.grow_mask(cube.mask.include(), step=2, axis=(0,1,2))
		'''
		if isinstance(axis, int): axis=(axis,)
		for s in range(step):
			for ax in axis:
				mask = mask | np.roll(mask,-1, axis=ax) | np.roll(mask, +1, axis=ax)
		return mask


	#  ____    ____       _      _____   ___      _      _     
	# / ___|  |  _ \     / \    |_   _| |_ _|    / \    | |    
	# \___ \  | |_) |   / _ \     | |    | |    / _ \   | |    
	#  ___) | |  __/   / ___ \    | |    | |   / ___ \  | |___ 
	# |____/  |_|     /_/   \_\   |_|   |___| /_/   \_\ |_____|	
	'''
	PVSLICE

	AVERAGE
	'''

	def pvslice(self, path_coord, path_frame='galactic', keep_width=False,
		width=1, step=0.5, spline=False, kernel=None):
		'''
		Extract position-velocity map from a data cube

		Notes
		-----
		Only accept datacube in celestial and Galactic system with CDELT1=CDELT2

		Parameters
		----------
		path_coord: SkyCoord, str or array-like
			if path_coord is a str, it should be the file name of the path, in "x y" column format,
			or it should be a SkyCoord object or a list/tuple/ndarray in [[x1,x2,x3,...],[y1,y2,y3,...]] format.
		path_frame: str, optional
			the coordinate frame of path_coord if path_coord is not a SkyCoord object
		keep_width: bool, optional
			whether to average the width axis in the slice
			if set to True, return a cube of (width x path x velcity).
		width: int, optional
			the width (in pixel unit) of the belt to be averaged, default is 1.
		step: float, optional
			resample the input slice path with interval of step (in pixel unit), default is 0.5.
		spline: bool, optional
			set to True to do spline interpolation to the input path, the path will become smooth spline. default is False.
		kernel: function object, optional
			kernel=f(x,y), normalization is not needed.
			if not provided, the pixels on pvslice image are bilinear interpolation of the cube.
			when a 2-dimensional function is given, the pixels will be the convolution result of cube and kernel.
			***using kernel keyword will be very time-consuming.

		Returns
		-------
		a HDUList with pvslice map in the primary HDU, and path in the second HDU.

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
		Jul,07,2023,v2.2
			modify to a method in DataCube class

		Examples
		--------
		>>> #slicing along the spline of path
		>>> path = [[83.5,-2.5], [84.5,-2], [85.2,0.5], [86.5,1.5]]
		>>> pvhdulist = cube.with_spectral_unit('km/s').pvslice(path, step=1, width=5, spline=True)
		>>> pvmap, pvpath = pvhdulist

		>>> #output map/path to file
		>>> pvhdulist.writeto('pvslice.fits', overwrite=True)

		>>> #show the slicing path
		>>> m0 = cube.moment(order=0)
		>>> ax = plt.subplot(projection=m0.wcs)
		>>> ax.imshow(m0.data)
		>>> ax.plot(pvpath.data[...,0], pvpath.data[...,1], '.', transform=ax.get_transform('galactic'))
		>>> plt.show()

		>>> #show the pvslice map
		>>> ax = plt.subplot(projection=WCS(pvmap.header))
		>>> ax.imshow(pvmap.data)
		>>> ax.coords[0].set_format_unit('km/s')
		>>> ax.coords[0].set_major_formatter('x')
		>>> plt.show()

		pvslice with kernel
		>>> #define a gauss2d kernel
		>>> g2d = lambda dx,dy: np.exp(-(dx**2+dy**2)/2/2**2)
		>>> pv, path = cube.pvslice(path, step=2, width=5, kernel=g2d)

		Todo list
		---------
		give position with emission more weight when averaging
		multi-kernel
		'''
		if step <= 0: raise ValueError('step must be positive float')
		if width <= 0: raise ValueError('width must be positive int.')
		wcs = self.wcs.sub(2)

		###load path coord
		if not isinstance(path_coord, SkyCoord):
			if isinstance(path_coord, str):
				path_coord = np.loadtxt(path_coord)
			try:
				path_coord = SkyCoord(path_coord, frame=path_frame, unit='deg')
			except:
				raise ValueError('path_coord must be a SkyCoord or an array-like object of shape (N, 2).')
		###coord to pixel
		x, y = path_coord.to_pixel(wcs, origin=0)

		#convert polyline to spline
		if spline and (x.size >=4):
			from scipy.interpolate import splprep, splev
			tck, u = splprep([x, y], s=0)
			x, y = splev(np.linspace(0,1,len(x)*20), tck)
		elif x.size < 4:
			warnings.warn('path nodes < 4, skip SPLINE interpolation.')

		#resample path with interval of step
		nodnum = len(x)
		px = x[0]	#give px, py the starting point
		py = y[0]	#then px, py will be a point on the path
		node = 1	#index of ending node on the segment
		path = [[px, py]]
		while node <= nodnum-1:
			rest = np.sqrt((x[node]-px)**2+(y[node]-py)**2)
			if rest >= step:	#the rest of current segment is enough for taking another step
				px += (x[node]-px)/rest*step
				py += (y[node]-py)/rest*step
				path.append([px,py])
			else:	#not enough, find from the next segment or the next next one...except the last node
				if node == nodnum-1:	#the last node, extend a little
					last = np.sqrt((x[node]-px)**2+(y[node]-py)**2)
					if last > step/1e3:
						px+=(x[node]-px)/last*step
						py+=(y[node]-py)/last*step
						path.append([px, py])
					break
				while node < nodnum-1:
					ptop2 = np.sqrt((x[node+1]-px)**2 + (y[node+1]-py)**2)
					if ptop2 < step:
						node += 1
					else:
						p1top2 = np.sqrt((x[node+1]-x[node])**2 + (y[node+1]-y[node])**2)
						dx21 = x[node+1] - x[node]
						dy21 = y[node+1] - y[node]
						crossx = x[node]*y[node+1]*dy21 - x[node+1]*y[node]*dy21 + dx21**2*px + dx21*dy21*py
						crossx /= dx21**2 + dy21**2
						crossy = x[node+1]*y[node]*dx21 - x[node]*y[node+1]*dx21 + dx21*dy21*px + dy21**2*py
						crossy /= dx21**2 + dy21**2
						crosstonext = np.sqrt((step**2 -(px-crossx)**2 -(py-crossy)**2)) #>0?
						px = crossx + dx21/p1top2*crosstonext
						py = crossy + dy21/p1top2*crosstonext
						path.append([px,py])
						node += 1
						break
		path = np.array(path)
		pathx = path[:,0]
		pathy = path[:,1]

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
			pathx = pathx + (np.arange(width)-(width-1)/2)[:,np.newaxis] * dx
			pathy = pathy + (np.arange(width)-(width-1)/2)[:,np.newaxis] * dy
		else:
			width=1
			pathx = pathx[np.newaxis]
			pathy = pathy[np.newaxis]

		#extract value from cube
		channum = self.header['NAXIS3']
		pvdata = np.ndarray([width, stepnum, channum], dtype=self._data.dtype)
		if callable(kernel):
			###convolve with kernel
			gx, gy = np.meshgrid(np.arange(self.shape[2]), np.arange(self.shape[1]))
			dx = gx[:,:,np.newaxis,np.newaxis] - pathx
			dy = gy[:,:,np.newaxis,np.newaxis] - pathy
			k = kernel(dx, dy)
			k = k / np.nansum(k, axis=(0,1))
			for i in tqdm(range(channum)):
				pvdata[:,:,i] = np.nansum(self._data[i,:,:,np.newaxis,np.newaxis] * k, axis=(0,1))
		else:
			###simple interpolation
			from scipy.interpolate import interpn
			points = (np.arange(self.header['NAXIS2']), np.arange(self.header['NAXIS1']))
			for i in tqdm(range(channum)):
				pvdata[:,:,i] = interpn(points, self._data[i], (pathy, pathx), bounds_error=False, fill_value=np.nan)
		#average width axis
		if not keep_width:
			pvdata = np.nanmean(pvdata, axis=0)

		#output pvdata and resampled path
		pvmap = fits.PrimaryHDU(pvdata)
		pvmap.header['CTYPE1'] = 'VELO'
		pvmap.header['CRPIX1'] = self.header['CRPIX3']
		pvmap.header['CRVAL1'] = self.header['CRVAL3']
		pvmap.header['CDELT1'] = self.header['CDELT3']
		pvmap.header['CUNIT1'] = self._spectral_unit.to_string()
		pvmap.header['CTYPE2'] = 'P'
		pvmap.header['CRPIX2'] = 1
		pvmap.header['CRVAL2'] = 0
		pvmap.header['CDELT2'] = step*abs(self.header['CDELT1'])
		pvmap.header['CUNIT2'] = 'deg'
		if keep_width:
			pvmap.header['CTYPE3'] = 'WIDTH'
			pvmap.header['CRPIX3'] = (width-1)/2
			pvmap.header['CRVAL3'] = 0
			pvmap.header['CDELT3'] = abs(self.header['CDELT1'])
			pvmap.header['CUNIT3'] = 'deg'
		pvmap.header['HISTORY'] = 'PV path:'
		for i in range(len(path_coord)):
			pvmap.header['HISTORY'] = path_coord[i].to_string()
		pvmap.header['HISTORY'] = str('PV width: %d pixels' %  width)

		pvpath = SkyCoord(0,0,unit='deg').from_pixel(pathx, pathy, wcs)
		pvpath = np.dstack((pvpath.data.lon.value, pvpath.data.lat.value))
		pvpath = fits.ImageHDU(pvpath)
		pvpath.header['CTYPE1'] = 'WIDTH'
		pvpath.header['CTYPE2'] = 'STEPS'
		pvpath.header['CTYPE3'] = 'LON/LAT'

		pvlist = fits.HDUList([pvmap, pvpath])
		return pvlist


	def average(self):
		'''
		Average spectra in the DataCube and return a spectrum

		Return
		----------
		a OneDSpectrum

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> equWeiAverSpec = cube.average()
		>>> rmsWeiAverSpec = cube.with_rms(rms).average()
		>>> plt.step(cube.velocity(), equWeiAverSpec._data, label='Equal weighting')
		>>> plt.step(cube.velocity(), rmsWeiAverSpec._data, label='RMS weighting')
		>>> plt.legend()
		>>> plt.show()
		'''
		if self._rms is None:
			print('DataCube.average - RMS not specified, use equal weight.')
			return self.mean(axis=(1,2))
		else:
			print('DataCube.average - average with RMS weighting.')
			return (self * self.rms._data).mean(axis=(1,2)) / np.nanmean(self.rms._data)

	#  ____                          _                    _ 
	# / ___|   _ __     ___    ___  | |_   _ __    __ _  | |
	# \___ \  | '_ \   / _ \  / __| | __| | '__|  / _` | | |
	#  ___) | | |_) | |  __/ | (__  | |_  | |    | (_| | | |
	# |____/  | .__/   \___|  \___|  \__| |_|     \__,_| |_|
	#         |_|                                           
	'''
	CHANNEL

	VELOCITY

	WITH_WINDOW

	BASLINE

	GET_RMS

	REBINVELOCITY

	RESAMPLE
	'''

	def channel(self, velocity=None, vunit='km/s'):
		'''
		Return 0-base channel axis, or convert velocity to channel

		Parameters
		----------
		velocity: array or Quantity, optional
			velcity to be converted to channel.
			if not specified, return the channel axis
		vunit: str or Unit, optional
			unit of velocity
			ignored if velocity is Quantity

		Return
		----------
		channel
			array of channel axis, or converted channel

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> #Return 0-base channel axis of a DataCube
		>>> caxis = cube.channel()
		>>> #Convert velcity to channel
		>>> channel = cube.channel(velocity = [-15, 20], vunit = 'km/s')
		'''
		if velocity is None: channel = np.arange(self.header['NAXIS3'])
		else:
			if isinstance(velocity, Quantity): velocity = velocity.to(self._spectral_unit).value
			else: velocity = np.array(velocity) * u.Unit(vunit).to(self._spectral_unit)
			channel = (velocity - self.header['CRVAL3']) / self.header['CDELT3'] + self.header['CRPIX3'] - 1
		return channel


	def velocity(self, channel=None, vunit='km/s'):
		'''
		Return velocity axis, or convert channel to velocity.

		Parameters
		----------
		channel: array, optional
			channel to be converted to velocity.
			if not specified, return the velocity axis.
			channel is in 0-base.
		vunit: str or Unit, optional
			unit of output velocity.

		Return
		----------
		channel
			array of channel axis, or converted channel

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> #Return velocity axis of a DataCube
		>>> vaxis = cube.velocity()
		>>> #Convert channel to velcity
		>>> velocity = cube.velocity(channel = [30, 220], vunit = 'km/s')
		'''
		wcsc = self.wcs[2]
		if channel is None: channel = np.arange(self.header['NAXIS3'])
		else: channel = np.array(channel)
		velocity = (channel - self.header['CRPIX3'] + 1) * self.header['CDELT3'] + self.header['CRVAL3']
		return (velocity*self._spectral_unit).to(vunit)


	def with_window(self, *window, modex=None, vunit='km/s'):
		'''
		Return a DataCube with mask including the baseline

		Parameters
		----------
		window: floats, functions, optional
			velocity ranges ignored in baseline fitting.
		modex: array of shape (2,), optional
			velocity range included in baseline fitting.
		vunit: str or Unit, optional
			unit of window and modex.

		Return
		----------
		a DataCube with mask

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> maskedCube = cube.with_window(-60, -40, -20, 20, modex = [-200, 200], vunit='km/s')
		*see DataCube.baseline for more infomation.
		'''
		zg, yg, xg = self.world[:]
		###set mode x
		if modex is None: mask = np.ones(self.shape, dtype=bool)
		elif len(modex) == 2:
			lo, hi = modex
			if callable(lo): lo = lo(xg, yg)
			elif not isinstance(lo, Quantity): lo = lo*u.Unit(vunit)
			if callable(hi): hi = hi(xg, yg)
			elif not isinstance(hi, Quantity): hi = hi*u.Unit(vunit)
			mask = (zg > lo) & (zg < hi)
		else: raise ValueError('modex must has two elements')
		###set window
		if window is None: return self.with_mask(mask)
		elif (len(window) % 2) == 0:
			for i in range(0, len(window)-1, 2):
				lo = window[i]
				hi = window[i+1]
				if callable(lo): lo = lo(xg, yg)
				elif not isinstance(lo, Quantity): lo = lo*u.Unit(vunit)
				if callable(hi): hi = hi(xg, yg)
				elif not isinstance(hi, Quantity): hi = hi*u.Unit(vunit)
				mask[(zg > lo) & (zg < hi)] = False
			return self.with_mask(mask)
		else: raise ValueError('number of window must be odd.')


	def baseline(self, deg=0):
		'''
		Do a polynomial baseline fitting

		Parameters
		----------
		deg: int, optional
			the degree of polynomial fitting

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> #set mode and window
		>>> fittedCube = cube.with_window(-60, -40, -20, 20, modex = [-200, 200], vunit='km/s')
		>>> #do the baseline fitting
		>>> fittedCube.baseline(deg=1)
		>>> #calculate RMS
		>>> rms = fittedCube.rms()
		'''
		if deg == 0:
			p0 = self.mean(axis=0)._data
			self._data -= p0
		else:
			channel = self.channel()
			p = np.ndarray((deg+1, self.shape[1], self.shape[2]))
			it = np.nditer(self._data[0], flags=['multi_index'])
			while not it.finished:
				spectrum = self._data[:, it.multi_index[0], it.multi_index[1]]
				mask = self.mask.include()[:, it.multi_index[0], it.multi_index[1]] & np.isfinite(spectrum)
				if mask.sum()>5:
					p[:, it.multi_index[0], it.multi_index[1]] = np.polyfit(channel[mask], spectrum[mask], deg)
				it.iternext()
			for d in range(deg+1):
				self._data -= p[-d-1,:,:]*self.channel()[:,np.newaxis,np.newaxis]**d


	def get_rms(self):
		'''
		return a RMS map in current mask

		Return
		----------
		a fits HDU like map

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> cube.with_window(-30, 30, modex=(-100,100), vunit='km/s')
		>>> rmsMap = cube.get_rms()
		'''
		return ((self**2).mean(axis=0))**0.5


	def rebinvelocity(self, start, stop, num, vunit='km/s', largecube=True):
		'''
		Rebin the velocity axis of a DataCube

		Parameters
		----------
		start, stop: float
			the center velocity of first and last channel in the rebinned cube.
		num: int
			the total number of channels.
		vunit: str or Unit, optional
			unit of start and stop.
		largecube: bool, optional
			if cube is small, use method of SpectralCube
			if cube is large, the fraction-integration method is much quicker.

		Return
		----------
		the binned DataCube

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> RebinnedCube = cube.rebinvelocity(-10, 10, 21, vunit='km/s', largecube=True)
		'''
		if isinstance(start, Quantity): start = start.to(vunit).value
		if isinstance(stop, Quantity): stop = stop.to(vunit).value
		if num < 2: raise ValueError('channel number "num" must be greater than 1.')

		if not largecube:
			step = (stop - start) / (num - 1)
			new_axis = np.linspace(start, stop, num) * u.Unit(vunit)
			kernel = Box1DKernel(step / self.header['CDELT3'] / self._spectral_unit.to(vunit))
			return self.spectral_smooth(kernel).spectral_interpolate(new_axis)
		else:
			#velocity to channel
			chn_start, chn_stop = self.channel(velocity=[start, stop], vunit=vunit)
			#channels and boundary
			step = (chn_stop - chn_start) / (num - 1)
			chn_bound = np.linspace(chn_start-step/2, chn_stop+step/2, num+1)

			#new datacube header
			newheader = self.header.copy()
			newheader['NAXIS3'] = num
			newheader['CRVAL3'] = start
			newheader['CRPIX3'] = 1	#fits is 1-base
			newheader['CDELT3'] = (stop-start)/(num-1)
			newheader['CUNIT3'] = vunit

			def _IntegralArray(ndarray, idxrange, axis=0):
				#integral an ndarray within range along an axis	
				step = 1 if idxrange[1]>idxrange[0] else -1
				intrange = np.round(idxrange).astype(np.int)
				fltrange = idxrange - (intrange-0.5)
				
				head = np.take(ndarray, intrange[0], axis=axis) * fltrange[0]
				rear = np.take(ndarray, intrange[1], axis=axis) * fltrange[1]
				intarray = np.nansum(np.take(ndarray, range(intrange[0], intrange[1], step), axis = axis), axis = axis)
				integral = intarray + rear - head
				return integral*step

			#rebin
			dim = (num, self.shape[1], self.shape[2])
			newdata = np.ndarray(dim, dtype = self._data.dtype)
			newmask = np.ndarray(dim, dtype = bool)
			for c in range(num):
				if ((chn_bound[c:c+2] < -0.5) | (chn_bound[c:c+2] > self.shape[0]-0.5)).any():
					newdata[c] = np.nan
					newmask[c] = False
				else:
					newdata[c] = _IntegralArray(self._data, chn_bound[c:c+2], axis=0) / step
					newmask[c] = _IntegralArray(self.mask.include(), chn_bound[c:c+2], axis=0) / step >0.5

			newhdu = fits.PrimaryHDU(data=newdata, header=newheader)
			return DataCube.read(newhdu).with_mask(newmask)


	def resample(self, *Xpar, vunit='km/s', largecube=False):
		'''
		Resample velocity of a DataCube

		Parameters
		----------
		Xpar: Header or 4 floats
			if a cube header is provided, load resample parameter from it.
			else provide 4 numbers as Nx, Xref, Xval, Xinc
		vunit: str or Unit, optional
			unit of start and stop.
		largecube: bool, optional
			cf. this keyword in DataCube.rebinvelocity

		Return
		----------
		the resampled DataCube

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> #resample to align with header
		>>> ResampleDataCube = cube.resample(referenceheader, vunit='km/s')
		>>> #resample with values
		>>> ResampleDataCube = cube.resample(201, 100, 0.0, 1.0, vunit = 'km/s')
		'''
		if isinstance(Xpar[0], fits.Header):
			Nx, Xref, Xval, Xinc = Xpar[0]['NAXIS3'], Xpar[0]['CRPIX3']-1, Xpar[0]['CRVAL3'], Xpar[0]['CDELT3']
			if 'CUNIT3' in Xpar[0] and Xpar[0]['CUNIT3'] != '': vunit = Xpar[0]['CUNIT3']
		elif len(Xpar) >= 4: Nx, Xref, Xval, Xinc = Xpar[:4]
		else: raise TypeError('Input should be a FITS header (fits.Header), or 4 numbers indicate Nx, Xref, Xval, Xinc.')

		start = (0 - Xref) * Xinc + Xval
		stop = (Nx-1 - Xref) * Xinc + Xval
		num = Nx
		self.rebinvelocity(start, stop, num, vunit=vunit, largecube=largecube)


	#  ____    _        ___    _____ 
	# |  _ \  | |      / _ \  |_   _|
	# | |_) | | |     | | | |   | |  
	# |  __/  | |___  | |_| |   | |  
	# |_|     |_____|  \___/    |_|  
	'''
	CHANNELMAP

	GRIDMAP

	PEAKVELOCITY

	HUEVELOCITY
	'''

	def channelmap(self, nrows = None, ncols = None, figureheight=6, vunit = 'km/s', subplot_kws = None, \
		imshow_kws = None, contour_kws = None, text_kws = None, colorbar_kws = None, tick_kws = None):
		'''
		Plot channel map of a DataCube

		Parameters
		----------
		nrows, ncols : int, optional
			number of panels in row, column.
			panels will be arranged automatically if both are not specified.
		figureheight: float, optional
			the height of figure, change nrows/ncols to adjust figure size.
		vunit: str, optional
			unit of velocity to show in each panel.
		subplot_kws: dict, bool, optional
			keywords to adjust the subplots.
			accept left, right, bottom, top, hspace, wspace, and other keyword for subplots_adjust.
		imshow_kws: dict, bool, optional
			keywords to plot the background image.
			accept vmin, vmax, cmap, and other keyword for imshow.
			set to False to remove them.
		contour_kws: dict, bool, optional
			keywords to plot the contours.
			accept levels, colors, linewidths, and other keyword for contour
			set to False to remove them.
		text_kws: dict, bool, optional
			keywords to label the velocity.
			accept fontsize, ha, va, and other keyword for text
			set to False to remove them.
		colorbar_kws: dict, bool, optional
			keywords to plot the colorbar.
			accept keywords for colorbar
			set to False to remove it.
		tick_kws: dict, bool, optional
			keywords of ticks on axis.
			accept size, direction, and other keyword for tick_params

		Return
		----------
		fig, axes
			figure and axes in grid shape.

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> #rebin velocity first
		>>> cube = cube.rebinvelocity(-10, 10, 21, vunit='km/s')
		>>> fig, ax = cube.channelmap(nrows=3, ncols=7, vunit='km/s', figureheight=8, \
			imshow_kws = dict(vmin=-0.1, vmax=5, cmap='rainbow'),\
			contour_kws = dict(levels=[1,2,3,4]),\
			text_kws = dict(x=0.1, y=0.1),\
			tick_kws = dict(size=10, direction='in'),\
			colorbar_kws = None)
		
		>>> #ax[-1,0] is the lower left pannel with tick labels
		>>> ax[-1,0].set_xlabel('glon')
		>>> ax[-1,0].set_ylabel('glat')
		>>> ax[-1,0].coords[0].set_major_formatter('d.d')
		>>> ax[-1,0].coords[1].set_major_formatter('d.d')
		>>> plt.show()
		'''
		#get number of rows and columns
		def _AutoArrange(npanels, panel_ratio, page_ratio=1.0):
			#decide nrows and ncols for several panels of certain aspect ratio
			#ratio = width / height
			nrows = np.arange(1,50)
			ncols = np.ceil(npanels/nrows)
			ratio = ncols/nrows*panel_ratio
			best = np.abs(ratio-page_ratio).argmin()
			return nrows[best].astype(np.int), ncols[best].astype(np.int)

		vnaxis = self.shape[0]
		if nrows is None and ncols is None:
			nrows, ncols = _AutoArrange(vnaxis, self.header['NAXIS1']/self.header['NAXIS2'], page_ratio=4/3)
		elif nrows is None:
			nrows = np.ceil(vnaxis / ncols).astype(np.int)
		elif ncols is None:
			ncols = np.ceil(vnaxis / nrows).astype(np.int)
		elif nrows*ncols < vnaxis:
				warnings.warn('nrows x ncols is less than the number of channels, extra channels will be ignored')
				vnaxis = nrows*ncols

		#create figure
		page_ratio = (self.header['NAXIS1']*ncols) / (self.header['NAXIS2']*nrows)

		fig = plt.figure(figsize=(figureheight*page_ratio, figureheight))

		subplot_params = dict(left=0.05, right=0.95, bottom=0.05, top=0.95, wspace=0, hspace=0)
		if subplot_kws: subplot_params.update(subplot_kws)
		fig.subplots_adjust(**subplot_params)

		#add subplot axes
		proj = self.wcs.sub(2)
		axes = [plt.subplot(nrows, ncols, 1, projection=proj)]
		for i in range(2, vnaxis+1):
			axes.append(plt.subplot(nrows, ncols, i, projection=proj, sharex=axes[0], sharey=axes[0]))

		#plot parameters
		vmin = np.nanmin(self._data)
		vmax = np.nanmax(self._data)
		imshow_params = dict(origin='lower', vmin=vmin, vmax=vmax, cmap=plt.get_cmap('Greys'))
		if imshow_kws: imshow_params.update(imshow_kws)
		contour_params = dict(levels=np.arange(0.1, 1, 0.2)*vmax, colors='black', linewidths=1)
		if contour_kws: contour_params.update(contour_kws)
		text_params = dict(x=0.02, y=0.98, fontsize=10, horizontalalignment='left', verticalalignment='top')
		if text_kws: text_params.update(text_kws)
		vtext = self.velocity(vunit = vunit)
		tick_params = dict(size=10, direction = 'in')
		if tick_kws: tick_params.update(tick_kws)

		#plot panels
		#self._fill_value = 0
		for panel, ax in enumerate(axes):
			img = self.filled_data[panel].value
			if imshow_kws != False:
				im = ax.imshow(img, **imshow_params)  # Display the image slice
			if contour_kws != False:
				ax.contour(img, **contour_params)
			#draw velocity TEXT
			if vnaxis > 1 and text_kws != False:
				ax.text(s=vtext[panel].to_string(precision=1), transform=ax.transAxes, **text_params)
			#settings
			ax.tick_params(which='both', axis='both', **tick_params)
			ax.coords[0].display_minor_ticks(True)
			ax.coords[1].display_minor_ticks(True)
			#hide ticklabels
			if panel != ncols*(nrows-1):
				ax.coords[0].set_ticklabel_visible(False)
				ax.coords[1].set_ticklabel_visible(False)

		if (imshow_kws != False) & (colorbar_kws != False):
			axins = inset_axes(axes[-1],
				width="5%",  # width = 10% of parent_bbox width
				height="90%" if vnaxis > 1 else '100%',  # height : 50%
				loc='lower left',
				bbox_to_anchor=(1.05, 0., 1, 1),
				bbox_transform=axes[-1].transAxes,
				borderpad=0)
			colorbar_params = dict(cax=axins)
			if colorbar_kws:
				colorbar_params.update(colorbar_kws)
			#add colorbar
			cb = fig.colorbar(im, **colorbar_params)

		#return fig, axes array
		while len(axes)<nrows*ncols: axes.append(None)
		axes = np.array(axes).reshape(nrows, ncols)
		return fig, axes


	def gridmap(self, modey=None, unit='K'):
		'''
		show spectra in grid as class> map /g

		Parameters
		----------
		modey: int, optional
			the y limit of panels
		unit: str, astropy.units.Unit, optional
			unit of modey

		Return
		----------
		fig, axes
			figure and axes in grid shape.

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> #extract required region
		>>> cube = cube[:, 10:20, 15:25]
		>>> fig, ax = cube.gridmap()
		>>> plt.show()
		'''
		dim = self.shape
		fig, ax = plt.subplots(nrows=dim[1], ncols=dim[2], figsize=[6*dim[-1]/dim[-2],6])
		fig.subplots_adjust(wspace=0, hspace=0)
		vaxis = self.velocity()

		if modey is None: modey = np.array([self._data.min(), self._data.max()])
		elif not isinstance(modey, Quantity): modey = np.array(modey)*u.Unit(unit).to(self.unit)

		for ix in range(dim[2]):
			for iy in range(dim[1]):
				ax[iy,ix].step(vaxis, self._data[..., iy, ix])
				ax[iy,ix].set_ylim(modey)
				if ix != 0 or iy != dim[1]-1:
					ax[iy,ix].set_xticklabels('')
					ax[iy,ix].set_yticklabels('')
		ax[-1,0].set_xlabel('$v_{lsr}$ (%s)' % (self._spectral_unit.to_string()))
		ax[-1,0].set_ylabel('T (%s)' % (self.unit.to_string()))
		return fig, ax


	def peakvelocity(self):
		'''
		return a peakvelocity map, where each pixel represents the velocity of peak along the line of sight.

		Return
		----------
		a fits HDU like map

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> #masking noise
		>>> mask = cube.v_mask()
		>>> peakv = cube.with_mask(mask).peakvelocity()
		>>> #show the map
		>>> ax = plt.subplot(projection=peakv.wcs)
		>>> ax.imshow(peakv.data)
		>>> plt.show()
		'''
		peakc = self.argmax(axis=0)
		peakv = self.velocity(peakc)
		allnan = self.mask.exclude().all(axis=0)
		peakv[allnan] = np.nan
		wcs = self.wcs.sub(2)
		header = self.header.copy()
		return Projection(peakv, wcs=self.wcs.sub(2), header=header)


	def huevelocity(self):
		'''
		return a N*M*3 map with each velocity cooresponding to a Hue value.

		Return
		----------
		a RGB image and wcs

		Versions
		--------
		Jul,07,2023,v1.0

		Examples
		--------
		>>> #get a huevelocity map
		>>> img, wcs = cube.huevelocity()
		>>> #show the map
		>>> ax = plt.subplot(projection=wcs)
		>>> ax.imshow(img)
		>>> plt.show()
		'''
		import colorsys
		rgbshape = (self.shape[1], self.shape[2], 3)
		rgbdata = np.zeros(rgbshape)
		nc = self.shape[0]
		step = 1 if self.header['CDELT3']>0 else -1
		for c, channelmap in enumerate(self._data[::step]):
			hue = (1-c/(nc-1)) * 0.666667	##0.6667 scales hue from red to blue but not back to red
			r,g,b = colorsys.hls_to_rgb(hue, 0.5, 0.9)
			rgbdata[...,0] += channelmap*r
			rgbdata[...,1] += channelmap*g
			rgbdata[...,2] += channelmap*b
		rgbdata /= np.nanmax(rgbdata)
		rgbwcs = self.wcs.sub(2)
		return rgbdata, rgbwcs



if __name__ == '__main__':
	#init
	old = 'example/old.fits'
	ref = 'example/ref.fits'
	rms = 'example/rms.fits'
	cxy = 'example/cxy.fits'
	ycx = 'example/ycx.fits'
	cube = DataCube.openMWISP(old)
	print(cube,'\n')
	#print(help(cube))

	###test rms
	if 0:
		print('new and old are using the same data')
		new = cube.with_rms(rms)
		new._data[0,0,0] = 42
		print('data[0,0,0] in oldcube is ', cube._data[0,0,0])
		new = cube.with_rms(None)
		print('1.RMS in newcube is ', new.rms)
		new = cube.with_rms(rms)
		print('2.RMS in oldcube is ', cube.rms)
		print('  RMS in newcube is ', new.rms.shape)
		new = cube.with_rms(fits.open(rms)[0])
		print('3.RMS in newcube is ', new.rms.shape)
		new = cube.with_rms(np.ones(cube.shape[1:]))
		print('4.RMS in newcube is ', new.rms.shape)
		new = cube.with_rms(np.ones((3,4)))

	###test slice rms
	if 0:
		#newcube = cube[0:10,10:25,20:45]
		new = cube.with_rms(rms)
		sub = new.subcube(xlo=85*u.deg, xhi=85.5*u.deg, ylo=-0.5*u.deg, yhi=0*u.deg, zlo=-2*u.km/u.s, zhi=2*u.km/u.s)
		print('cube with rms:', new)
		print('rms shape:', new.rms.shape)
		print('clip cube with rms:', sub)
		print('rms shape:', sub.rms.shape)

	###test c/v convert
	if 0:
		new = cube[20:30]
		print('Channels are ', new.channel())
		print(new.channel([-10, 10], 'km/s'))
		print(new.channel(np.array([10,-10])*u.Unit('km/s')))
		print('Velocities are ', new.velocity())
		print(new.velocity([53,179]))
		print(new.velocity([53,179], vunit='m/s'))

	###test average
	if 0:
		new = cube.with_rms(None)
		osp = new.average()
		new = cube.with_rms(rms)
		nsp = new.average()
		plt.step(new.velocity(), osp._data, label='Equal weighting')
		plt.step(new.velocity(), nsp._data, label='RMS weighting')
		plt.legend()
		plt.show()

	###test lv/bv map
	if 0:
		new = cube.with_spectral_unit(u.km/u.s).with_mask(cube>0*u.K)
		lv = new.subcube(ylo=-0.5*u.deg, yhi=0*u.deg).moment(order=0, axis=1)
		bv = new.subcube(xlo=85*u.deg, xhi=85.5*u.deg).moment(order=0, axis=2)
		fig = plt.figure(figsize=(10, 8))
		ax = fig.add_subplot(121, projection=lv.wcs, slices=('y', 'x'))
		im = ax.imshow(lv.T._data)
		ax.coords[1].set_format_unit('km/s')
		ax.coords[1].set_major_formatter('x.x')
		ax.set_xlabel('v$_{lsr}$')
		ax.set_ylabel('glon')

		ax = fig.add_subplot(122, projection=bv.wcs, slices=('y', 'x'))
		im = ax.imshow(bv.T._data)
		ax.coords[1].set_format_unit('km/s')
		ax.coords[1].set_major_formatter('x.x')
		ax.set_xlabel('v$_{lsr}$')
		ax.set_ylabel('glat')
		plt.show()

	###test mask
	if 0:
		lf = lambda b,v: -b+85*u.deg
		cube.with_mask(~cube.y_mask(lf)).moment(order=0).quicklook()
		bf = lambda l,v: l-86*u.deg
		cube.with_mask(~cube.y_mask(bf)).moment(order=0).quicklook()

		vflo = lambda l,b: ((l/u.deg-84.5)**2*1.5-12.5)*u.km/u.s
		vfhi = 10*u.km/u.s
		mask = cube.z_mask(vflo) & ~cube.z_mask(vfhi)
		cube.with_mask(mask).moment(order=0, axis=0).quicklook()
		plt.show()

		###show v_mask
		cube.moment(order=0).quicklook()
		cube.with_rms(rms).with_mask(cube.v_mask(threshold=5, consecutive=5)).moment(order=0).quicklook()
		plt.show()

	#test baseline
	if 0:
		new = cube[:,100:110, 120: 135].with_window(-1, 13, modex=[-19, 19])

		plt.step(new.velocity(), new._data[:,5,-2], 'k-')
		plt.step(new.velocity(), new.mask.include()[:,5,-2])
		
		new.baseline(deg=0)
		plt.step(new.velocity(), new._data[:,5,-2], 'r-')
		plt.plot(new.velocity(), np.zeros(new.velocity().size))
		plt.show()

		new.get_rms().quicklook()
		plt.show()

	#test spectral smooth
	if 0:
		from astropy.convolution import Gaussian1DKernel, Box1DKernel
		new = cube[:,100:110, 120:135]

		reb = new.rebinvelocity(-15,15,31)
		reh = new.resample(reb.header)
		res = new.resample(31, 15, 0.0, 1.0)
		ref = new.rebinvelocity(-15,15,31, largecube=True)

		plt.plot(new.velocity(), new._data[:,5,13], 'k.-')
		plt.plot(reb.velocity(), reb._data[:,5,13], 'r.-')
		plt.plot(reh.velocity(), reh._data[:,5,13]+0.1, 'g.-')
		plt.plot(res.velocity(), res._data[:,5,13]+0.2, 'b.-')
		plt.plot(ref.velocity(), ref._data[:,5,13]+0.3, 'y.-')
		plt.show()

	###test gridmap
	if 0:
		new = cube[:, 100:110, 120:132]
		fig, ax = new.gridmap()
		plt.show()

	###test channel map
	if 0:
		new = cube[:, 100:170, 100:150].rebinvelocity(-10,10,21, largecube=True)
		new._data[:,50:] = np.nan
		fig, ax = new.channelmap(nrows=3, ncols=7, figureheight=8, \
			imshow_kws=dict(cmap='gist_rainbow_r'))#,vmax=new.filled_data[:].max().value/3))
		ax[-1,0].set_xlabel('glon')
		ax[-1,0].set_ylabel('glat')
		ax[-1,0].coords[0].set_major_formatter('d.d')
		ax[-1,0].coords[1].set_major_formatter('d.d')
		plt.show()

	###test hue velocity map
	if 0:
		rgb3c = cube.rebinvelocity(-7,7,3,largecube=True)
		ax1 = plt.subplot(121, projection=rgb3c.wcs, slices=('x', 'y', 3))
		ax1.imshow(rgb3c._data[::-1].transpose(1,2,0)/np.nanmax(rgb3c)*2)
		ax1.set_title('3 components of rgb')
		rgbhue, wcs = cube.subcube(zlo=-10.6*u.km/u.s, zhi=10.6*u.km/u.s).huevelocity()
		ax2 = plt.subplot(122, projection=wcs)
		ax2.imshow(rgbhue*2)
		ax2.set_title('velocity as hue')
		plt.show()

	###test peakvelocity
	if 0:
		new = cube.subcube(zlo=-10.5*u.km/u.s, zhi=10.5*u.km/u.s)
		peakv = new.with_mask(new.v_mask()).peakvelocity()
		print(peakv.header)
		ax = plt.subplot(projection=peakv.wcs)
		im = ax.imshow(peakv.data, cmap='gist_rainbow_r')
		plt.colorbar(im, ax=ax)
		plt.show()

	###test pvslice
	if 0:
		path = [[83.5,-2.5],[84.5,-2],[85.2,0.5],[86.5,1.5]]
		#path = SkyCoord(path, frame='galactic', unit='deg')
		#f = open('pvpath.txt', 'w')
		#for p in path: f.write('%f  %f\n' % tuple(p))
		#f.close()
		#path = 'pvpath.txt'
		from astropy.convolution import Gaussian2DKernel

		pvhdulist = cube.with_spectral_unit('km/s').pvslice(path, step=5, width=5, keep_width=True,
			spline=True)#, kernel=lambda dx,dy: ((dx**2+dy**2)<9).astype(int))
		pvhdulist.writeto('testpv.fits', overwrite=True)
		pvmap, path = pvhdulist
		
		m0 = cube.moment(order=0)
		ax = plt.subplot(projection=m0.wcs)
		ax.imshow(m0.data)
		ax.plot(path.data[...,0].T,path.data[...,1].T,'-',transform=ax.get_transform('galactic'))
		plt.show()

		ax = plt.subplot(projection=WCS(pvhdu.header,naxis=2))
		ax.imshow(pvhdu.data.mean(axis=0))
		ax.contour(pvhdu.data[-1])
		plt.show()
