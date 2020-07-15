from astropy.io import fits
import astropy.wcs as wcs
from reproject import reproject_interp
from numpy import flip

def reproject_fits(oldfits, reffits, newfits='reproject.fits', order = 'bilinear'):
	'''
	reproject a fits to a new ref wcs

	Parameters
	----------
	oldfits : str
		the file name of the fits needed to be reproject.
	reffits : str
		the file name of the reference fits to reproject to.
	newfits : str, optional
		the file name for output.

	Examples
	--------
	>>> reproject_fits('old.fits','ref.fits','new.fits')
	'''
	#modify header as new fits header
	oldhdr = fits.getheader(oldfits)
	newhdr = oldhdr.copy()
	refhdr = fits.getheader(reffits)
	refhdr2 = wcs.WCS(refhdr, naxis=2).to_header()	#only keep celestial axis
	refhdr2.remove('WCSAXES')
	for kw in refhdr2:
		newhdr[kw] = refhdr2[kw]
	newhdr['NAXIS1'] = refhdr['NAXIS1']
	newhdr['NAXIS2'] = refhdr['NAXIS2']

	#get wcs
	olddat = fits.getdata(oldfits)
	oldwcs = wcs.WCS(oldhdr, naxis=olddat.ndim)
	newwcs = wcs.WCS(newhdr, naxis=olddat.ndim)

	#reproject
	newdat, coverage = reproject_interp((olddat, oldwcs), newwcs, shape_out = flip(newwcs._naxis), 
		order = order, independent_celestial_slices=True)

	#output
	newhdu = fits.PrimaryHDU(newdat, header=newhdr)
	newhdu.writeto(newfits, overwrite=True)

if __name__ == '__main__':
	#test example
	old = '/Users/shaobo/Work/mwips/L935/mosaic_L_m0.fits'
	ref = '/Users/shaobo/Work/script/mylib/python/demo/data/img_icrs.fits'
	reproject_fits(old,ref,'new.fits')