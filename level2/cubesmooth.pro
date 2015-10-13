pro cubesmooth, oldfits, oldresolution, newresolution, newfits
;smooth a datacube fits to a new resolution
if n_params() lt 3 then begin
	print, 'Syntax - CUBESMOOTH, oldfits, old_resolution, new_resolution, [ newfits ]'
	return
endif
if n_elements(newfits) lt 1 then newfits = 'cubesmooth.fits'
if total(oldresolution gt newresolution) gt 0 then begin
	print, 'Error - new resolution is smaller than the old one'
	return
endif

fits_read,oldfits,olddat,oldhdr
pixelsize = abs([sxpar(oldhdr,'CDELT1'), sxpar(oldhdr,'CDELT2')])
kwidth = sqrt(newresolution^2-oldresolution^2) / pixelsize
npixel = fix(kwidth*2.5)
npixel += ~(npixel mod 2)
kernel = psf_gaussian(npixel = npixel, fwhm = kwidth, /normalize)

for i=0,sxpar(oldhdr,'NAXIS3')-1 do olddat[*,*,i] = convol(olddat[*,*,i],kernel)

sxaddhist,'Smooth spatial resolution from '+string(oldresolution, format='(F0)')+' to '+string(newresolution, format='(F0)'), oldhdr
fits_write, newfits, olddat, oldhdr
end

