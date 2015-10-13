pro reproject, oldfits, reffits, newfits
;reproject a datacube fits according to ref fits
;you may need to rebin the oldfits as same as reffits first!!!
if n_params() lt 2 then begin
	print, 'Syntax - REPROJECT, oldfits, reffits, [ newfits ]'
	return
endif
if n_elements(newfits) lt 1 then newfits = 'reproject.fits'

fits_read,oldfits,olddat,oldhdr
fits_read,reffits,refdat,refhdr,/header_only
oldprj = strmid(sxpar(oldhdr,'CTYPE1'),0,2) eq 'RA'
refprj = strmid(sxpar(refhdr,'CTYPE1'),0,2) eq 'RA'

;resample?
;nx = sxpar(oldhdr, 'NAXIS1')
;ny = sxpar(oldhdr, 'NAXIS2')
nv = sxpar(oldhdr, 'NAXIS3')
;resdat = make_array(nx, ny, nv, type = size(olddat, /type))


;reproject
nx = sxpar(refhdr, 'NAXIS1')
ny = sxpar(refhdr, 'NAXIS2')
newdat = make_array(nx, ny, nv ,type = size(olddat, /type))
x = dindgen(nx)#replicate(1,ny)
y = replicate(1,nx)#dindgen(ny)
xyad, refhdr, x, y, a, d
if oldprj and ~refprj then begin
	glactc,at,dt,2000,a,d,2,/degree
	a=at & d=dt
end
if ~oldprj and refprj then begin
	glactc,a,d,2000,at,dt,1,/degree
	a=at & d=dt
endif
adxy, oldhdr, a, d, x, y
for i=0, nv-1 do newdat[*,*,i] = interpolate(olddat[*,*,i],x,y,missing=!values.f_nan)
sxaddpar,refhdr,'NAXIS',3
sxaddpar,refhdr,'NAXIS3',sxpar(oldhdr,'NAXIS3')
sxaddpar,refhdr,'CTYPE3',sxpar(oldhdr,'CTYPE3')
sxaddpar,refhdr,'CRVAL3',sxpar(oldhdr,'CRVAL3')
sxaddpar,refhdr,'CDELT3',sxpar(oldhdr,'CDELT3')
sxaddpar,refhdr,'CRPIX3',sxpar(oldhdr,'CRPIX3')
sxaddpar,refhdr,'CROTA3',sxpar(oldhdr,'CROTA3')
fits_write, newfits, newdat, refhdr
end

