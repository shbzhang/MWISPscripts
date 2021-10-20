;Procedure for Milky Way Imaging Scroll Painting (MWISP) Project
;NAME:
;	grid2.pro
;PURPOSE:
;	Generates the template grid file
;RESPONSIBLE PERSON
;	Shaobo Zhang
;CURRENT VERSION:
;	1.1
;REVISION HISTORY:
;	V1.1
;	new procedures to generate grid for a cell / in a range / from a catalog
;	add cos(b) factor for high latitude survey
;	V1.0
;	Written by Shaobo Zhang, 2011

function _hdr
	;make a header
	mkhdr,hdr,1,[1,1]
	;;;position information
	sxaddpar,hdr,'NAXIS1',91               ;dimension
	sxaddpar,hdr,'NAXIS2',91
	sxaddpar,hdr,'CTYPE1','GLON-GLS    '   ;projection type
	sxaddpar,hdr,'CRVAL1',0                ;coordinate of reference pixel
	sxaddpar,hdr,'CDELT1',-30d/3600.       ;pixelsize
	sxaddpar,hdr,'CRPIX1',46               ;referencepixel
	sxaddpar,hdr,'CROTA1',0
	sxaddpar,hdr,'CRVAL1',0.
	sxaddpar,hdr,'CTYPE2','GLAT-GLS    '
	sxaddpar,hdr,'CRVAL2',0
	sxaddpar,hdr,'CDELT2',30d/3600.
	sxaddpar,hdr,'CRPIX2',46
	sxaddpar,hdr,'CROTA2',0
	sxaddpar,hdr,'CRVAL2',0.
	sxaddpar,hdr,'GLON',0.
	sxaddpar,hdr,'GLAT',0.
	;;;other keywords
	sxaddpar,hdr,'EQUINOX',2000
	sxaddpar,hdr,'OBJECT','DLHSURVEY'
	return, hdr
end

function _dat, hdr
	;make the grid fits data
	return, bytarr(sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2'))
end

function _id, l, b
	;make cell id from l,b
	return, string(round(l*10),format='(i04)')+string(round(b*10),format='(i+04)')
end

pro _grid, dat, hdr, l, b
	;write a grid
	id = _id(l, b)
	sxaddpar,hdr,'CRVAL1',l
	sxaddpar,hdr,'CRVAL2',b
	sxaddpar,hdr,'GLON',l
	sxaddpar,hdr,'GLAT',b
	sxaddpar,hdr,'CRPIX2',46
	;for cell beyond |b|<=5deg
	if abs(b) gt 5.1 then sxaddpar,hdr,'CDELT1',-30d/3600.*cos(b*!dtor)
	;for region -1<gb<0
	if (b gt -1) and (b lt 0) then begin
		sxaddpar,hdr,'CRVAL2',0d
		sxaddpar,hdr,'CRPIX2',round(-b*120+46)
	endif
	fits_write,id+'_grid1.fits',dat,hdr
	print,'Grid for "'+id+'" has been generated successfully!'
end

pro grid_in_cell, l, b
	l = round(l*2)/2.
	b = round(b*2)/2.
	hdr=_hdr()
	dat=_dat(hdr)
	_grid, dat, hdr, l ,b
end

pro grid_in_range, l0, l1, b0, b1
	l0 = ceil(l0*2)/2.
	l1 = floor(l1*2)/2.
	b0 = ceil(b0*2)/2.
	b1 = floor(b1*2)/2.
	hdr=_hdr()
	dat=_dat(hdr)
	for l = l0, l1, 0.5 do begin
		for b = b0, b1, 0.5 do begin
			_grid, dat, hdr, l, b
		endfor
	endfor
end

pro grid_in_cat, cat
	if ~arg_present(cat) then cat = 'mymap_v4.cat'
	if ~file_test(cat) then begin
		print,'cat not exist!'
		return
	endif
	readcol,cat,num,id,gl,gb,format = 'i,a,d,d',skipline = 1
	hdr=_hdr()
	dat=_dat(hdr)
	for i = 0, 4 do _grid, dat, hdr, gl[i], gb[i]
	print, 'Done'
end
