;by ShaoboZhang (shbzhang@pmo.ac.cn)
;MWISPMERGE
;  merge 2-dimensional file into the galactic plane
;Optional input:
;  imagefile: a string array that contains the name of 2-dimensional files (e.g. generate by CUBEMOMENT, CUBERMS)
;  coveragefile: a string array that contains the name of coverage files of the same region (e.g. generate by MOSAIC)
;Optional input keywords:
;  mergefile: a string scalar that contains the name of galactic plane file for both input and output.
;Example:
;  put a moment 0 file into a new file of galactic plane:
;    mwispmerge, 'mosaic_U_m0.fits', 'mosaic_U_coverage.fits'
;  put several moment 0 files into an existing file of galactic plane:
;    mwispmerge, ['a_m0.fits','b_m0.fits'], ['a_coverage.fits','b_coverage.fits'], mergefile='mwisp_m0.fits'
;  generate a file of empty galactic plane:
;    mwispmerge, mergefile='mwisp.fits'
;History:
;Sep,25,2015,v1.0
pro mwispmerge, imagefile, coveragefile, mergefile=mergefile
	if n_params() lt 1 and ~keyword_set(mergefile) then begin
	    print, 'Syntax - MWISPMERGE, imagefile, [coveragefile, [mergefile=]]'
	    print, 'Syntax - MWISPMERGE, mergefile='
	    return
	endif
	if ~keyword_set(mergefile) then mergefile='mwisp_image.fits'
	if file_test(mergefile) then fits_read, mergefile, bgdat, bghdr else emptygplane, bgdat, bghdr
	if n_params() lt 1 then begin
		fits_write, mergefile, bgdat, bghdr
		return
	endif
	bgdim = size(bgdat,/dimension)

	for i=0,n_elements(imagefile)-1 do begin
		fits_read, imagefile[i], mdat, mhdr
		offsetx = sxpar(bghdr,'CRPIX1')-sxpar(mhdr,'CRPIX1')
		offsety = sxpar(bghdr,'CRPIX2')-sxpar(mhdr,'CRPIX2')

		if n_params() ge 2 and i lt n_elements(coveragefile) then fits_read, coveragefile[i], cov else begin
			cov = byte(mdat)
			cov[*] = 1b
		endelse
		idx = where(cov, count)
		if count eq 0 then begin
			print, 'Warning - region does not cover any pixel: '+imagefile[i]
			continue
		endif
		ii = mdat[idx]
		xy = array_indices(cov, idx)
		x = xy[0,*] + offsetx	;xxx
		y = xy[1,*] + offsety	;xxx
		in = where(y ge 0 and y le bgdim[1]-1, count)
		if count eq 0 then begin
			print, 'Warning - region outside the range of current survey: '+imagefile[i]
			continue
		endif
		x = (x[in]+bgdim[0]) mod bgdim[0]
		y = y[in]
		ii = ii[in]
		bgdat[x,y] = ii
	endfor
	fits_write, mergefile, bgdat, bghdr
end


pro emptygplane, bgdat, bghdr, header_only=header_only
;generate a empty galactic plane
	xdim = 360l*3600/30	;180~0~-180
	ydim = 12l*3600/30+1	;-6~6
	mkhdr, bghdr, 4, [xdim, ydim]
	sxaddpar,bghdr,'BUNIT','K (T_MB) * km/s'
	sxaddpar,bghdr,'CTYPE1','GLON-CAR'
	sxaddpar,bghdr,'CRVAL1',0d
	sxaddpar,bghdr,'CRPIX1',180l*3600/30+1
	sxaddpar,bghdr,'CDELT1',-30d/3600
	sxaddpar,bghdr,'CROTA1',0d
	sxaddpar,bghdr,'CTYPE2','GLAT-CAR'
	sxaddpar,bghdr,'CRVAL2',0d
	sxaddpar,bghdr,'CRPIX2',6l*3600/30+1
	sxaddpar,bghdr,'CDELT2',30d/3600
	sxaddpar,bghdr,'CROTA2',0d
	sxaddpar,bghdr,'EQUINOX',0d
	if keyword_set(header_only) then return
	bgdat = fltarr(xdim,ydim)
	bgdat[*] = !values.f_nan
end

pro damegplane
;resample Dame2001 data to MWISP sample
	fits_read,'Wco_DHT2001.fits',dat,hdr
	emptygplane, bgdat, bghdr, /header_only
	bgdim = [sxpar(bghdr,'NAXIS1'),sxpar(bghdr,'NAXIS2')]
	adxy,bghdr,0,0,x0,y0
	xyad,bghdr,findgen(bgdim[0]),replicate(y0,bgdim[0]),glx,gbx
	xyad,bghdr,replicate(x0,bgdim(1)),findgen(bgdim[1]),gly,gby
	adxy,hdr,glx,gbx,x,t
	adxy,hdr,gly,gby,t,y
	xx = round(x)#replicate(1,bgdim[1])
	yy = replicate(1,bgdim[0])#round(y)
	bgdat = dat[xx, yy]
	fits_write,'Wco_DHT2001_mwisp.fits',bgdat,bghdr
end
