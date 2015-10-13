function averspec, cubefile, rmsfile=rmsfile
;average spectra of a cube
if n_params() lt 1 then begin
	print, 'Syntax - REPROJECT, cubefile, [rmsfile=, range=]'
	return, 0
endif
if ~file_test(cubefile) then begin
	print, 'Error - datacube file does not exist!'
	return, 0
endif
if keyword_set(rmsfile) then begin
	if ~file_test(rmsfile) then begin
	    print, 'Error - RMS file does not exist!'
	    return, 0
	endif
	fits_read,rmsfile,rms
endif else rms = 1
wei = 1/rms^2

fits_read, cubefile, dat, hdr
nc=sxpar(hdr,'NAXIS3')
spec=dblarr(nc)

if keyword_set(rmsfile) then begin
	for i=0,nc-1 do spec[i] = total(dat[*,*,i]*wei, /nan)/total(finite(dat[*,*,i])*wei, /nan)
endif else begin
	for i=0,nc-1 do spec[i] = total(dat[*,*,i], /nan)/total(finite(dat[*,*,i]))
endelse

c=findgen(nc)
c2v,hdr,c,v

return,[[v],[spec]]
end
