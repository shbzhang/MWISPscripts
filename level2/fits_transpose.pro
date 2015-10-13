pro fits_transpose,inputfile,outputfile,p
;transpose a datacube
if n_params() lt 1 then begin
	print,'Syntax - FITS_TRANSPOSE, inputfilename, outputfilename, [p]'
	return
endif
if ~file_test(inputfile) then begin
	print,'Error: fits file not exist!'
	return
endif
fits_read,inputfile,dat,hdr
;dat[where(dat eq max(dat))]=!values.f_nan
naxis=size(dat,/n_dimensions)
if n_params() lt 3 then p = reverse(indgen(naxis))
dat=transpose(dat,p)
nc={naxis:0, ctype:'', crval:0d, cdelt:0d, crpix:0d, crota:0d}
nc=replicate(nc,naxis)
for i=1,naxis do begin
	nc[i-1].naxis=sxpar(hdr,'NAXIS'+string(i,format='(I0)'))
	nc[i-1].ctype=sxpar(hdr,'CTYPE'+string(i,format='(I0)'))
	nc[i-1].crval=sxpar(hdr,'CRVAL'+string(i,format='(I0)'))
	nc[i-1].cdelt=sxpar(hdr,'CDELT'+string(i,format='(I0)'))
	nc[i-1].crpix=sxpar(hdr,'CRPIX'+string(i,format='(I0)'))
	nc[i-1].crota=sxpar(hdr,'CROTA'+string(i,format='(I0)'))
endfor
nc = nc[p]
for i=1,naxis do begin
        sxaddpar,hdr,'NAXIS'+string(i,format='(I0)'),nc[i-1].naxis
        sxaddpar,hdr,'CTYPE'+string(i,format='(I0)'),nc[i-1].ctype
        sxaddpar,hdr,'CRVAL'+string(i,format='(I0)'),nc[i-1].crval
        sxaddpar,hdr,'CDELT'+string(i,format='(I0)'),nc[i-1].cdelt
        sxaddpar,hdr,'CRPIX'+string(i,format='(I0)'),nc[i-1].crpix
        sxaddpar,hdr,'CROTA'+string(i,format='(I0)'),nc[i-1].crota
endfor
fits_write,outputfile,dat,hdr
end
