pro cubeclip, cubefile, rmsfile=rmsfile, coveragefile=coveragefile, outputfile=outputfile
;clip noise in a datecube, only keep 3 continuous channels over 3sigma
if n_params() lt 1 then begin
    print, 'Syntax - CUBECLIPNOISE, cubefile, [rmsfile=, coveragefile=, outputfile=]'
    return
endif
if ~file_test(cubefile) then begin
    print, 'Error - fits datacube file does not exist!'
    return
endif
if ~keyword_set(outputfile) then outputfile=file_basename(cubefile,'.fits',/fold_case)+'_clip.fits'

fits_read, cubefile, dat, hdr
if keyword_set(rmsfile) then begin
	fits_read, rmsfile, rms
	rms = sqrt(total(dat^2,3,/nan)/total(finite(dat),3))
endif else fits_read, rmsfile, rms

nv = sxpar(hdr,'NAXIS3')
cleast = 3
msk = bytarr(size(dat,/dimension)+[0,0,cleast-1])
for i=0,nv-1 do msk[*,*,i] = dat[*,*,i] gt rms*3
msk = msk and shift(msk,0,0,1) and shift(msk,0,0,2)
msk = msk or shift(msk,0,0,-1) or shift(msk,0,0,-2)

dat *= msk[*,*,0:nv-1]

if keyword_set(coveragefile) then begin
	fits_read, coveragefile, cov
	cov = float(cov)
	cov[where(~cov)]=!values.f_nan
	for i=0,nv-1 do dat[*,*,i] *= cov
endif
fits_write,outputfile,dat,hdr
end
