;by ShaoboZhang (shbzhang@pmo.ac.cn)
;CUBERMS
;Calculate rms for a datacube
;Usage: cuberms, fitsfile, range, [window=[w1, w2, ..., wn]]
;Input:
;  cubefile: a string scalar that contains datacube filename for calculation
;  range: 2-element vector. Range over which the rms should be computed.
;Optional input keyword:
;  window: ranges avoided by baseline fitting.
;  silent: if set, then the display of the process description on the terminal will be suppressed 
;Example:
;    cuberms, 'cube.fits', [-200,200]
;    cuberms, 'cube.fits', [-100,100], window=[-50,-30,-10,10]
;History:
;Dec,23,2013,v1.0
;Oct,17,2014,v1.1
;    consider the situation that velocity range in cube is less then the given range.

function stdrange, r, rmax
;standardize the range
r = r[sort(r)]
if finite(r[0]) then r[0] = round(r[0])
if finite(r[1]) then r[1] = round(r[1])
if (r[0] gt rmax) or (r[1] lt 0) then begin
    r=-1
    return,r
endif
r = r >0 <rmax
return,r
end


pro cuberms, cubefile, range, window=window, silent=silent

if n_params() lt 1 then begin
    print, 'Syntax - CUBERMS, cubefile, [range, window=]'
    return
endif

if ~file_test(cubefile) then begin
    print, 'Error: fits datacube file does not exist!'
    return
endif
if n_elements(range) lt 2 then begin
    print, 'Error: invalid range!'
    return
endif else range=range[0:1]
if keyword_set(silent) then silent=0b else silent=1b
outname=file_basename(cubefile,'.fits',/fold_case)

if silent then print,'RMS calculation for '+cubefile+' begins:'
if silent then print,'  Loading fits file......(this may take a while!)'
fits_read, cubefile, dat, hdr
if sxpar(hdr,'BITPIX') gt 0 then nan=where(dat eq max(dat)) $
else nan=where(dat eq -1000)
if nan[0] ne -1 then dat[nan] = !values.d_nan


if silent then print,'  Checking range'
if n_elements(range) lt 1 then range=[-!values.f_infinity,!values.f_infinity] else range=range * 1000d		;km/s to m/s
nc = sxpar(hdr,'NAXIS3')
cv = sxpar(hdr,'CRVAL3')
cp = sxpar(hdr,'CRPIX3')
cd = sxpar(hdr,'CDELT3')
;v = ((findgen(nc)-cp+1)*cd+cv)
c = (range-cv)/cd+cp-1
c=stdrange(c,sxpar(hdr,'NAXIS3')-1)
if c[0] eq -1 then begin
    print,'Error: Range is out of the given datacube.'
    return
endif
flag=bytarr(nc)
flag[*]=0b
flag[c[0]:c[1]]=1b
if keyword_set(window) then begin
    window=window * 1000d	;km/s to m/s
    cw = (window-cv)/cd+cp-1
    cw = round(cw) > 0 < (nc-1)
    for i=0, n_elements(cw)/2-1 do begin
	w=cw[i*2:i*2+1]
	w=w[sort(w)]
	if w[0] ne w[1] then flag[w[0]:w[1]]=0
    endfor
endif
;stop
if total(flag) lt 3 then begin
    print,'Error: No enough range for baseline fitting.'
    return
endif


if silent then print,'  Calculating RMS'
;flag=where(flag)
dat = dat[*,*,where(flag)]
rms=sqrt(total(dat^2,3,/nan)/total(finite(dat),3))
;nan=where(~finite(rms))
;if nan[0] ne -1 then rms[nan] = -1000.

if silent then print,'  Writing FITS'
sxaddpar, hdr, 'NAXIS', 2
sxaddhist,'CUBERMS: '+systime(0),hdr
sxaddhist,'Original datacube: '+cubefile,hdr
sxaddhist,'Velocity range: '+strjoin(string(range,format='(F0)'),' '),hdr
if keyword_set(window) then sxaddhist,'Window: '+strjoin(string(window,format='(F0)'),' '),hdr
fits_write,outname+'_rms.fits',rms,hdr
if silent then print,'Done!'
end
