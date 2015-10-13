;by ShaoboZhang (shbzhang@pmo.ac.cn)
;CUBEMOMENT
;Calculate moment for a datacube
;Usage: cubemoment, fitsfile, range [, direction='v', threshold=-10, /zeroth_only, outname='Region' $
;                [, /goodlooking, rmsfile='rms.fits', coveragefile='coverage.fits']]
;Input:
;  cubefile: a string scalar that contains datacube filename for calculation
;  range: 2-element vector. Range over which the moments should be computed.
;        For 'L' or 'B' direction, the unit should be arcdeg; for 'V', the unit should be km/s
;Optional input keyword:
;  direction: default calculation is done along the "Velocity" direction. 
;        Other directions as "L" or "B" could be provided to calculate along the gl or gb
;  threshold: the minimum intensity under which a pixel should no longer be considered.
;  zeroth_only: set this keyword if you only need the 0th moments
;Optional goodlooking mode:
;  In the "goodlooking" mode, only keep those pixels with at least 3 channels above 3 sigma,
;  and only use pixels above 3 sigma for integration.
;  This mode will make the 1st and 2nd moments more reliable.
;Input keyword:
;  goodlooking: set to 1 to active the goodlooking mode
;  rmsfile: a string scalar that contains the filename of sigma map for the datacube, 
;  coveragefile: a string scalar that contains the filename of coverage map for the datacube 
;         default is cubefile(no suffix)+'_coverage.fits'
;Example:
;  for general use:
;    cubemoment, 'cube.fits', [-10,10]
;    cubemoment, 'cube.fits', [-10,10], threshold=1, /zeroth_only
;  integrate along b direction to create an l-v map:
;    cubemoment, 'cube.fits', [-2,-1], direction='b'
;  goodlooking mode:
;    cubemoment, 'cube.fits', [-10,10], /goodlooking, rmsfile='rms.fits', coveragefile='coverage.fits'
;  for mosaiced PMO survey data produce by mosaic.pro
;    cubemoment, 'mosaic_U.fits', [-10,10], /goodlooking
;History:
;May,08,2012,v1.0
;Dec,25,2013,v1.1
;  goodlooking mode: only keep those pixels with at least 3 "continuous" channels above 3 sigma
;Jun,24,2014,v1.2
;  fix minor error when apply to all FITS file

function stdrange, r, rmax
;standardize the range
r = r[sort(r)]
r = [round(r[0]),round(r[1])]
if (r[0] gt rmax) or (r[1] lt 0) then begin
    r=-1
    return,r
endif
r = r >0 <rmax
return,r
end


pro cubemoment, cubefile, range, direction=direction, threshold=threshold, zeroth_only=zeroth_only, $
    coveragefile=coveragefile, goodlooking=goodlooking, rmsfile=rmsfile, outname=outname

if n_params() lt 2 then begin
    print, 'Syntax - CUBEMOMENT, cubefile, range, [direction=, threshold=, /zeroth_only,'
    print, '             coveragefile=, /goodlooking, rmsfile=, outname=]'
    return
endif

if ~file_test(cubefile) then begin
    print, 'Error: fits datacube file does not exist!'
    return
endif
if ~keyword_set(outname) then outname=file_basename(cubefile,'.fits',/fold_case)
if n_elements(range) lt 2 then begin
    print, 'Error: invalid range!'
    return
endif else range=range[0:1]

if ~keyword_set(direction) then direction='V'
direction = where(strcmp(direction,['L','B','V'],1,/fold_case))+1
if direction eq 0 then begin
    print, "Error: Invalid direction!"
    return
endif
;  if need threshold and only zeroth moment
if n_elements(threshold) eq 0 then threshold=-!values.f_infinity
if ~keyword_set(zeroth_only) then zeroth_only=0
if ~keyword_set(coveragefile) then coveragefile=outname+'_coverage.fits'
if ~file_test(coveragefile) then begin
        print, 'Warning: no valid coverage file! No coverage will be used.'
        cov=0
endif else begin
        print, 'Note: use coveragefile '+coveragefile
        cov=1
endelse
if keyword_set(goodlooking) then begin
    if ~keyword_set(rmsfile) then begin
        rmsfile=outname+'_rms.fits'
        if ~file_test(rmsfile) then begin
            print, "Error: rmsfile is not provided!"
            return
        endif else print, "Note: Use rmsfile "+rmsfile
    endif else begin
        if ~file_test(rmsfile) then begin
            print, "Error: rmsfile does not exist!"
            return
        endif else print, "Note: Use rmsfile "+rmsfile
    endelse
endif else goodlooking=0

print,'Moment begins:'
print,'  Loading fits file......(this may take a while!)'
fits_read, cubefile, dat, hdr
;dat[where(dat eq max(dat))] = !values.f_nan
invalid = where(dat lt threshold, /l64)
if invalid[0] ne -1 then dat[invalid] = 0
if cov then begin
    fits_read,coveragefile,cov
    if total(size(dat[*,*,0],/dim) eq size(cov,/dim)) ne 2 then $
        print, "Warning: dimension of coverage file does not match! No coverage will be used." $
    else begin
        cov = float(cov ne 0)
        nan = where(cov eq 0)
        if nan[0] ne -1 then cov[nan]=!values.f_nan
        for i=0,n_elements(dat[0,0,*])-1 do dat[*,*,i] *= cov
    endelse
endif
sxaddpar, hdr, 'NAXIS', 2
sxaddhist,'CUBEMOMENT: '+systime(0),hdr
sxaddhist,'Original datacube: '+cubefile,hdr

case direction of
    1:begin
    print,'  Creating b-v map.'
    adxy,hdr,range,replicate(sxpar(hdr,'CRVAL2'),2),x,y
    x=stdrange(x,sxpar(hdr,'NAXIS1')-1)
    if x[0] eq -1 then begin
        print,'Error: Range is out of the given datacube.'
        return
    endif
    xyad,hdr,[x[0]-0.5,x[1]+0.5],y,realrange,d
    dat = temporary(dat[x[0]:x[1],*,*])
    flag = total(finite(dat),1) ne 0
    out = where(flag eq 0)
    map = total(dat,1,/nan) * abs(sxpar(hdr,'CDELT1'))
    if out[0] ne -1 then map[out] = !values.f_nan
    sxaddpar,hdr,'BUNIT','K * arcdeg'
    sxaddpar,hdr,'NAXIS1',sxpar(hdr,'NAXIS2')
    sxaddpar,hdr,'CTYPE1',sxpar(hdr,'CTYPE2')
    sxaddpar,hdr,'CRVAL1',sxpar(hdr,'CRVAL2')
    sxaddpar,hdr,'CRPIX1',sxpar(hdr,'CRPIX2')
    sxaddpar,hdr,'CDELT1',sxpar(hdr,'CDELT2')
    sxaddpar,hdr,'NAXIS2',sxpar(hdr,'NAXIS3')
    sxaddpar,hdr,'CTYPE2',sxpar(hdr,'CTYPE3')
    sxaddpar,hdr,'CRVAL2',sxpar(hdr,'CRVAL3')
    sxaddpar,hdr,'CRPIX2',sxpar(hdr,'CRPIX3')
    sxaddpar,hdr,'CDELT2',sxpar(hdr,'CDELT3')
    sxdelpar,hdr,['NAXIS','CTYPE','CRVAL','CDELT','CRPIX','CROTA']+'3'
    sxdelpar,hdr,['NAXIS','CTYPE','CRVAL','CDELT','CRPIX','CROTA']+'4'
    prompt='Integrate L between '+strjoin(string(realrange,format='(F0)'),' and ')
    print,'  '+prompt
    sxaddhist,prompt,hdr
    fits_write,outname+'_bvmap.fits',map,hdr
    end

    2:begin
    print,'  Creating l-v map.'
    adxy,hdr,replicate(sxpar(hdr,'CRVAL1'),2),range,x,y
    y=stdrange(y,sxpar(hdr,'NAXIS2')-1)
    if y[0] eq -1 then begin
        print,'Error: Range is out of the given datacube.'
        return
    endif    
    xyad,hdr,x,[y[0]-0.5,y[1]+0.5],a,realrange
    dat = temporary(dat[*,y[0]:y[1],*])
    flag = total(finite(dat),2) ne 0
    out = where(flag eq 0)
    map = total(dat,2,/nan) * abs(sxpar(hdr,'CDELT2'))
    if out[0] ne -1 then map[out] = !values.f_nan 
    sxaddpar,hdr,'BUNIT','K * arcdeg'
    sxaddpar,hdr,'NAXIS2',sxpar(hdr,'NAXIS3')
    sxaddpar,hdr,'CTYPE2',sxpar(hdr,'CTYPE3')
    sxaddpar,hdr,'CRVAL2',sxpar(hdr,'CRVAL3')
    sxaddpar,hdr,'CRPIX2',sxpar(hdr,'CRPIX3')
    sxaddpar,hdr,'CDELT2',sxpar(hdr,'CDELT3')
    sxdelpar,hdr,['NAXIS','CTYPE','CRVAL','CDELT','CRPIX','CROTA']+'3'
    sxdelpar,hdr,['NAXIS','CTYPE','CRVAL','CDELT','CRPIX','CROTA']+'4'
    prompt='Integrate B between '+strjoin(string(realrange,format='(F0)'),' and ')
    print,'  '+prompt
    sxaddhist,prompt,hdr
    fits_write,outname+'_lvmap.fits',map,hdr
    end

    3:begin
    print,'  Creating 0th moment map.'
    range=range * 1000d		;km/s to m/s
    nc = sxpar(hdr,'NAXIS3')
    cv = sxpar(hdr,'CRVAL3')
    cp = sxpar(hdr,'CRPIX3')
    cd = sxpar(hdr,'CDELT3')
    v = ((findgen(nc)-cp+1)*cd+cv)
    c = (range-cv)/cd+cp-1
    c=stdrange(c,sxpar(hdr,'NAXIS3')-1)
    if c[0] eq -1 then begin
        print,'Error: Range is out of the given datacube.'
        return
    endif
    realrange = interpolate(v,[c[0]-0.5,c[1]+0.5]) / 1000d
    dat = temporary(dat[*,*,c[0]:c[1]])
    v = v[c[0]:c[1]]
    ;find positions with all channels equal nan
    flag = total(finite(dat),3) ne 0
    ;goodlooking
    if goodlooking then begin
        fits_read,rmsfile,rms
        if total(size(dat[*,*,0],/dim) eq size(rms,/dim)) ne 2 then begin
            print, "Error: dimension rms file does not match!"
            return
        endif	;check for consistent

        ;only keep those pixels with at least 3 channels above 3 sigma
;        msk = intarr(size(rms,/dim))
;        for i=0, n_elements(v)-1 do msk += (dat(*,*,i) gt rms*3)	;s/n gt 3
;        msk = msk ge 3	;at least 3 channels
;        for i=0, n_elements(v)-1 do dat[*,*,i] *= msk
;
;        ;only use pixel above 3 sigma for integration
;        for i=0, n_elements(v)-1 do dat[*,*,i] *= (dat[*,*,i] ge rms*3)

        ;only keep those pixels with at least 3 #continuous# channels above 3 sigma
        rmscube=make_array(size=size(dat))
        for i=0, n_elements(v)-1 do rmscube[*,*,i]=rms
        mskcube = dat gt rmscube*3
;	nchannel=total(mskcube,3)
        mskcube = mskcube and shift(mskcube,0,0,1) and shift(mskcube,0,0,2)
        msk = total(mskcube[*,*,2:*],3) ne 0
;	fits_write,'detection.fits',msk,hdr
;	fits_write,'nchannel.fits',nchannel*msk,hdr
        for i=0, n_elements(v)-1 do dat[*,*,i] *= msk

	;only use pixel above 3 sigma for integration
        dat *= (dat gt rmscube*3)
    endif
    out = where(flag eq 0)
    sumi = total(dat,3,/nan)
    map = sumi * abs(sxpar(hdr,'CDELT3')) / 1000d
    if out[0] ne -1 then map[out] = !values.f_nan
    sxaddpar,hdr,'BUNIT','K * km / s'
    sxdelpar,hdr,['NAXIS','CTYPE','CRVAL','CDELT','CRPIX','CROTA']+'3'
    sxdelpar,hdr,['NAXIS','CTYPE','CRVAL','CDELT','CRPIX','CROTA']+'4'
    prompt='Channel between '+strjoin(string(c+1,format='(I0)'),' and ')
    print,'    '+prompt
    sxaddhist,prompt,hdr
    prompt='V between '+strjoin(string(realrange,format='(F0)'),' and ')
    print,'    '+prompt
    sxaddhist,prompt,hdr
    fits_write,outname+'_m0.fits',map,hdr
    if ~zeroth_only then begin
        print,'  Creating 1st moment map.'
        temp=dat
        for i=0, n_elements(v)-1 do temp[*,*,i] *= v[i]
        nan = where(sumi eq 0)
        if nan[0] ne -1 then sumi[nan] = !values.f_nan
        map = total(temporary(temp),3,/nan) / sumi / 1000d
        if out[0] ne -1 then map[out] = !values.f_nan
        sxaddpar,hdr,'BUNIT','km / s'
        fits_write,outname+'_m1.fits',map,hdr
        print,'  Creating 2nd moment map.'
        for i=0, n_elements(v)-1 do dat[*,*,i] *= (v[i]-map*1000d)^2
        map = sqrt(total(dat,3,/nan) / sumi * 8*alog(2)) / 1000d
        if out[0] ne -1 then map[out] = !values.f_nan
        sxaddpar,hdr,'BUNIT','km / s'
        fits_write,outname+'_m2.fits',map,hdr
    endif
    end
endcase
print,'Done!'

end

