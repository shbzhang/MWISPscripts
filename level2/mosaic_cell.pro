;by ShaoboZhang
;Mosaic fits file of DLH survey
;Usage: for USB:mosaic, gl, gb    or    for LSB:mosaic, gl, gb, /lower
;Input:
;  gl,gb: two float scalar. Galactic coordinate for the mosaic center.
;History:
;Nov,09,2011,v1.0

function getcellname, gl, gb
;get cell name from its gl and gb
gl = round(gl*2)/2d
gb = round(gb*2)/2d
gl = gl mod 360
return,string(fix(gl*10),format='(I04)')+string(fix(gb*10),format='(I+04)')
end



pro clipv, dat, hdr, v1, v2
;clip a velocity range (v in km/s) in data according to fits header
nc = sxpar(hdr,'NAXIS3')
cv = sxpar(hdr,'CRVAL3')
cp = sxpar(hdr,'CRPIX3')
cd = sxpar(hdr,'CDELT3')
n = ([v1,v2]*1000d -cv)/cd+cp-1
n = n[sort(n)]
n = [floor(n[0]),ceil(n[1])]
if (n[0] gt nc) or (n[1] lt 0) then begin
        print,'Velocity out of range.'
	hdr='OUT'
        return
endif
n = n >0 <nc
print,'Clip channel between:',n
dat = dat[*,*,n[0]:n[1]]
sxaddpar,hdr,'NAXIS3',n[1]-n[0]+1
sxaddpar,hdr,'CRPIX3',cp-n[0]
end



pro mosaic_cell, gl, gb, lower=lower
;mosaic a 30'*30' cell
if keyword_set(lower) then sb = 'L' else sb = 'U'
gl = round(gl*2)/2d
gb = round(gb*2)/2d
print,'Mosaic '+sb+'SB cell data at '+string(gl,format='(F05.1)')+', '+string(gb,format='(F+5.1)')

name=getcellname(gl,gb)
if ~file_test(name+sb+'.fits') then begin
  print, '  Center fits file "'+name+sb+'.fits" does not exist!'
  return
endif
fits_read,name+sb+'.fits',dat,hdr
nx=sxpar(hdr,'NAXIS1')
ny=sxpar(hdr,'NAXIS2')
cx=(nx+1)/2-1
cy=(ny+1)/2-1
dat[where(dat eq max(dat))] = !values.d_nan
dat=dat[(cx-30):(cx+30),(cy-30):(cy+30),*]
fits_read,name+sb+'_rms.fits',rms
rms[where(rms eq max(rms))] = !values.d_infinity
rms=rms[(cx-30):(cx+30),(cy-30):(cy+30)]
wei=1d/rms^2
for i=0l,16383l do dat[*,*,i]=dat[*,*,i]*wei

sxaddpar,hdr,'NAXIS1',61
sxaddpar,hdr,'NAXIS2',61
sxaddpar,hdr,'CRPIX1',31d
sxaddpar,hdr,'CRPIX2',31d
sxaddpar,hdr,'CRVAL1',gl
sxaddpar,hdr,'CRVAL2',gb
sxaddpar,hdr,'CTYPE1','GLON-CAR    '
sxaddpar,hdr,'CTYPE2','GLAT-CAR    '

nearx=[1,1,1,0,0,-1,-1,-1]
neary=[1,0,-1,1,-1,1,0,-1]
for r=0,7 do begin
  mname=getcellname(gl+nearx[r]*0.5,gb+neary[r]*0.5)
  if ~file_test(mname+sb+'.fits') then begin
    print, '['+string(r+1,format='(I1)')+'/8]Neighbor fits file "'+mname+sb+'.fits" does not exist!'
    continue
  endif
  print, '['+string(r+1,format='(I1)')+'/8]Sum neighbor fits file "'+mname+sb+'.fits'
  fits_read,mname+sb+'.fits',mdat
  mdat[where(mdat eq max(mdat))] = 0
  mdat=shift(mdat, -nearx[r]*60, neary[r]*60, 0)
  mdat=mdat[(cx-30):(cx+30),(cy-30):(cy+30),*]
  fits_read,name+sb+'_rms.fits',mrms
  mrms[where(mrms eq max(mrms))] = !values.d_infinity
  mwei=dblarr(300,300)
  mwei[0:(nx-1),0:(ny-1)] = 1d/mrms^2
  mwei = shift(mwei, -nearx[r]*60, neary[r]*60)
  mwei = mwei[0:(nx-1),0:(ny-1)]
  mwei=mwei[(cx-30):(cx+30),(cy-30):(cy+30)]

  for i=0l,16383l do mdat[*,*,i]=mdat[*,*,i]*mwei

  dat += mdat
  wei += mwei
endfor
tvscl,congrid(wei,300,300)
for i=0l,16383l do dat[*,*,i]=dat[*,*,i]/wei
fits_write,name+sb+'_cell.fits',temporary(dat),hdr
fits_write,name+sb+'_weight.fits',wei,hdr
print,'Create mosaic cell fits '+name+sb+'_cell.fits'
end

pro mosaic, lmin,lmax, bmin,bmax, vmin,vmax, skip_cell=skip_cell

;l1 = l1 mod 360 - 360*(l1 gt 300)
;l2 = l2 mod 360 - 360*(l2 gt 300)
;nx = long(l2*120) long(l1*120)+1

lmin = (ceil(lmin*2)/2. mod 360) - 360*(lmin gt 300)
lmax = (floor(lmax*2)/2. mod 360) - 360*(lmax gt 300)
bmin = ceil(bmin*2)/2.
bmax = floor(bmax*2)/2.
if lmin lt -10.5 or lmax gt 250.5 $
or abs(bmin) gt 5.5 or abs(bmax) gt 5.5 then begin
  print,'Galactic latitude is out of DLH survey range!'
  print,'l = 350 ~ 0 ~ 250 deg, b = -5 ~ 5 deg'
  return
endif
print,'Mosaic begins:'
print,'from '+string(lmin mod 360)+' to '+string(lmax mod 360)
print,'from '+string(bmin)+' to '+string(bmax)

;mosaic cell
if ~keyword_set(skip_cell) then begin
print,'Mosaic each cell in this region:'
for l=lmax,lmin,-0.5 do begin
for b=double(bmin),bmax,0.5 do begin
mosaiccell,l,b
mosaiccell,l,b,/lower
endfor
endfor
endif

;mosaic
count=0l
for l=double(lmax),lmin,-0.5 do begin
for b=double(bmin),bmax,0.5 do begin
name = getcellname(l,b)
sb = 'U'
if file_test(name+sb+'_mosaic.fits') then begin
   count += 1
   fits_read,name+sb+'_mosaic.fits',dat,hdr
   clipv, dat, hdr, vmin, vmax
   if hdr[0] eq 'OUT' then return
   if count eq 1 then begin
      nx=60*((lmax-lmin)/0.5+1)+1
      ny=60*((bmax-bmin)/0.5+1)+1
      nv=sxpar(hdr,'NAXIS3')
      mdat=dblarr(nx,ny,nv)
      mhdr=hdr
      sxaddpar,mhdr,'NAXIS1',nx
      sxaddpar,mhdr,'NAXIS2',ny
      sxaddpar,mhdr,'NAXIS3',nv
      sxaddpar,mhdr,'CRPIX1',(lmax-l)*120+31
      sxaddpar,mhdr,'CRPIX2',(b-bmin)*120+31
   endif
   left=(lmax-l)*120
   down=(b-bmin)*120
help,dat,mdat
   mdat[left:(left+60),down:(down+60),*]=temporary(dat)
endif
endfor
endfor
if count eq 0 then print,'No fits file in this region.' else $
fits_write,'MOSAIC.fits',mdat,mhdr
end
