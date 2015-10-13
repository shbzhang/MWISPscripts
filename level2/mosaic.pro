;by ShaoboZhang
;Mosaic fits file of DLH survey
;Usage: mosaic, gl1, gl2, gb1, gb2, velo1, velo2 [, sb='L', path='./', /display]
;Input:
;  gl1,gl2,gb1,gb2,velo1,velo2: float scalar. Galactic coordinate range and velocity range.
;Input keyword:
;  sb: sideband. for 12CO, sb='U'; for 13CO, sb='L'; for C18O, sb='L2'
;  path: array of the directories containing datacube and rms files, ended with '/'
;  dispaly: whether to display the process on a window. defalut is 1.
;Note:
;  all the input fits file (data and rms) must be the same as the grid define by template.
;History:
;Nov,09,2011,v1.0
;Nov,28,2011,v1.1
;  add the fitspath and rmspath keywords.
;Jan,16,2012,v1.2
;  fix a error when mosaic region with gl > 180
;Apr,12,2012,v1.3
;  change the output file name
;  output a mask fits file for the mosaic
;May,09,2012,v1.4
;  complete the header of the result fits file
;  the mask file is now called 'coverage' file to fit its function
;  add the display keywords for users who do not support window
;Dec,25,2013,v1.5
;  accept rms file created by procedure cuberms
;Oct,13,2015,v1.6
;  remove keyword "fitspath" and "rmspath"
;  new keyword "path" now accepts string array, procedure will 
;    automatically search datacube and rms file in these directories

function getcellname, gl, gb
;get cell name from its gl and gb
gl = round(gl*2)/2d
gb = round(gb*2)/2d
gl = gl mod 360
return,string(fix(gl*10),format='(I04)')+string(fix(gb*10),format='(I+04)')
end

pro lowhigh, a, b
t=max([a,b])
a=min([a,b])
b=t
end

pro clipv, dat, hdr, v1, v2
;clip a velocity range (v in km/s) in data according to fits header
nc = sxpar(hdr,'NAXIS3')
cv = sxpar(hdr,'CRVAL3')
cp = sxpar(hdr,'CRPIX3')
cd = sxpar(hdr,'CDELT3')
if cd lt 0 then begin
	dat = reverse(dat,3)
	cd = abs(cd)
	cp = nc+1-cp
endif
n = ([v1,v2]*1000d -cv)/cd+cp-1
n = n[sort(n)]
n = [floor(n[0]),ceil(n[1])]
if (n[0] gt nc) or (n[1] lt 0) then begin
        print,'Velocity out of range.'
	hdr='OUT'
        return
endif
n = n >0 <nc
;print,'Clip channel between:',n
dat = dat[*,*,n[0]:n[1]]
sxaddpar,hdr,'NAXIS3',n[1]-n[0]+1
sxaddpar,hdr,'CRPIX3',cp-n[0]
sxaddpar,hdr,'CDELT3',cd
end



pro mosaic, l1, l2, b1, b2, v1, v2, sb=sb, path=path, display=display

if n_params() lt 6 then begin
    print, 'Syntax - MOSAIC, l1, l2, b1, b2, v1, v2, [sb=, path=, /display ]'
    return
endif
if ~keyword_set(sb) then sb='U'
if ~keyword_set(path) then path=''
if n_elements(display) eq 0 then display=1

l1 = l1 mod 360 - 360*(l1 gt 300)
l2 = l2 mod 360 - 360*(l2 gt 300)
lowhigh, l1, l2
lowhigh, b1, b2
lowhigh, v1, v2
if l1 lt -10.5 or l2 gt 250.5 or b1 lt -7.5 or b2 gt 7.5 then begin
  print,'Galactic coordinate is out of the range of DLH survey!'
  print,'l = 350 ~ 0 ~ 250 deg, b = -5 ~ 5 deg'
  return
endif

nx = long((l2-l1)*120) + 1
cx = long(l2*120) + 1
ny = long((b2-b1)*120) + 1
cy = -long(b1*120) + 1

fitsl1 = floor(l1*2)/2.;+0.5*((l1*2-ceil(l1*2)) gt 0.25)
fitsl2 = ceil(l2*2)/2.
fitsb1 = floor(b1*2)/2.
fitsb2 = ceil(b2*2)/2.
print,'Mosaic begins:'
print,'GL from '+string(l1 mod 360)+' to '+string(l2 mod 360)
print,'GB from '+string(b1)+' to '+string(b2)
if display then erase

;mosaic
count=0l
num=0l
fitsnum=string(((fitsl2-fitsl1)/0.5+1)*((fitsb2-fitsb1)/0.5+1), format='(I0)')
for l=double(fitsl2),fitsl1,-0.5 do begin
for b=double(fitsb1),fitsb2,0.5 do begin
   count += 1
   name = getcellname(l,b)
   ;search for cube file
   cubcan = path+name+sb+'.fits'	;candidate
   idx = where(file_test(cubcan), ccan)
   if ccan eq 0 then begin
      print, '['+string(count,format='(I0)')+'/'+fitsnum $
      +']File "'+name+sb+'.fits" does not exist!'
      continue
   endif
   cubfile = cubcan[idx[0]]
   ;search for rms file
   rmscan = path+name+sb+'_rms.fits'	;candidate
   idx = where(file_test(rmscan), ccan)   
   if ccan eq 0 then begin
      print, '['+string(count,format='(I0)')+'/'+fitsnum $
      +']RMS file for "'+name+sb+'.fits" does not exist'
      continue
      ;fits_read,'rmsmodel.fits',wei
   endif
   rmsfile = rmscan[idx[0]]
   ;read cube files 
   print, '['+string(count,format='(I0)')+'/'+fitsnum $
   +']Mosaic file "'+cubfile
   fits_read,cubfile,dat,hdr
   if sxpar(hdr,'BITPIX') gt 0 then nan = where(dat eq max(dat), L64) $
   else nan = where(dat eq -1000, L64)
   if nan[0] ne -1 then dat[nan] = 0
   ;read rms file
   fits_read,rmsfile,rms
   if sxpar(hdr,'BITPIX') gt 0 then nan=where(rms eq max(rms)) $
   else nan=where(rms eq -1000 or ~finite(rms))
   if nan[0] ne -1 then rms[nan] = !values.d_infinity
   wei = 1d/rms^2

   num += 1
   clipv, dat, hdr, v1, v2
   if hdr[0] eq 'OUT' then begin
      print,'Velocity is out of range!'
      return
   endif
   if num eq 1 then begin
      nv=sxpar(hdr,'NAXIS3')
      mdat=dblarr(nx,ny,nv)
      mwei=dblarr(nx,ny)
      mcov=bytarr(nx,ny)
      mcov[*]=0
      mhdr=hdr
   endif
   if sxpar(hdr,'NAXIS3') ne nv or sxpar(hdr,'CRVAL3') ne sxpar(mhdr,'CRVAL3') then begin
      print,'Inconsistent velosity dimension!'
      return
   endif
   for i=0l,nv-1 do dat[*,*,i]=dat[*,*,i]*wei
   fcx = cx-long(l*120)-1 ;idl 0-based
   fcy = cy+long(b*120)-1
   fnx = sxpar(hdr,'NAXIS1')
   fny = sxpar(hdr,'NAXIS2')
   x1=fcx-(fnx-1)/2 & x2=fcx+(fnx-1)/2	;subscript on the mosaic array, within x1:x2
   y1=fcy-(fny-1)/2 & y2=fcy+(fny-1)/2
   if x2 lt 0 or x1 gt (nx-1) or y2 lt 0 or y1 gt (ny-1) then continue
   xcut1=(x1 lt 0)*(-x1) & xcut2=fnx-1-(x2 gt nx-1)*(x2-nx+1)
   ycut1=(y1 lt 0)*(-y1) & ycut2=fny-1-(y2 gt ny-1)*(y2-ny+1)
   dat=dat[xcut1:xcut2,ycut1:ycut2,*]
   wei=wei[xcut1:xcut2,ycut1:ycut2]
   mdat[max([x1,0]):min([x2,nx-1]),max([y1,0]):min([y2,ny-1]),*] += temporary(dat)
   mwei[max([x1,0]):min([x2,nx-1]),max([y1,0]):min([y2,ny-1])] += temporary(wei)
   if display then tvscl,mwei
   fnx=61
   fny=61
   x1=fcx-(fnx-1)/2 & x2=fcx+(fnx-1)/2
   y1=fcy-(fny-1)/2 & y2=fcy+(fny-1)/2
   if x2 lt 0 or x1 gt (nx-1) or y2 lt 0 or y1 gt (ny-1) then continue
   xcut1=(x1 lt 0)*(-x1) & xcut2=fnx-1-(x2 gt nx-1)*(x2-nx+1)
   ycut1=(y1 lt 0)*(-y1) & ycut2=fny-1-(y2 gt ny-1)*(y2-ny+1)
   mcov[max([x1,0]):min([x2,nx-1]),max([y1,0]):min([y2,ny-1]),*]=1
endfor
endfor
if num eq 0 then begin
   print,'No file in the selected region!'
   return
endif
for i=0l,nv-1 do mdat[*,*,i]=mdat[*,*,i]/mwei
sxaddpar,mhdr,'NAXIS1',nx
sxaddpar,mhdr,'NAXIS2',ny
sxaddpar,mhdr,'CTYPE1','GLON-CAR'
sxaddpar,mhdr,'CTYPE2','GLAT-CAR'
sxaddpar,mhdr,'CRVAL1',0d
sxaddpar,mhdr,'CRVAL2',0d
sxaddpar,mhdr,'CRPIX1',(cx gt 21601l)?(cx-43200l):cx
sxaddpar,mhdr,'CRPIX2',cy
sxdelpar,mhdr,'OBJECT'
sxdelpar,mhdr,'GLAT'
sxdelpar,mhdr,'GLON'
sxaddpar,mhdr,'BUNIT','K (T_MB)'
sxaddhist,'MOSAIC: '+systime(0),mhdr
sxaddhist,'Contain '+string(num,format='(I0)')+' cells.',mhdr
fits_write,'mosaic_'+sb+'.fits',mdat,mhdr
fits_write,'mosaic_'+sb+'_weight.fits',mwei,mhdr
fits_write,'mosaic_'+sb+'_coverage.fits',mcov,mhdr
mrms=1/sqrt(mwei)
nan=where(~finite(mrms))
if nan[0] ne -1 then mrms[nan]=!values.d_nan
fits_write,'mosaic_'+sb+'_rms.fits',mrms,mhdr
print,'Mosaic succesfully'
end
