;by ShaoboZhang
;Mosaic fits file of DLH survey
;Usage: mosaic, gl1, gl2, gb1, gb2, velo1, velo2 [, sb='L', path='./', /display]
;Input:
;  gl1,gl2,gb1,gb2,velo1,velo2: float scalar. Galactic coordinate range and velocity range (in km/s).
;Input keyword:
;  sb: sideband, for 12CO, sb='U'; for 13CO, sb='L'; for C18O, sb='L2'
;  path: array of the directories containing datacube and rms files, ended with '/'
;  dispaly: whether to display the process on a window. Default is 1.
;  output: prefix of file names
;  silent: whether to suppress message in the terminal. Default is 0
;Note:
;  all the input fits file (data and rms) must be the same as the grid define by template.
;History:
;Nov,09,2011,v1.0
;Nov,28,2011,v1.1
;  add the fitspath and rmspath keywords.
;Jan,16,2012,v1.2
;  fix a error when mosaic region with gl > 180
;Apr,12,2012,v1.3
;  change the output file name.
;  output a mask fits file for the mosaic.
;May,09,2012,v1.4
;  complete the header of the result fits file.
;  the mask file is now called 'coverage' file to fit its function.
;  add the keywords "display" for users who do not support window.
;Dec,25,2013,v1.5
;  accept rms file created by procedure cuberms.
;Oct,13,2015,v1.6
;  remove keyword "fitspath" and "rmspath".
;  new keyword "path" now accepts string array, procedure will
;    automatically search datacube and rms file in these directories.
;Dec,09,2016,v1.7
;  add history to indicate the completeness of the mosaic.
;  add keyword 'output' for the names of the output files.
;  add keyword 'silent' to suppress any message in terminal.
;  fix spelling mistakes.
;May,17,2018,v1.8
;  fix a minor error that certain input range does not match the output range.

function getcellname, gl, gb
;get cell name from its gl and gb
cl = round(gl*2)/2d
cb = round(gb*2)/2d
cl = cl mod 360
if cl lt 0 then cl += 360
return,string(fix(cl*10),format='(I04)')+string(fix(cb*10),format='(I+04)')
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
if (n[0] gt nc-1) or (n[1] lt 0) then begin
;        print,'Velocity out of range.'
	hdr='OUT'
        return
endif
n = n >0 <(nc-1)
;print,'Clip channel between:',n
dat = dat[*,*,n[0]:n[1]]
sxaddpar,hdr,'NAXIS3',n[1]-n[0]+1
sxaddpar,hdr,'CRPIX3',cp-n[0]
sxaddpar,hdr,'CDELT3',cd
end



pro mosaic, l1, l2, b1, b2, v1, v2, sb=sb, path=path, display=display, output=output, silent=silent

if n_params() lt 6 then begin
    print, 'Syntax - MOSAIC, l1, l2, b1, b2, v1, v2, [sb=, path=, /display, output=, /silent]'
    return
endif
if ~keyword_set(sb) then sb='U'
if ~keyword_set(path) then path=''
if n_elements(display) eq 0 then display=1
if ~keyword_set(output) then output='mosaic'
if ~keyword_set(silent) then silent=0b

l1 = l1 mod 360 - 360*(l1 gt 300)
l2 = l2 mod 360 - 360*(l2 gt 300)
lowhigh, l1, l2
lowhigh, b1, b2
lowhigh, v1, v2
if l1 lt -10.5 or l2 gt 250.5 or b1 lt -7.5 or b2 gt 7.5 then begin
  if ~silent then print,'Galactic coordinate is out of the range of DLH survey!'
  if ~silent then print,'l = 350 ~ 0 ~ 250 deg, b = -5 ~ 5 deg'
  return
endif

nx = floor(l2*120) - ceil(l1*120) + 1
cx = floor(l2*120) + 1
ny = floor(b2*120) - ceil(b1*120) + 1
cy = -ceil(b1*120) + 1

;nx = long((l2-l1)*120) + 1
;cx = long(l2*120) + 1
;ny = long((b2-b1)*120) + 1
;cy = -long(b1*120) + 1

fitsl1 = floor(l1*2)/2.;+0.5*((l1*2-ceil(l1*2)) gt 0.25)
fitsl2 = ceil(l2*2)/2.
fitsb1 = floor(b1*2)/2.
fitsb2 = ceil(b2*2)/2.
if ~silent then print,'Mosaic begins:'
if ~silent then print,'GL from '+string(l1 mod 360)+' to '+string(l2 mod 360)
if ~silent then print,'GB from '+string(b1)+' to '+string(b2)
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
   ;idx = where(file_test(cubcan), ccan)
   cubfile=[]
   for i=0,n_elements(cubcan)-1 do begin
      fs = file_search(cubcan[i],count=fc)
      if fc gt 0 then cubfile=[cubfile,fs]
   endfor
   if n_elements(cubfile) eq 0 then begin
      if ~silent then print, '['+string(count,format='(I0)')+'/'+fitsnum $
      +']File "'+name+sb+'.fits" does not exist!'
      continue
   endif
   cubfile = cubfile[0]

   ;search for rms file
   rmscan = path+name+sb+'_rms.fits'	;candidate
   ;idx = where(file_test(rmscan), ccan)   
   rmsfile=[]
   for i=0,n_elements(rmscan)-1 do begin
      fs = file_search(rmscan[i],count=fc)
      if fc gt 0 then rmsfile=[rmsfile,fs]
   endfor
   if n_elements(rmsfile) eq 0 then begin
      if ~silent then print, '['+string(count,format='(I0)')+'/'+fitsnum $
      +']RMS file for "'+name+sb+'.fits" does not exist'
      continue
      ;fits_read,'rmsmodel.fits',wei
   endif
   rmsfile = rmsfile[0]
   
   ;read cube files 
   if ~silent then print, '['+string(count,format='(I0)')+'/'+fitsnum $
   +']Mosaic file "'+cubfile+'"'
   fits_read,cubfile,dat,hdr
   if sxpar(hdr,'BITPIX') gt 0 then nan = where(dat eq max(dat) or ~finite(dat), /L64) $
   else nan = where(dat eq -1000 or ~finite(dat), /L64)
   if nan[0] ne -1 then dat[nan] = 0
  
   ;read rms file, get weight
   fits_read,rmsfile,rms
   if sxpar(hdr,'BITPIX') gt 0 then nan=where(rms eq max(rms) or ~finite(dat[*,*,8192])) $
   else nan=where(rms eq -1000 or ~finite(rms) or ~finite(dat[*,*,8192]))
   if nan[0] ne -1 then rms[nan] = !values.d_infinity
   wei = 1d/rms^2

   ;clip dat
   clipv, dat, hdr, v1, v2
   if hdr[0] eq 'OUT' then begin
      print,'Error - Velocity is out of the channel range!'
      return
   endif

   num += 1
   ;generate a mosaic cube and header
   if num eq 1 then begin
      nv=sxpar(hdr,'NAXIS3')
      mdat=dblarr(nx,ny,nv)
      mwei=dblarr(nx,ny)
      mcov=bytarr(nx,ny)
      mcov[*]=0
      mhdr=hdr
      mvran = ([0,sxpar(mhdr,'NAXIS3')-1]-sxpar(mhdr,'CRPIX3')+1)*sxpar(mhdr,'CDELT3')+sxpar(mhdr,'CRVAL3')
   endif

   vran = ([0,sxpar(hdr,'NAXIS3')-1]-sxpar(hdr,'CRPIX3')+1)*sxpar(hdr,'CDELT3')+sxpar(hdr,'CRVAL3')
   if sxpar(hdr,'NAXIS3') ne nv or abs(vran[0]-mvran[0]) gt sxpar(mhdr,'CDELT3')/1000 or abs(vran[1]-mvran[1]) gt sxpar(mhdr,'CDELT3')/1000 then begin	;velocity range error should be less than 0.1% of channel width
      print,'Error - Inconsistent velocity dimension!'
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
   print,'Error - No file in the selected region!'
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
sxaddhist,'L from '+string(l1)+' to '+string(l2),mhdr
sxaddhist,'B from '+string(b1)+' to '+string(b2),mhdr
sxaddhist,'V from '+string(v1)+' to '+string(v2),mhdr
sxaddhist,'Mosaic Completeness: '+'['+string(num,format='(I0)')+'/'+fitsnum+']',mhdr
;sxaddhist,'Contain '+string(num,format='(I0)')+' cells.',mhdr
;output
fits_write,output+'_'+sb+'.fits',mdat,mhdr
fits_write,output+'_'+sb+'_weight.fits',mwei,mhdr
fits_write,output+'_'+sb+'_coverage.fits',mcov,mhdr
mrms=1/sqrt(mwei)
nan=where(~finite(mrms))
if nan[0] ne -1 then mrms[nan]=!values.d_nan
fits_write,output+'_'+sb+'_rms.fits',mrms,mhdr
if ~silent then print,'Mosaic successfully'
end
