pro checkoffset, fitsfile
;smooth the fits file of a cell to ~1km/s
if n_params() lt 1 then begin
	print,'Syntax - CHECKOFFSET, fitsfile'
	return
endif
fits_read,fitsfile,dat,hdr
off = fltarr(sxpar(hdr,'NAXIS1'), sxpar(hdr,'NAXIS2'), 401)
start = sxpar(hdr,'CRPIX3')-1-3-200*7
for i=0,400 do off[*,*,i] = total(dat[*,*,(start+i*7):(start+i*7+6)],3)/7
sxaddpar,hdr,'NAXIS3',401
sxaddpar,hdr,'CRPIX3',201
sxaddpar,hdr,'CDELT3',sxpar(hdr,'CDELT3')*7
fits_write,file_basename(fitsfile,'.fits')+'_off.fits',off,hdr

;collapse offset cube
area = total(off gt 0.5, 3) gt 0
area += total(off gt 1.0, 3) gt 0
area += total(off gt 1.5, 3) gt 0
area += total(off gt 2.0, 3) gt 0
area /= 2.
fits_write,file_basename(fitsfile,'.fits')+'_offarea.fits',area,hdr

;find possible offset in area map
x=findgen(11)#replicate(1,11)-5
y=findgen(11)##replicate(1,11)-5
c=sqrt(x^2+y^2) le 5
mask = convol(area ge 1, c) gt 0
fits_write,file_basename(fitsfile,'.fits')+'_offposs.fits',mask,hdr
end

pro beamposition,gl,gb,y,m,d,lst
;display beams in ds9
jdcnv,y,m,d,12-8,jd      ;LCT -> UTC -> JULDAY
ct2lst,lst12,ten(97,33.6),8,jd,d,m,y
dlst=lst-lst12
jdcnv,y,m,d,12-8+dlst,jd      ;LCT -> UTC -> JULDAY
x=[[-172.1, 0, 175.8],$
   [-175.1, 0, 175.4],$
   [-177.0, 0, 173.4]]
y=[[177.3, 175.6, 174.5],$
   [0, 0, 0],$
   [-175.1, -177.2, -178.3]]
s=[[58.6, 59.2, 59.7],$
   [57.8, 59.8, 60.2],$
   [58.4, 59.8, 60.0]]
glactc,ra,dec,2000,gl,gb,2,/degree
eq2hor,ra,dec,jd,el,az,lat=ten(37,22.4),lon=ten(97,33.6),altitude=3200d
az += x/3600d/cos(el*!dtor)
el += y/3600d
hor2eq,el,az,jd,ra,dec,lat=ten(37,22.4),lon=ten(97,33.6),altitude=3200d
glactc,ra,dec,2000,l,b,1,/degree
print,'PA = ',atan((l[1]-l[4])/(b[1]-b[4]))/!dtor,'positive means clockwise'
print,'circle('+transpose(string(l[*],form='(f0)')+','+string(b[*],form='(f0)')+','+string(s[*],form='(f0)'))+'") # color=red' 
end
