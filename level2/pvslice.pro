;pvslice, pvshow
;by ShaoboZhang
;History:
;Apr,26,2012,v1.0
;May,16,2012,v1.1
;  PVSHOW: delete vrange keyword, add keyword: line,transpose,arcmin,log_scale.
;          now can plot l-v or b-v map created by cubemoment.
;  PVSLICE: move the reference pixel of position to the center of the axis


pro pvslice, fitsfile, a, d, gal=gal, step=step;, center=
;Extract position-velocity map from a data cube
;Only accept celestial and Galactic system
;Input:
;	fitsfile: data cube file name
;	a,d:	coordinate, a=[x1,x2,x3,...], d=[y1,y2,x3,...]
;Input keyword:
;	gal:	set to 1 if the given a,d are in galactic system
;	step:	resample step in pixel unit, default is 0.5
;Usage: pvslice, 'XXX.fits',[x1,x2],[y1,y2],/gal
if n_params() lt 3 then begin
    print, 'Syntax - PVSLICE, fitsfile, a, d, [/gal, step= ]'
    return
endif
if ~file_test(fitsfile) then begin
	print,'Error: Fits file not exist!'
	return
end
fits_read,fitsfile,dat,hdr
;extast,hdr,ast

if keyword_set(gal) then gal = 1b else gal = 0b
fitsgal = strcmp(sxpar(hdr,'CTYPE1'), 'GL', 2, /fold_case)
if fitsgal eq gal then at=a & dt=d
if fitsgal and ~gal then glactc,a,d,2000,at,dt,1,/deg
if ~fitsgal and gal then glactc,at,dt,2000,a,d,2,/deg
sxaddpar,hdr,'CTYPE1',repstr(sxpar(hdr,'CTYPE1'),'GLS','SFL')
sxaddpar,hdr,'CTYPE2',repstr(sxpar(hdr,'CTYPE2'),'GLS','SFL')
adxy,hdr,at,dt,x,y

defstep=0.5	;unit pixel
if ~keyword_set(step) then step = defstep
if step le 0 then step = defstep
length = sqrt( ((x-shift(x,1))[1:*])^2+((y-shift(y,1))[1:*])^2 )
nstep = ceil(total(length)/step)
step = total(length)/nstep
print,'Use step '+string(step)+' * pixelsize'
xs = x[0]+(x[1]-x[0])*findgen(nstep+1)/nstep
ys = y[0]+(y[1]-y[0])*findgen(nstep+1)/nstep

intlen = length
for i=n_elements(length)-1,0,-1 do intlen[i]=total(length[0:i])
;xs = fltarr(nstep+1)
;ys = fltarr(nstep+1)
for i=0, nstep do begin
	way = i*step
	node = (where(way le intlen))[0]
	xs[i] = x[node+1]-(intlen[node]-way)/length[node]*(x[node+1]-x[node])
	ys[i] = y[node+1]-(intlen[node]-way)/length[node]*(y[node+1]-y[node])
endfor

nx1 = sxpar(hdr,'NAXIS3')
nx2 = n_elements(xs)
slice = make_array(nx1, nx2, type=size(dat,/type))
for i=0,nx1-1 do slice[i,*] = interpolate(dat[*,*,i],xs,ys,missing=0)
mkhdr,pvhdr,slice
sxaddpar,pvhdr,'CTYPE1','VELOCITY'
sxaddpar,pvhdr,'CRPIX1',sxpar(hdr,'CRPIX3')
sxaddpar,pvhdr,'CRVAL1',sxpar(hdr,'CRVAL3')/1000d
sxaddpar,pvhdr,'CDELT1',sxpar(hdr,'CDELT3')/1000d
sxaddpar,pvhdr,'CTYPE2','POSITION'
sxaddpar,pvhdr,'CRPIX2',1
sxaddpar,pvhdr,'CRVAL2',0d
sxaddpar,pvhdr,'CDELT2',step*abs(sxpar(hdr,'CDELT1'))
sxaddhist,'PV file: '+fitsfile,pvhdr
sxaddhist,'PV path:',pvhdr
for i=0,n_elements(a)-1 do sxaddhist,string(a[i])+' '+string(d[i]),pvhdr
sxaddhist,'Position in Degree',pvhdr
sxaddhist,'Velocity in km/s',pvhdr
fits_write,'pvslice.fits',slice,pvhdr

xyad,hdr,xs,ys,as,ds
;velo = sxpar(pvhdr,'CRVAL1')+(dindgen(sxpar(pvhdr,'NAXIS1'))-sxpar(pvhdr,'CRPIX1')+1)*sxpar(pvhdr,'CDELT1')
posi = sxpar(pvhdr,'CRVAL2')+(dindgen(sxpar(pvhdr,'NAXIS2'))-sxpar(pvhdr,'CRPIX2')+1)*sxpar(pvhdr,'CDELT2')
;mmnt =  make_array(nx2, type=size(dat,/type))
;for i=0,nx2-1 do mmnt[i] = total(slice[*,i]*velo)/total(slice[*,i])
openw,lun,'pvslice.track',/get_lun
printf,lun,string(lindgen(n_elements(as))+1)+' '+string(as)+' '+string(ds)+' '+string(posi);+' '+string(mmnt)
close,lun
free_lun,lun
end


pro pvshow, pvfile, levels=levels, range=range, grey=grey, line=line, $
	transpose=transpose, arcmin=arcmin, log_scale=log_scale, _extra=_extra
;show a pvslice fits
;Input:
;	pvfile: file created by pvslice.pro
;Input keyword:
;	grey: set 1 to draw a grey background, default is 0
;	line: set 0 to not draw the contour line, default is 1
;	range: level range, default 10 levels. eg. range=[1,10]
;	levels: contour levels. eg. levels=[1,3,5,7,9]
;	transpose: transpose the velocity and position axis
;	arcmin: set 1 to change position unit to arcmin, default is 0 (deg)
;	log_scale: scale the data with log
;	accept other plot keywords
;Usage:
;	pvshow,'pvslice.fits',xrange=[-200,200],yrange=[0,1],range=[1,10]
;	pvshow,'pvslice.fits',xrange=[0,60],yrange=[-200,200],levels=findgen(10)*2+1,/grey,line=0,/transpose,/arcmin
if n_params() lt 1 then begin
    print, 'Syntax - PVSHOW, pvfile, [levels=, range=, /line, /grey, /transpose, /arcmin, /logs_cale ]'
    return
endif
if ~file_test(pvfile) then begin
	print,'PVslice file not exist!'
	return
end
if keyword_set(transpose) then transpose=1b else transpose=0b
if keyword_set(arcmin) then begin
  unit=" ' "
  factor=60d
endif else begin
  unit=' !Uo!N '
  factor=1d
endelse
fits_read,pvfile,slice,hdr
if keyword_set(log_scale) then slice = alog(slice)
axis1 = sxpar(hdr,'CRVAL1')+(dindgen(sxpar(hdr,'NAXIS1'))-sxpar(hdr,'CRPIX1')+1)*sxpar(hdr,'CDELT1')
axis2 = sxpar(hdr,'CRVAL2')+(dindgen(sxpar(hdr,'NAXIS2'))-sxpar(hdr,'CRPIX2')+1)*sxpar(hdr,'CDELT2')
title1 = strtrim(gettok(sxpar(hdr,'CTYPE1'),'-'),2)
title2 = strtrim(gettok(sxpar(hdr,'CTYPE2'),'-'),2)
;find velocity axis
axisv = where(strcmp([title1,title2], 'v', 1, /fold_case))

case axisv of	;change unit and title
  0:begin
;  axis1=axis1/1000d	;v to km/
  title1=title1+' (km/s)'
  axis2=axis2*factor
  title2=title2+' ('+unit+')'
  end
  1:begin
;  axis2=axis2/1000d	;v to km/s
  title2=title2+' (km/s)'
  axis1=axis1*factor
  title1=title1+' ('+unit+')'
  end
endcase
if transpose then begin	;transpose the axis
  slice=transpose(slice)
  t=axis1
  axis1=axis2
  axis2=t
  t=title1
  title1=title2
  title2=t
endif

if n_elements(levels) lt 1 then begin
  if n_elements(range) lt 2 then begin
    high=max(slice[where(finite(slice))], min=low)
    range=[low,high]
  endif else range=range(sort(range))
  levels= range[0]+(range[1]-range[0])*(findgen(10))/9
endif
levels=[-1000,levels[sort(levels)]]
levels=levels[uniq(levels)]
if n_elements(grey) lt 1 then grey=0
if n_elements(line) lt 1 then line=1
color=findgen(n_elements(levels))/(n_elements(levels)-1)*200
if !p.color eq 0 then color = 255-color	;for white background device, eg. PS

loadct,0
plot,[min(axis1),max(axis1)],[min(axis2),max(axis2)],/nodata,/ynozero, _extra=_extra
print,[min(axis1),max(axis1)],[min(axis2),max(axis2)]
loadct,0
case grey of
1:begin
  nan = where(finite(slice,/nan))
  if nan[0] ne -1 then slice[nan]=-2000
  contour, slice, axis1, axis2, levels=levels, /overplot, /fill, c_color=color
  if nan[0] ne -1 then slice[nan]=!values.f_nan
end
2:begin
  glevel = [-1000,range[0]+(range[1]-range[0])*(findgen(200))/199]
  gcolor = findgen(n_elements(levels))/(n_elements(levels)-1)*200
  nan = where(finite(slice,/nan))
  if nan[0] ne -1 then slice[nan]=-2000
  contour,gslice,gaxis1,gaxis2, levels=levels, /overplot, /fill, c_color=color
end
else:
endcase

if line then contour, slice, axis1, axis2, levels=levels, /overplot
loadct,0
plot,[min(axis1),max(axis1)],[min(axis2),max(axis2)],/nodata, /ynozero, /noerase, xtitle=title1, ytitle=title2, _extra=_extra
end
