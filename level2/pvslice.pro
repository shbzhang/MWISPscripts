;pvslice, pvshow
;by ShaoboZhang
;History:
;Apr,26,2012,v1.0
;May,16,2012,v1.1
;  PVSHOW: delete vrange keyword, add keyword: line,transpose,arcmin,log_scale.
;          now can plot l-v or b-v map created by cubemoment.
;  PVSLICE: move the reference pixel of position to the center of the axis
;Dec,01,2015,v1.2
;  PVSLICE: encounter CTYPE error while calling EXTAST in astrolib,
;	    change CTYPE from 'VELOCITY' to 'VEL', from 'POSITION' to 'POS'
;Jan,04,2016,v2.0
;  PVSLICE: rewrite the procedure to accept arbitary path with WIDTH, correct the path resampling error


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





pro pvslice, fitsfile, a, d, width=width, step=step, gal=gal, spline=spline, kernel=kernel, ds9=ds9
;Extract position-velocity map from a data cube
;Only accept fits in celestial and Galactic system with CDELT1=CDELT2
;Input:
;	fitsfile: a string indicates the file name of datacube
;	a:	if 'd' is not present, 'a' is a string indicates the file name of the path, in "x y" format
;	a,d:	coordinate of the path, a=[x1,x2,x3,...], d=[y1,y2,x3,...]
;Input keyword:
;	width:	the width (in pixel unit) of the belt to be averaged, default is 1
;	step:	resample the input slice path with interval of step (in pixel unit), default is 0.5
;	gal:	set to 1 if the given a, d are in galactic system
;	spline:	set to 1 to do spline interpolation to the input path, the path will become smooth spline
;	kernel: string of user defined function kernel=f(x,y), normalization is not need.
;		normally, the pixels on pvslice image are bilinear interpolation of the cube
;		when a function name is given to this keyword, the pixels will be the convolution result of cube and kernel
;		at the end of this file, there is an example of gaussian kernel function.
;		WARNING - using kernel keyword will be very time-consuming.
;Output:
;	"pvslice.fits": pv map
;	"pvslice.path": the resampled path
;	"pvslice.belt": the outline of the belt
;Usage: pvslice, 'XXX.fits','path.cat',/gal
;	pvslice, 'XXX.fits',[x1,x2,x3],[y1,y2,y3],width=5,/spline
;Todo:	give position with emission more weight when averaging?
;	multi-kernel?
;	

;regular input check
if n_params() le 1 then begin
	print, 'Syntax - PVSLICE, fitsfile, catalog, [width=, step=, /gal, /spline]'
	print, 'Syntax - PVSLICE, fitsfile, a, d, [width=, step=, /gal, /spline]'
	return
endif
if ~file_test(fitsfile) then begin
	print,'Error: fits file does not exist!'
	return
end

;read fits
fits_read, fitsfile, dat, hdr

;read catalog
if n_params() eq 2 then begin
	catalog = a
	readcol, catalog, a, d
endif

;convert coordinate
if keyword_set(gal) then pathgal = 1b else pathgal = 0b	;path in gal or not?
fitsgal = strcmp(sxpar(hdr,'CTYPE1'), 'GL', 2, /fold_case)	;
if fitsgal eq pathgal then at=a & dt=d
if fitsgal and ~pathgal then glactc,a,d,2000,at,dt,1,/deg
if ~fitsgal and pathgal then glactc,at,dt,2000,a,d,2,/deg
adxy,hdr,at,dt,x,y	;convert wcs path to xy path on datacube

;convert polyline to spline
if keyword_set(spline) then spline_p,x,y,x,y	;if a spline resample is need

;resample path with interval of step
defstep=0.5	;unit pixel
if ~keyword_set(step) then step = defstep
if step le 0 then step = defstep
x = double(x) & y = double(y)
nodnum = min([n_elements(x),n_elements(y)])

px=x[0]
py=y[0]	;xy point on path
nod = 1	;xy of nod on the segment
path = [px, py]
while nod le nodnum-1 do begin
rest = sqrt((x[nod]-px)^2+(y[nod]-py)^2)
if rest ge step then begin ;the rest of current segment is enough for another step
	px += (x[nod]-px)/rest*step
	py += (y[nod]-py)/rest*step
	path = [path,px,py]
endif else begin	;not enough, find from the next segment or the next next one...except the last nod
	if nod eq nodnum-1 then begin	;the last nod, extend a little
		last = sqrt((x[nod]-px)^2+(y[nod]-py)^2)
		if last gt step/1e3 then begin
			px+=(x[nod]-px)/last*step
			py+=(y[nod]-py)/last*step
		path = [path, px, py]
		endif
		break
	endif
	while nod lt nodnum-1 do begin
		ptop2 = sqrt((x[nod+1]-px)^2 + (y[nod+1]-py)^2)
		if ptop2 lt step then nod++ else begin
			p1top2 = sqrt((x[nod+1]-x[nod])^2 + (y[nod+1]-y[nod])^2)
			dx21 = x[nod+1] - x[nod]
			dy21 = y[nod+1] - y[nod]
			crossx = x[nod]*y[nod+1]*dy21 - x[nod+1]*y[nod]*dy21 + dx21^2*px + dx21*dy21*py
			crossx /= dx21^2+dy21^2
			crossy = x[nod+1]*y[nod]*dx21 - x[nod]*y[nod+1]*dx21 + dx21*dy21*px + dy21^2*py
			crossy /= dx21^2+dy21^2
			crosstonext = sqrt((step^2 -(px-crossx)^2 -(py-crossy)^2) >0)
			px = crossx + dx21/p1top2*crosstonext
			py = crossy + dy21/p1top2*crosstonext
			path = [path,px,py]
			nod++
			break
		endelse
	endwhile
endelse
endwhile
path = reform(path,2,n_elements(path)/2)
pathx = (path[0,*])[*]
pathy = (path[1,*])[*]

;export the path
xyad,hdr,pathx,pathy,patha,pathd
if fitsgal and ~pathgal then glactc,patha,pathd,2000,patha,pathd,2,/deg
if ~fitsgal and pathgal then glactc,patha,pathd,2000,patha,pathd,1,/deg
openw,lun,'pvslice.path',/get_lun
if ~keyword_set(ds9) then $
	printf,lun,transpose(string(patha)+' '+string(pathd)) $
else begin
	openw,lun2,'pvslice2.path',/get_lun
	ds9color=['White','Red','Yellow','Green','Cyan','Blue','Magenta','Black']
	printf,lun,'wcs;'
	printf,lun2,'image;'
	for i=0,n_elements(patha)-2 do begin
		printf,lun,'line('+string(patha[i])+','+string(pathd[i])+','+string(patha[i+1])+','+string(pathd[i+1])+ $
				') # width=2 color='+ds9color[i/10 mod n_elements(ds9color)]
		printf,lun2,'line(1,'+string(i+1)+',1,'+string(i+2)+ $
				') # width=4 color='+ds9color[i/10 mod n_elements(ds9color)]
	endfor
	close,lun2
	free_lun,lun2
endelse
close,lun
free_lun,lun

;expand path to belt
stepnum = n_elements(pathx)
if keyword_set(width) then width=fix(width) else width=1
if width gt 1 then begin
	dy = shift(pathy, 1) - shift(pathy, -1)
	dy[0] = pathy[0]-pathy[1]
	dy[stepnum-1] = pathy[stepnum-2]-pathy[stepnum-1]
	dx = shift(pathx, 1) - shift(pathx, -1)
	dx[0] = pathx[0]-pathx[1]
	dx[stepnum-1] = pathx[stepnum-2]-pathx[stepnum-1]
	pa = atan(dy/dx)+(dx lt 0)*!pi+!pi/2
	dx = cos(pa)
	dy = sin(pa)
	pathx = pathx#replicate(1,width)+dx#(findgen(width)-(width-1)/2.)
	pathy = pathy#replicate(1,width)+dy#(findgen(width)-(width-1)/2.)
endif else width=1

;export the outline of the belt
xyad,hdr,pathx,pathy,patha,pathd
if fitsgal and ~pathgal then glactc,patha,pathd,2000,patha,pathd,2,/deg
if ~fitsgal and pathgal then glactc,patha,pathd,2000,patha,pathd,1,/deg
ola = [patha[*,0], reverse(patha[*,width-1]), patha[0,0]]
old = [pathd[*,0], reverse(pathd[*,width-1]), pathd[0,0]]
openw,lun,'pvslice.belt',/get_lun
printf,lun,transpose(string(ola)+' '+string(old))
close,lun
free_lun,lun

;extract value from cube
channum = sxpar(hdr, 'NAXIS3')
slice = make_array(stepnum, width, channum, type=size(dat,/type))
if ~keyword_set(kernel) then $
	for i=0,channum-1 do slice[*,*,i] = interpolate(dat[*,*,i],pathx,pathy) $
else	for i=0,channum-1 do begin
		print,'convolving: ',i+1,'/',channum
		slice[*,*,i] = convolve(dat[*,*,i],pathx,pathy,kernel)
	endfor
slice = mean(slice, dimension=2, /nan)
;slice = total(slice,2,/nan)/total(finite(slice),2)	;for idl7
slice = transpose(slice)

;output slice and resampled path
mkhdr,pvhdr,slice
sxaddpar,pvhdr,'CTYPE1','V'
sxaddpar,pvhdr,'CRPIX1',sxpar(hdr,'CRPIX3')
sxaddpar,pvhdr,'CRVAL1',sxpar(hdr,'CRVAL3')/1e3	;convert m/s to km/s
sxaddpar,pvhdr,'CDELT1',sxpar(hdr,'CDELT3')/1e3
sxaddpar,pvhdr,'CTYPE2','P'
sxaddpar,pvhdr,'CRPIX2',1
sxaddpar,pvhdr,'CRVAL2',0d
sxaddpar,pvhdr,'CDELT2',step*abs(sxpar(hdr,'CDELT1'))
sxaddhist,'PV file: '+fitsfile,pvhdr
sxaddhist,'PV path:',pvhdr
for i=0,n_elements(a)-1 do sxaddhist,string(a[i])+' '+string(d[i]),pvhdr
sxaddhist,'PV width: '+string(width)+' pixels',pvhdr
sxaddhist,'Position in Degree',pvhdr
sxaddhist,'Velocity in km/s',pvhdr
fits_write,'pvslice.fits',slice,pvhdr

end





pro pvcircle, fitsfile, a0, d0, radius, width=width, step=step, gal=gal
;a circle slice
if n_params() le 3 then begin
	print, 'Syntax - PVCIRCLE, fitsfile, a0, d0, radius, [width=, step=, /gal]'
	return
endif
fits_read,fitsfile,dat,hdr

if keyword_set(gal) then pathgal = 1b else pathgal = 0b	;path in gal or not?
fitsgal = strcmp(sxpar(hdr,'CTYPE1'), 'GL', 2, /fold_case)	;
if fitsgal eq pathgal then at=a0 & dt=d0
if fitsgal and ~pathgal then glactc,a0,d0,2000,at,dt,1,/deg
if ~fitsgal and pathgal then glactc,at,dt,2000,a0,d0,2,/deg
adxy,hdr,at,dt,x0,y0	;convert wcs path to xy path on datacube

defstep=0.5	;unit pixel
if ~keyword_set(step) then step = defstep
if step le 0 then step = defstep

;get path
ang = 2*asin(step/2./radius)
stepnum = round(2*!dpi/ang)
ang = 2*!dpi/stepnum
step = 2*radius*sin(!dpi/stepnum)
pa = 2*!dpi/stepnum*dindgen(stepnum+1)	;PA
pathx = x0-radius*sin(pa)
pathy = y0+radius*cos(pa)

;expand to belt
if keyword_set(width) then width=fix(width) else width=1
if width gt 1 then begin
	dx = -sin(pa)
	dy = cos(pa)
	pathx = pathx#replicate(1,width)+dx#(findgen(width)-(width-1)/2.)
	pathy = pathy#replicate(1,width)+dy#(findgen(width)-(width-1)/2.)
endif else width=1

;extract value
channum = sxpar(hdr, 'NAXIS3')
slice = make_array(stepnum+1, width, channum, type=size(dat,/type))
for i=0,channum-1 do slice[*,*,i] = interpolate(dat[*,*,i],pathx,pathy,missing=0)
slice = mean(slice, dimension=2, /nan)
slice = transpose(slice)

;output slice and resampled path
mkhdr,pvhdr,slice
sxaddpar,pvhdr,'CTYPE1','V'
sxaddpar,pvhdr,'CRPIX1',sxpar(hdr,'CRPIX3')
sxaddpar,pvhdr,'CRVAL1',sxpar(hdr,'CRVAL3')/1e3	;convert m/s to km/s
sxaddpar,pvhdr,'CDELT1',sxpar(hdr,'CDELT3')/1e3
sxaddpar,pvhdr,'CTYPE2','PA'
sxaddpar,pvhdr,'CRPIX2',1
sxaddpar,pvhdr,'CRVAL2',0d
sxaddpar,pvhdr,'CDELT2',ang/!dtor
sxaddhist,'PV file: '+fitsfile,pvhdr
sxaddhist,'PVCIRCLE:'+string(a0)+' '+string(d0),pvhdr
sxaddhist,'PV radius:'+string(radius),pvhdr
sxaddhist,'PV width: '+string(width)+' pixels',pvhdr
sxaddhist,'Velocity in km/s',pvhdr
fits_write,'pvslice.fits',slice,pvhdr

end





pro pathtest, x,y,pathx,pathy,width,step
x = double(x) & y = double(y)
nodnum = min([n_elements(x),n_elements(y)])

px=x[0]
py=y[0]	;xy point on path
nod = 1	;xy of nod on the segment
path = [px, py]
while nod le nodnum-1 do begin
rest = sqrt((x[nod]-px)^2+(y[nod]-py)^2)
if rest ge step then begin ;the rest of current segment is enough for another step
	px += (x[nod]-px)/rest*step
	py += (y[nod]-py)/rest*step
	path = [path,px,py]
endif else begin	;not enough, find from the next segment or the next next one...except the last nod
	if nod eq nodnum-1 then begin	;the last nod, extend a little
		last = sqrt((x[nod]-px)^2+(y[nod]-py)^2)
		if last gt step/1e3 then begin
			px+=(x[nod]-px)/last*step
			py+=(y[nod]-py)/last*step
		path = [path, px, py]
		endif
		break
	endif
	while nod lt nodnum-1 do begin
		ptop2 = sqrt((x[nod+1]-px)^2 + (y[nod+1]-py)^2)
		if ptop2 lt step then nod++ else begin
			p1top2 = sqrt((x[nod+1]-x[nod])^2 + (y[nod+1]-y[nod])^2)
			dx21 = x[nod+1] - x[nod]
			dy21 = y[nod+1] - y[nod]
			crossx = x[nod]*y[nod+1]*dy21 - x[nod+1]*y[nod]*dy21 + dx21^2*px + dx21*dy21*py
			crossx /= dx21^2+dy21^2
			crossy = x[nod+1]*y[nod]*dx21 - x[nod]*y[nod+1]*dx21 + dx21*dy21*px + dy21^2*py
			crossy /= dx21^2+dy21^2
			crosstonext = sqrt((step^2 -(px-crossx)^2 -(py-crossy)^2) >0)
			px = crossx + dx21/p1top2*crosstonext
			py = crossy + dy21/p1top2*crosstonext
			path = [path,px,py]
			nod++
			break
		endelse
		;edge1=sqrt((x[nod]-px)^2+(y[nod]-py)^2)
		;edge2=sqrt((x[nod+1]-px)^2+(y[nod+1]-py)^2)
		;;if edge2 ge step then nod ++ continue
		;edge3=sqrt((x[nod+1]-x[nod])^2+(y[nod+1]-y[nod])^2)
		;a1 = acos(((edge1^2+edge3^2-edge2^2)/(2*edge1*edge3)) >(-1) < 1)
		;a2 = asin(edge1/step*sin(a1) >(-1) <1)
		;len = step/sin(a1)*sin(!pi-a1-a2)
		;tpx = x[nod]+(x[nod+1]-x[nod])/edge3*len
		;tpy = y[nod]+(y[nod+1]-y[nod])/edge3*len
		;stop
		;if tpx gt max(x[nod:nod+1]) or tpx lt min(x[nod:nod+1]) then nod++ else begin
		;	px = tpx & py = tpy
		;	path = [path, px, py]
		;	nod++
		;	break
		;endelse
	endwhile
endelse
endwhile
path = reform(path,2,n_elements(path)/2)
pathx = (path[0,*])[*]
pathy = (path[1,*])[*]

;expand path to belt
stepnum = n_elements(pathx)
if keyword_set(width) then width=fix(width) else width=1
if width gt 1 then begin
	dy = shift(pathy, 1) - shift(pathy, -1)
	dy[0] = pathy[0]-pathy[1]
	dy[stepnum-1] = pathy[stepnum-2]-pathy[stepnum-1]
	dx = shift(pathx, 1) - shift(pathx, -1)
	dx[0] = pathx[0]-pathx[1]
	dx[stepnum-1] = pathx[stepnum-2]-pathx[stepnum-1]
	pa = atan(dy/dx)+(dx lt 0)*!pi+!pi/2
	dx = cos(pa)
	dy = sin(pa)
	pathx = pathx#replicate(1,width)+dx#(findgen(width)-(width-1)/2.)
	pathy = pathy#replicate(1,width)+dy#(findgen(width)-(width-1)/2.)
endif else width=1
end


function convolve, img, x, y, kfunction
;x,y must have the same size
isz = size(img, /dim)
osz = size(x, /dim)

o = make_array(osz, type=size(img,/type))

gx = findgen(isz[0])#replicate(1,isz[1])
gy = replicate(1,isz[0])#findgen(isz[1])
for i=0, n_elements(x)-1 do begin
	k = call_function(kfunction, gx-x[i], gy-y[i])
	k /= total(k)
	o[i] = total(img * k, /nan)
endfor
return,o
end

function gaussian, x, y
d=sqrt((x)^2+(y)^2)
f=exp(-d^2/(2.*(0.70776818d /2)^2))
return,f
end
