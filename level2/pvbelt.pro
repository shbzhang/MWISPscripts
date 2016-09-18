pro pvbelt, fitsfile, a, d, width, gal=gal, step=step, pcoordiante=pcoordinate
;OBSOLETE, replaced by new version of PVSLICE
;Collapse position-velocity map of a belt from a data cube
;Only accept celestial and Galactic system
;Input:
;       fitsfile: data cube file name
;       a,d: 2-element vectors indicate coordinate, a=[x1,x2], d=[y1,y2]
;	width: belt width, in pixel
;Input keyword:
;       gal:    set to 1 if the given a,d are in galactic system
;       step:   resample step in pixel unit, default is 0.5
;	pcoordinate: whether keep position info in the cube. 0=use degree, 1=use CTYPE1, 2=use CTYPE2
;Usage: pvbelt, 'XXX.fits',[x1,x2],[y1,y2],width,/gal
;Revise history

if n_params() lt 4 then begin
    print, 'Syntax - PVBELT, fitsfile, a, d, width, [/gal, step=, pcoordiante= ]'
    return
endif
if ~file_test(fitsfile) then begin
        print,'Error: Fits file not exist!'
        return
end
fits_read,fitsfile,dat,hdr
;dat[where(dat eq max(dat))] = 0

if keyword_set(gal) then gal = 1b else gal = 0b
fitsgal = strcmp(sxpar(hdr,'CTYPE1'), 'GL', 2, /fold_case)
if fitsgal eq gal then at=a & dt=d
if fitsgal and ~gal then glactc,a,d,2000,at,dt,1,/deg
if ~fitsgal and gal then glactc,at,dt,2000,a,d,2,/deg
sxaddpar,hdr,'CTYPE1',repstr(sxpar(hdr,'CTYPE1'),'GLS','SFL')
sxaddpar,hdr,'CTYPE2',repstr(sxpar(hdr,'CTYPE2'),'GLS','SFL')
adxy,hdr,at,dt,x,y

defstep=0.5     ;unit pixel
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
xg = fltarr(nstep+1,width)
yg = fltarr(nstep+1,width)
width = fix(width)
temp = sqrt((xs[0]-xs[1])^2+(ys[0]-ys[1])^2)
xg[0,*] = xs[0]+(findgen(width)-(width-1)/2.)/temp*(ys[1]-ys[0])
yg[0,*] = ys[0]-(findgen(width)-(width-1)/2.)/temp*(xs[1]-xs[0])
for i=1, nstep do begin
temp = sqrt((xs[i]-xs[i-1])^2+(ys[i]-ys[i-1])^2)
xg[i,*] = xs[i]+(findgen(width)-(width-1)/2.)/temp*(ys[i]-ys[i-1])
yg[i,*] = ys[i]-(findgen(width)-(width-1)/2.)/temp*(xs[i]-xs[i-1])
endfor

;plot,x,y,/iso,/yno
;oplot,xg,yg,psym=3

nx1 = sxpar(hdr,'NAXIS3')
nx2 = n_elements(xs)
slice = make_array(nx1, nx2, type=size(dat,/type))
for i=0,nx1-1 do begin
	g = interpolate(dat[*,*,i],xg,yg,missing=0)
	slice[i,*] = total(g,2)/width
endfor
mkhdr,pvhdr,slice
factor = abs(sxpar(hdr,'CDELT3') lt 100)?(1d):(1000d)
sxaddpar,pvhdr,'CTYPE1','VELOCITY'
sxaddpar,pvhdr,'CRPIX1',sxpar(hdr,'CRPIX3')
sxaddpar,pvhdr,'CRVAL1',sxpar(hdr,'CRVAL3')/factor
sxaddpar,pvhdr,'CDELT1',sxpar(hdr,'CDELT3')/factor

if n_elements(pcoordinate) lt 1 then pcoordinate = 0 else pcoordinate = fix(abs(pcoordinate))
xyad,hdr,xs,ys,as,ds
da = mean((as-shift(as,1))[1:*])
dd = mean((ds-shift(ds,1))[1:*])
case pcoordinate of
	1:begin
		sxaddpar,pvhdr,'CTYPE2',sxpar(hdr,'CTYPE1')
		sxaddpar,pvhdr,'CRPIX2',1
		sxaddpar,pvhdr,'CRVAL2',as[0]
		sxaddpar,pvhdr,'CDELT2',da
	end
	2:begin
		sxaddpar,pvhdr,'CTYPE2',sxpar(hdr,'CTYPE1')
		sxaddpar,pvhdr,'CRPIX2',1
		sxaddpar,pvhdr,'CRVAL2',ds[0]
		sxaddpar,pvhdr,'CDELT2',dd
	end
	else:begin
		sxaddpar,pvhdr,'CTYPE2',sxpar(hdr,'POSITION')
		sxaddpar,pvhdr,'CRPIX2',1
		sxaddpar,pvhdr,'CRVAL2',0d
		sxaddpar,pvhdr,'CDELT2',step*abs(sxpar(hdr,'CDELT1'))
	end
endcase
sxaddhist,'PV file: '+fitsfile,pvhdr
sxaddhist,'PV path:',pvhdr
for i=0,n_elements(a)-1 do sxaddhist,string(a[i])+' '+string(d[i]),pvhdr
sxaddhist,'Position in Degree',pvhdr
sxaddhist,'Velocity in km/s',pvhdr
fits_write,'pvbelt.fits',slice,pvhdr
end
