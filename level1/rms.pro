function getname, gl, gb
gl=gl mod 360
return,string(fix(gl*10),format='(I04)')+string(fix(gb*10),format='(I+04)')
end

pro rmslist, gl, gb
l = round(gl*2)/2d
b = round(gb*2)/2d

sb='U'
name=getname(l,b)
if ~file_test(name+sb+'_rms.fits') then begin
  print, 'Fits file "'+name+sb+'.fits" not exist!'
endif else begin
fits_read,name+sb+'_rms.fits',dat,hdr
dat=dat[15:75,15:75]
m=median(dat)
mc=median(dat[6:53,6:53])
v=variance(dat)
print,'rms of '+name+sb+': '+strjoin(string([m,v,mc],format='(F0.4)'),' ')+(m gt 0.55?'  ?':'')
endelse

sb='L'
name=getname(l,b)
if ~file_test(name+sb+'_rms.fits') then begin
  print, 'Fits file "'+name+sb+'.fits" not exist!'
endif else begin
fits_read,name+sb+'_rms.fits',dat,hdr
dat=dat[15:75,15:75]
m=median(dat)
mc=median(dat[6:53,6:53])
v=variance(dat)
print,'rms of '+name+sb+': '+strjoin(string([m,v,mc],format='(F0.4)'),' ')+(m gt 0.32?'  ?':'')
endelse
end

pro rms
for gl=104.5,100.5,-0.5 do begin
for gb=-2.0,0.5,0.5 do begin
rmslist,gl,gb
endfor
endfor
end
