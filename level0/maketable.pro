pro maketable, x, y, t
;make a obs table
if n_params() lt 2 then begin
  print,'Syntax - maketable, gl, gb, [scanrate]'
  return
endif
;offset='81.492802 -2.9057733'
;offset='85.786559 -0.31681212'
;offset='86.611375 -0.65384447'
;offset='87.654871 -0.98327155'
;offset='87.980669 -1.786982'
;offset='88.232434 -0.73621614'
;offset='88.836972 -2.3347563'

offset='89.949759 5.5212777'

;offset='94.835056 -4.2888827'
;offset='92.169345 -4.8215622'
;offset='93.256539 -3.4517661'
;offset='94.087813 -2.5690013'
off=double(strsplit(offset,/extract))
;	d_el=(azel_pmo(x,y)).el-(azel_pmo(off[0],off[1])).el
;print,d_el

sou=string(x*10,format='(I04)')+string(y*10,format='(I+04)')
if y eq -0.5 then begin
	y=0
	oy='-1800'
endif else oy='0'
pos=string(x,format='(D9.5)')+"  "+string(y,format='(D+9.5)')
if n_params() le 2 then begin
        r='50'
	t='0.3'
endif else begin
	r=string(15./t,format='(i0)')
	t=string(t,format='(f3.1)')
endelse

openw,lun,sou+'.table',/get_lun
printf,lun,'1  '+sou+'     Galactic   '+pos+'    00.0   1  1  0  '+oy+'  '+offset
printf,lun,'REPEAT      2'
printf,lun,'THETAROW    15'
printf,lun,'WIDTH       1800'
printf,lun,'HEIGHT      1800'
printf,lun,'RAMP        10'
printf,lun,'SCANRATE    '+r
printf,lun,'TDUMP       '+t
printf,lun,'CALRATE     4'
printf,lun,'OFFINTTIME  2'
printf,lun,'CALINTTIME  2'
printf,lun,'ROTANG      360'
printf,lun,''
printf,lun,'OnceSweepTime(Min)   85'
printf,lun,'Regrid 30 arcsec RMS(K)   0.22'
close,lun
free_lun,lun
print,'maketable done!'
end
