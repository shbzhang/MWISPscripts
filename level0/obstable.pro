pro set_offset, gl, gb
;	gl = 81.492802d & gb = -2.9057733d
;	gl = 85.786559d & gb = -0.31681212d
;	gl = 86.611375d & gb = -0.65384447d
;	gl = 87.654871d & gb = -0.98327155d
;	gl = 87.980669d & gb = -1.786982d
;X	gl = 88.232434d & gb = -0.73621614d
;	gl = 88.836972d & gb = -2.3347563d
	gl = 89.949759d & gb = 5.5212777d

;	gl = 94.835056d & gb = -4.2888827d
;	gl = 92.169345d & gb = -4.8215622d
;	gl = 93.256539d & gb = -3.4517661d
;	gl = 94.087813d & gb = -2.5690013d
end

pro set_date, y, m, d
	y = 2015
	m = 09
	d = 05
end

pro set_range, l, b
	l = [89.,91] & b = [4,5]
	l = double(l[sort(l)]) & b = double(b[sort(b)])
end

function lct
	p=10
	return,dindgen(24l*p+1)/p
end

function cell_name, x, y
	return, string(round(x*10), format='(I04)') + string(round(y*10), format='(I+04)')
end

pro lb2hor,gl,gb,jd,az,el
glactc,ra,dec,2000,gl,gb,2,/degree
eq2hor,ra,dec,jd,el,az,lat=!pmo.latitude,lon=!pmo.longitude,altitude=!pmo.altitude
end

function el2at,el
;from elevation to atmosphere thick
Re=6371.393d;radius of earth in km
Ra=1000d;radius of atmosphere
theta = (el+90d)*!dtor
th_ce = !dpi-theta-asin(sin(theta)/(Re+Ra)*Re) ;theta in the earth center
return, (Re+Ra)/sin(theta)*sin(th_ce)
end



pro obstime, ls, bs, t
;give the best time for observation
if n_params() lt 2 then begin
	print,'Syntax - obstime, gl, gb, time'
	return
endif

ctime = lct()
set_date,y,m,d
jdcnv,y,m,d,ctime-8,jd      ;LCT -> UTC -> JULDAY

set_offset, lo, bo
lb2hor,lo,bo,jd,ao,eo

lb2hor,ls,bs,jd,as,es

d_e = es-eo
d_a = as-ao
d_at = (el2at(es) - el2at(eo))/el2at(es)

best = es gt 40 and abs(d_at) lt 0.001;2
better = es gt 30 and abs(d_at) lt 0.005;4
good = es gt 30 and abs(d_at) lt 0.01;6

if n_params() ge 3 then begin
t = good#[0,1,1,1,1,1,1,1,0]
t += better#[0,0,1,1,1,1,1,0,0]
t += best#[0,0,0,1,1,1,0,0,0]
endif else begin
ct2lst,lst,!pmo.longitude,8,ctime,d,m,y
idx = sort(lst)
window,xsize=800,ysize=700
top=0.95
plot,lst[idx],abs(d_a[idx]),position=[0.1,top/4*3.2,0.95,top],psym=0, $
        /xst,xtickinterval=1,xminor=1,/yst,yrange=[0,0],ytitle='!7D!X AZ'
plot,lst[idx],es[idx] >0,position=[0.1,top/4*2.2,0.95,top/4*3],/noerase,psym=0, $
        /xst,xtickinterval=1,xminor=1,/yst,yrange=[0,90],ytitle='Elevation'
plot,lst[idx],abs(d_e[idx]),position=[0.1,top/4*1.2,0.95,top/4*2],/noerase,psym=0, $
        /xst,xtickinterval=1,xminor=1,/yst,yrange=[0,3],ytitle='!7D!X EL'
plot,lst[idx],abs(d_at[idx]),position=[0.1,top/4*0.2,0.95,top/4*1],/noerase,psym=0, $
        /xst,xtickinterval=1,xminor=1,/yst,yrange=[0,20],xtitle='LST',ytitle='!7D!X AT'
oplot,lst[idx],good[idx]*11-1,psym=10,color='ff0000'x
oplot,lst[idx],better[idx]*9-1,psym=10,color='00ff00'x
oplot,lst[idx],best[idx]*7-1,psym=10,color='0000ff'x

endelse
end



pro obslist
;show all cells around offset, list the ideal time for observe each of them.
defsysv,'!pmo',{longitude:ten(97,33.6),latitude:ten(37,22.4),altitude:3200d}
set_offset, lo, bo
lo_cell = round(lo*2)/2.
bo_cell = round(bo*2)/2.
set_range,l_cell,b_cell

ctime = lct()
set_date,y,m,d
ct2lst,lst,!pmo.longitude,8,ctime,d,m,y
idx = sort(lst)

n_cell = (abs(l_cell[1]-l_cell[0])*2+1) * (abs(b_cell[1]-b_cell[0])*2+1)
;set plot background
window,xsize=800,ysize=700
plot, [0], xrange=[0,24], yrange=[0,n_cell], xstyle=5, ystyle=5
gtime = lst[idx] ge 19 or lst[idx] le 2
contour,gtime#[1,1], lst[idx], [0,n_cell], level=[1], c_colors=['ff5050'x], /fill, /overplot
contour,~gtime#[1,1], lst[idx], [0,n_cell], level=[1], c_colors=['a00000'x], /fill, /overplot
count=0
for x = l_cell[0], l_cell[1], 0.5 do begin
for y = b_cell[0], b_cell[1], 0.5 do begin
	name = cell_name(x,y)
	obstime, x, y, otime
	contour, otime[idx,*], lst[idx], count+findgen(9)/9, /fill, /overplot, $
		level=[0,1,2,3], c_color=['0'x,'5f'x,'8f00'x,'0000ff'x]
	xyouts,9,count,name
	count+=1
endfor
endfor
plot, [0], /noerase, xrange=[0,24], yrange=[0,n_cell], xstyle=1, ystyle=1, $
	xtitle='LST', ytitle='Cell', xtickinterval=3, xminor=3
end
