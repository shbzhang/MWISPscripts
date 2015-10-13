pro rmshist, name, psfile=psfile, binsize=binsize, _extra=_extra
;plot the distribution of the rms
;Usage: rmshist, '0800+010', psfile='rmshist.ps', binsize=0.05
if n_params() lt 1 then begin
	print,'Syntax - RMSHIST, name, [psfile=, binsize=, charsize=, charthick= ]'
	return
endif
if keyword_set(psfile) then ps=1b else ps=0b
if ~keyword_set(binsize) then binsize=0.01
if ~keyword_set(charsize) then !p.charsize=1.5 else !p.charsize=charsize
if ~keyword_set(charthick) then !p.charthick=2 else !p.charthick=charthick
if ~keyword_set(thick) then !p.thick=2 else !p.thick=thick
if ~keyword_set(xthick) then !x.thick=2 else !x.thick=xthick
if ~keyword_set(ythick) then !y.thick=2 else !y.thick=ythick
fits_read,name+'_rms.fits',rms
fits_read,name+'_coverage.fits',cov
good=where(cov,/l64)
if good[0] ne -1 then rms=rms[good]
good=where(finite(rms),/l64)
if good[0] ne -1 then rms=rms[good]
med=median(rms)
hist=histogram(rms,binsize=binsize,location=x)
if ps then begin
	svd=!d
	svp=!p
	set_plot,'ps'
	device,filename=psfile,/color,xsize=20,ysize=16
endif
plot,x,hist,psym=10,xtitle='rms (K)', _extra=_extra
oplot,[med,med],[0,n_elements(rms)],linestyle=2,thick=2
xyouts,0.6,0.7,'!6Median = '+string(med,format='(f0.3)')+' K!X',/normal
if ps then begin
	device,/close
endif
end
