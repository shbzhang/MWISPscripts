pro mosaic_bstrip, l, v1, v2, fitspath=fitspath, rmspath=rmspath
;l=83
fitspath='/home/shbzhang/Lcloud/check/checked/'
rmspath='/home/shbzhang/Lcloud/check/checked/'
v1=-100 & v2=100

for i=0,n_elements(l)-1 do begin
outname = 'BSTRIP' + string(l[i],format='(I03)')
print,outname

mosaic,l[i],l[i]+1,-6,6,v1,v2,sb='U',fitspath=fitspath, rmspath=rmspath, display=0
if file_test('mosaic_U.fits') then begin
	file_move,'mosaic_U.fits',outname+'_U.fits',/overwrite
	file_move,'mosaic_U_rms.fits',outname+'_U_rms.fits',/overwrite
	file_move,'mosaic_U_coverage.fits',outname+'_U_coverage.fits',/overwrite
endif

mosaic,l[i],l[i]+1,-6,6,v1,v2,sb='L',fitspath=fitspath, rmspath=rmspath, display=0
if file_test('mosaic_L.fits') then begin
	file_move,'mosaic_L.fits',outname+'_L.fits',/overwrite
	file_move,'mosaic_L_rms.fits',outname+'_L_rms.fits',/overwrite
	file_move,'mosaic_L_coverage.fits',outname+'_L_coverage.fits',/overwrite
endif

mosaic,l[i],l[i]+1,-6,6,v1,v2,sb='L2',fitspath=fitspath, rmspath=rmspath, display=0
if file_test('mosaic_L2.fits') then begin
	file_move,'mosaic_L2.fits',outname+'_L2.fits',/overwrite
	file_move,'mosaic_L2_rms.fits',outname+'_L2_rms.fits',/overwrite
	file_move,'mosaic_L2_coverage.fits',outname+'_L2_coverage.fits',/overwrite
endif

endfor
end


pro merge_bstrip, l, irange
mwispmerge,mergefile='mwips_m0.fits'
for i=0, n_elements(l)-1 do begin
outname = 'BSTRIP' + string(l[i],format='(I03)')
if ~file_test(outname+'_U_.fits') then continue
cubemoment,outname+'_U.fits',irange,/zeroth_only
mwispmerge,outname+'_U_m0',outname+'_U_coverage.fits',mergefile='mwips_m0.fits'
endfor
end
