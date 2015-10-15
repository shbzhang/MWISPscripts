pro mosaic_bstrip, l, v1, v2, path=path
;mosaic strips along b direction
if n_params() lt 3 then begin
    print, 'Syntax - MOSAIC_BSTRIP, l, v1, v2, [path=]'
    return
endif

for i=0,n_elements(l)-1 do begin
outname = 'BSTRIP_' + string(l[i],format='(I03)')
print,outname:

mosaic,l[i],l[i]+1,-6,6,v1,v2,sb='U',path=path, display=0
if file_test('mosaic_U.fits') then begin
	file_move,'mosaic_U.fits',outname+'_U.fits',/overwrite
	file_move,'mosaic_U_rms.fits',outname+'_U_rms.fits',/overwrite
	file_move,'mosaic_U_coverage.fits',outname+'_U_coverage.fits',/overwrite
endif

mosaic,l[i],l[i]+1,-6,6,v1,v2,sb='L',path=path, display=0
if file_test('mosaic_L.fits') then begin
	file_move,'mosaic_L.fits',outname+'_L.fits',/overwrite
	file_move,'mosaic_L_rms.fits',outname+'_L_rms.fits',/overwrite
	file_move,'mosaic_L_coverage.fits',outname+'_L_coverage.fits',/overwrite
endif

mosaic,l[i],l[i]+1,-6,6,v1,v2,sb='L2',path=path, display=0
if file_test('mosaic_L2.fits') then begin
	file_move,'mosaic_L2.fits',outname+'_L2.fits',/overwrite
	file_move,'mosaic_L2_rms.fits',outname+'_L2_rms.fits',/overwrite
	file_move,'mosaic_L2_coverage.fits',outname+'_L2_coverage.fits',/overwrite
endif

endfor
end

