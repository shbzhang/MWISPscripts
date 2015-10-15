pro merge_bstrip, l, v1, v2, mergefile=mergefile, _extra=_extra
;merge bstrips
if n_params() lt 2 then begin
    print, 'Syntax - MERGE_BSTRIP, l, v1, v2, [mergefile=]'
    return
endif
if ~keyword_set(mergefile) then mergefile='mwisp_m0.fits'

name = 'BSTRIP_' + string(l,format='(I03)')
idx = where(file_test(name+'_U.fits'), count)
if count eq 0 then begin
	print, 'Error - can not find bstrip file'
	return
endif
name = name[idx]
for i=0, n_elements(name)-1 do begin
	cubemoment,name[i]+'_U.fits',[v1,v2],/zeroth_only,_extra=_extra
endfor
mwispmerge,name+'_U_m0.fits',name+'_U_coverage.fits',mergefile=mergefile
end

