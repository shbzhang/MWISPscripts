pro v2c,hdr,v,c
        if n_params() lt 3 then begin
                print, 'Syntax - V2C, hdr, v, c'
                return
        endif
        nc = sxpar(hdr,'NAXIS3')
        v0 = sxpar(hdr,'CRVAL3')
        c0 = sxpar(hdr,'CRPIX3')
        dv = sxpar(hdr,'CDELT3')
        c = (v-v0)/dv+c0-1
end
