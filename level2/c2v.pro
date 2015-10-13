pro c2v,hdr,c,v
        if n_params() lt 3 then begin
                print, 'Syntax - C2V, hdr, c, v'
                return
        endif
        nc = sxpar(hdr,'NAXIS3')
        v0 = sxpar(hdr,'CRVAL3')
        c0 = sxpar(hdr,'CRPIX3')
        dv = sxpar(hdr,'CDELT3')
        v = (c-c0+1)*dv+v0
end

