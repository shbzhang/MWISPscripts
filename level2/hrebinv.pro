;hrebinv
;by ShaoboZhang
;History:
;Jul,23,2015,v1.0
;Oct,26,2015,v1.1
;  fix error while the 'cdelt3' is a minus value

pro hrebinv, oldcub, oldhdr, newcub, newhdr, velo1, velo2, nchannel
;Rebin velocity axis for channel map
;the given velo1 and velo2 is the CENTER velocity of the first and last channel
;Input:
;	oldcub, oldhdr: the data and header of old fits
;	velo1, velo2: scalars indicate the velocity range
;	nchannel: number of channel of result cube
;Output:
;	newcub, newhdr: the data and header for new fits
;Usage: hrebinv, oldcub, oldhdr, newcub, newhdr, v1, v2, nc
	if n_params() lt 2 then begin
	    print, 'Syntax - HREBINV, oldcub, oldhdr, newcub, newhdr, velocity1, velocity2, num_channels'
	    print, 'Provided velocity should have the same unit as that in oldhdr.'
	    return
	endif
	v2c,oldhdr,velo1,c1
	v2c,oldhdr,velo2,c2
	width = float(abs(c1-c2))/(nchannel-1)
	c = min([c1,c2])-width/2 + findgen(nchannel+1)*width
	olddim = size(oldcub,/dimension)
	newcub = fltarr(olddim[0],olddim[1],nchannel)

	for i=0,nchannel-1 do begin
		if c[i] lt -0.5 or c[i+1] gt olddim[2]-0.5 then begin
			newcub[*,*,i] = !values.f_nan
			print, 'Warning - given velocity range exceeds the old datacube'
			continue
		end
		intc1 = ceil(c[i]+0.5)
		fltc1 = intc1-c[i]-0.5
		intc2 = floor(c[i+1]-0.5)
		fltc2 = c[i+1]-intc2-0.5

		if fltc1 eq 0 then head = 0 else head = (oldcub[*,*,intc1-1] * fltc1)
		if fltc2 eq 0 then rear = 0 else rear = (oldcub[*,*,intc2+1] * fltc2)
		case intc1-intc2 of
			0:begin
				newcub[*,*,i] = head + oldcub[*,*,intc1] + rear
			end
			1:begin
				newcub[*,*,i] = head + rear
			end
			2:begin
				newcub[*,*,i] = head /fltc1 *(fltc1+fltc2-1)
			end
			else:begin
				newcub[*,*,i] = head + total(oldcub[*,*,intc1:intc2],3) + rear
			end
		endcase
	endfor
	dvold = sxpar(oldhdr,'CDELT3')
	newcub *= abs(dvold)	;convert sum(Tmb) to integrated intensity
	newhdr = oldhdr
	c2v,oldhdr,min([c1,c2]),v0
	sxaddpar, newhdr, 'NAXIS3', nchannel
	sxaddpar, newhdr, 'CRPIX3', 1
	sxaddpar, newhdr, 'CRVAL3', v0
	sxaddpar, newhdr, 'CDELT3', width*dvold
end
