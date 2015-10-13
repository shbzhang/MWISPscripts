;cubemask
;by ShaoboZhang
;History
;Aug,20,2015,v1.0


function cubemask, datacube, threshold, nchannels
;Mask a datacube, only keep pixels with nchannels over threshold
;Input:
;	datacube: 3d data cube with x-y-v axes
;	threshold: a scalar of array indicates the threshold to be masked
;	nchannels: number of channels to be masked over threshold
;Usage: mask = cubemask(data, sigma*3, nchannels)
	if n_params() lt 2 then begin
	    print, 'Syntax - mask = CUBEMASK( datacube, threshold, nchannels )'
	    return, 0
	endif
	dim = size(datacube,/dimension)
	mask = bytarr(dim)
	for i=0,dim[2]-1 do mask[*,*,i] = datacube[*,*,i] gt threshold

        temp = mask
        for i=1,nchannels-1 do temp = temp and shift(mask,0,0,i)
        mask = temp
        for i=1,nchannels-1 do mask = mask or shift(temp,0,0,-i)
	return, mask
end
