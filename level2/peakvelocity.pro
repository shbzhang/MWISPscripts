pro peakvelocity, fitsfile, outputfile, vrange=vrange
;derive a peak velocity map
;each pixel represents the velocity of the peak on the spectrum
if n_params() lt 2 then begin
	print, 'Syntax - PEAKVELOCITY, fitsfile, output, vrange=[v1, v2]'
	return
endif
fits_read, fitsfile, dat, hdr
datdim = size(dat,/dimension)

if ~keyword_set(vrange) then vrange=[-300,300]
v2c, hdr, vrange*1000., crange
crange = round(crange >0 <(datdim[2]-1))
crange = crange[sort(crange)]
if crange[0] eq crange[1] then return
if crange[0] ge 0 then dat[*,*,0:crange[0]]=0
if crange[1] le datdim[2]-1 then dat[*,*,crange[1]:*]=0

peak = max(dat, sub, dimension=3, /nan)
xyc = array_indices(datdim, sub, /dimensions)
chan = reform(xyc[2,*], datdim[0], datdim[1])
c2v, hdr, chan, velo
fits_write, outputfile, velo/1000., hdr
end
