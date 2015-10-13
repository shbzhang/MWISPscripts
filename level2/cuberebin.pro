pro cuberebin, filename
fits_read,filename+'.fits',d,h,/header_only
fits_read,filename+'_reb.fits',dr,hr
sxaddpar,h,'NAXIS3',sxpar(hr,'NAXIS3')
sxaddpar,h,'CRVAL3',sxpar(hr,'CRVAL3')
sxaddpar,h,'CDELT3',sxpar(hr,'CDELT3')
sxaddpar,h,'CRPIX3',sxpar(hr,'CRPIX3')
sxaddpar,h,'CROTA3',sxpar(hr,'CROTA3')
sxaddhist,'CUBEREBIN: '+systime(0),h
fits_write,'rebin.fits',dr,h
end
