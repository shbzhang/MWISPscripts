;Procedure for Milky Way Imaging Scroll Painting (MWISP) Project
;NAME:
;	grid.pro
;PURPOSE:
;	Generates the template grid file
;RESPONSIBLE PERSON
;	Shaobo Zhang
;CURRENT VERSION:
;	1.0
;REVISION HISTORY:
;	Written by Shaobo Zhang, 2011

;;;read cell ID
readcol,'mymap_v4.cat',num,id,gl,gb,format = 'i,a,d,d',skipline = 1
;id=['0940-055','0945-055']
;gl = [94,94.5]
;gb = [-5.5,-5.5]
;;;make a grid fits header
mkhdr,hdr,1,[1,1]
;;;position information
sxaddpar,hdr,'NAXIS1',91               ;dimension
sxaddpar,hdr,'NAXIS2',91
sxaddpar,hdr,'CTYPE1','GLON-GLS    '   ;projection type
sxaddpar,hdr,'CRVAL1',0                ;coordinate of reference pixel
sxaddpar,hdr,'CDELT1',-30d/3600.       ;pixelsize
sxaddpar,hdr,'CRPIX1',46               ;referencepixel
sxaddpar,hdr,'CROTA1',0
sxaddpar,hdr,'CTYPE2','GLAT-GLS    '
sxaddpar,hdr,'CRVAL2',0
sxaddpar,hdr,'CDELT2',30d/3600.
sxaddpar,hdr,'CRPIX2',46
sxaddpar,hdr,'CROTA2',0
;;;other keywords
sxaddpar,hdr,'EQUINOX',2000
sxaddpar,hdr,'OBJECT','DLHSURVEY'
;;;make the grid fits data
dat=bytarr(sxpar(hdr,'NAXIS1'),sxpar(hdr,'NAXIS2'))
;;;write grid
for i=0,n_elements(id)-1 do begin
sxaddpar,hdr,'CRVAL1',gl[i]
sxaddpar,hdr,'CRVAL2',gb[i]
sxaddpar,hdr,'CRPIX2',46
sxaddpar,hdr,'GLON',gl[i]
sxaddpar,hdr,'GLAT',gb[i]
;;;for region -1<gb<0
if (gb[i] gt -1) and (gb[i] lt 0) then begin
sxaddpar,hdr,'CRVAL2',0d
sxaddpar,hdr,'CRPIX2',round(-gb[i]*120+46)
endif
;;;write “LLLL±BBB_grid.fits”
fits_write,id[i]+'_grid.fits',dat,hdr
print,'Grid "'+id[i]+'.fits" has been generated successfully!'
endfor
end
