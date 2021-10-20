One-line descriptions of procedures

# level0
table preparation, data download/checking

>### IDL
>- OBSTABLE/OBSLIST: show suitable observing time for observing tables
>- MAKETABLE: use you own OFF point to write a observing table
>- DLBUR: download bur files from DLH database
>- DLFITS: download fits files from DLH database
>- MAPFILE: have quick look at the number and position of files in current directory
>### PYTHON
>- CELLMAP: have quick look at the number and position of files in certain directories


# level1
preliminary data reduction

>### GILDAS
>- AIO: MWISP pipeline - put all files in one file
>- BAS: MWISP pipeline - baseline fitting
>- SUM: MWISP pipeline - average spectra on each position
>- MODIFY: MWISP pipeline - modify frequency from 13CO to C18O
>- MAKECUBE: MWISP pipeline - regrid to datacube
>- AREA: quick look of integrated intensity
>- RESAMPLE: resample the spectra and make a new datacube
>### IDL
>- GRID: generate grid file for regirding
>- GRID2: generate MWISP-II grid file for regirding
>- CHECKOFFSET: search for possible OFF position from 12CO datacube


# level2
data reduction: datacube manipulation

>### IDL
>- CUBEREBIN: rebin the velocity dimension of a datacube
>- CUBESMOOTH: convolve a datecube to a specified resolution
>- CUBECLIP: clip noise channels in a datecube
>- FITS_TRANSPOSE: transpose a datacube, e.g. from (l-b-v) to (l-v-b) or (b-v-l)
>- HREBINV: rebin datacube for channel map
>- MOSAIC: mosaic cells of DLH survey
>- MOSAIC_CELL: mosaic a 30*30 arcmin cell
>- REPROJECT: reproject a datecube to a reference coordinate
>### PYTHON
>- MOSAIC: mosaic cells of DLH survey

data reduction: info extraction

>### IDL
>- AVERSPEC: derive the averaged spectrum of a datacube weighted with rms
>- CUBERMS: calculate rms noise for a datacube
>- CUBEMOMENT: calculate moment for a datacube
>- CUBEMASK: mask a datacube, only keep pixels with nchannels over threshold
>- PEAKVELOCITY: derive a peak velocity map of a datacube
>- PVSLICE: extract position-velocity map from a datacube
>- PVBELT: collapse position-velocity map within a belt of a datacube
>- RMSHIST: plot histogram of the rms distribution
>### PYTHON
>- CUBERMS: calculate rms noise for a datacube
>- CUBEMOMENT: calculate moment for a datacube
>- PVSLICE: extract position-velocity map from a datacube
>- REPROJECT_FITS: reproject an image to a reference coordinate

# level3
merge data to galactic plane

>### IDL
>- MOSAIC_BSTRIP: mosaic cells to a 1 degree width strip along galactic latitude
>- MERGE_BSTRIP: merge moment 0 of b-strips
>- MWISPMERGE: merge 2-dimensional data (e.g. moment map) to galactic plane map
>### PYTHON
>- TILE: tile images from separately mosaicked datacube.
