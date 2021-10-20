One-line descriptions and examples of procedures

# level0
table preparation, data download/checking

>### IDL
>- OBSTABLE/OBSLIST: show suitable observing time for observing tables
>- MAKETABLE: use you own OFF point to write a observing table
>- DLBUR: download bur files from DLH database
>- DLFITS: download fits files from DLH database
>- MAPFILE: have quick look at the number and position of files in current directory
> <pre><code>IDL> mapfile, 80 </code></pre>
>### PYTHON
>- CELLMAP: have quick look at the number and position of files in certain directories
> <pre><code>>>> from cellmap import CellMap
> >>> CellMap('data/*.fits', l=[10,30], b=[-5,6]) </code></pre>


# level1
preliminary data reduction

>### GILDAS
>- AIO: MWISP pipeline - put all files in one file
> <pre><code>LAS> @aio 0810+010 U </code></pre>
>- BAS: MWISP pipeline - baseline fitting
> <pre><code>LAS> @base 0810+010U -200 200 -100 -60 -30 30 </code></pre>
>- SUM: MWISP pipeline - average spectra on each position
> <pre><code>LAS> @sum 0810+010U </code></pre>
>- MAKECUBE: MWISP pipeline - regrid to datacube
> <pre><code>LAS> @makecube 0810+010 U </code></pre>
>- MODIFY: MWISP pipeline - modify frequency from 13CO to C18O
> <pre><code>LAS> @modify 0810+010 </code></pre>
>- AREA: quick look of integrated intensity
> <pre><code>LAS> @area 0810+010U </code></pre>
>- RESAMPLE: resample the spectra and make a new datacube
> <pre><code>LAS> @resample 0810+010 U </code></pre>
>### IDL
>- GRID: generate grid file for regirding
> <pre><code>IDL> .r grid </code></pre>
>- GRID2: generate MWISP-II grid file for regirding
> <pre><code>IDL> .compile grid2
> IDL> grid_in_range, 30, 39.5, 5.5, 10 </code></pre>
>- CHECKOFFSET: search for possible OFF position from 12CO datacube
> <pre><code>IDL> checkoffset, '0810+030U.fits' </code></pre>

# level2
data reduction: datacube manipulation

>### IDL
>- CUBEREBIN: rebin the velocity dimension of a datacube
> <pre><code>IDL> cuberebin, 'test.fits' </code></pre>
>- CUBESMOOTH: convolve a datecube to a specified resolution
> <pre><code>IDL> cubesmooth, 'old.fits', 52/3600d, 60/3600d 'new.fits' </code></pre>
>- CUBECLIP: clip noise channels in a datecube
> <pre><code>IDL> cubeclip, 'old.fits', rmsfiler='old_rms.fits', coveragefile='old_cov.fits', outputfile='new.fits' </code></pre>
>- FITS_TRANSPOSE: transpose a datacube, e.g. from (l-b-v) to (l-v-b) or (b-v-l)
> <pre><code>IDL> fits_transpose, 'old.fits', 'new.fits', [3,1,2] </code></pre>
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
