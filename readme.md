One-line descriptions and examples of procedures

# level0
table preparation, data download/checking

>### IDL
>- OBSTABLE/OBSLIST: show suitable observing time for observing tables
>- MAKETABLE: use you own OFF point to write a observing table
>- DLBUR: download bur files from DLH database
> <pre><code>IDL> dlbur, range=[80, 85, -1, 2], /append </code></pre>
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
>- AIO: [MWISP pipeline] put all files in one file
> <pre><code>LAS> @aio 0810+010 U </code></pre>
>- BAS: [MWISP pipeline] baseline fitting
> <pre><code>LAS> @base 0810+010U -200 200 -100 -60 -30 30 </code></pre>
>- SUM: [MWISP pipeline] average spectra on each position
> <pre><code>LAS> @sum 0810+010U </code></pre>
>- MAKECUBE: [MWISP pipeline] regrid to datacube
> <pre><code>LAS> @makecube 0810+010 U </code></pre>
>- MODIFY: [MWISP pipeline] modify frequency from 13CO to C18O
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
> <pre><code>IDL> fits_transpose, 'old.fits', 'new.fits', [0,2,1] </code></pre>
>- HREBINV: rebin datacube for channel map
> <pre><code>IDL> hrebinv, oldcub, oldhdr, newcub, newhdr, -200, 200, 401 </code></pre>
>- MOSAIC: mosaic cells of DLH survey
> <pre><code>IDL> mosaic, 80, 100, -3, 3, -200, 200, 'U', path=['./','data/'], /display </code></pre>
>- REPROJECT: reproject a datecube to a reference coordinate
> <pre><code>IDL> reproject, 'old.fits', 'ref.fits', 'new.fits' </code></pre>
>### PYTHON
>- MOSAIC: mosaic cells of DLH survey
> <pre><code>>>> from mosaic import mosaic
> >>> mosaic(80, 100, -3, 3, -200, 200, sb='U', path=['./','data/']) </code></pre>

data reduction: info extraction

>### IDL
>- AVERSPEC: derive the averaged spectrum of a datacube weighted with rms
> <pre><code>IDL> spec = averspec('test.fits', rmsfile = 'test_rms.fits') </code></pre>
>- CUBERMS: calculate rms noise for a datacube
> <pre><code>IDL> cuberms, 'cube.fits', [-100,100], window=[-50,-30,-10,10] </code></pre>
>- CUBEMOMENT: calculate moment for a datacube
> <pre><code>IDL> cubemoment, 'cube.fits', [-10,10], /goodlooking, rmsfile='rms.fits', coveragefile='coverage.fits' </code></pre>
>- CUBEMASK: mask a datacube, only keep pixels with nchannels over threshold
> <pre><code>IDL> cubemask, 'cube.fits', 1.5, 5 </code></pre>
>- PEAKVELOCITY: derive a peak velocity map of a datacube
> <pre><code>IDL> peakvelocity, 'cube.fits', 'peakvelocity.fits', [-30,30] </code></pre>
>- PVSLICE: extract position-velocity map from a datacube
> <pre><code>IDL> pvslice, 'cube.fits', [80,81,82], [-1,1,2], width=5, /spline </code></pre>
>- PVBELT: collapse position-velocity map within a belt of a datacube
>- RMSHIST: plot histogram of the rms distribution
> <pre><code>IDL> rmshist, '0800+010', psfile='rmshist.ps', binsize=0.05, xrange=[0,1] </code></pre>
>### PYTHON
>- CUBERMS: calculate rms noise for a datacube
> <pre><code>>>> from cuberms import cuberms
> >>> cuberms('test.fits', [-200, 200], window=[-100, -60, -30, 30]) </code></pre>
>- CUBEMOMENT: calculate moment for a datacube
> <pre><code>>>> from cubemoment import cubemoment
> >>> cubemoment('cube.fits', [-10,10], goodlooking=True, rmsfile='rms.fits', coveragefile='coverage.fits') </code></pre>
>- PVSLICE: extract position-velocity map from a datacube
> <pre><code>>>> from pvslice import pvslice
> >>> pvslice('test.fits', [[80,81,82], [-1,1,2]], width=5, step=0.2) </code></pre>
>- REPROJECT_FITS: reproject an image to a reference coordinate

# level3
merge data to galactic plane

>### IDL
>- MOSAIC_BSTRIP: mosaic cells to a 1 degree width strip along galactic latitude
>- MERGE_BSTRIP: merge moment 0 of b-strips
>- MWISPMERGE: merge 2-dimensional data (e.g. moment map) to galactic plane map
> <pre><code>IDL> mwispmerge, ['a_m0.fits','b_m0.fits'], ['a_coverage.fits','b_coverage.fits'], mergefile='mwisp_m0.fits' </code></pre>
>### PYTHON
>- TILE: tile images from separately mosaicked datacube.
> <pre><code>>>> from tile import tile
> >>> from glob import glob
> >>> tile(glob('*_U_lvmap.fits'), output='tile.fits') </code></pre>
