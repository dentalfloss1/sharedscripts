from astropy.io import fits
import glob
file = glob.glob('*_IClean.fits')
for f in file:
    hdu = fits.open(f)[0]
    hdu.header['CTYPE3']  = 'FREQ'
    hdu.header['RESTFRQ'] = hdu.header['CRVAL3']
    fits.writeto(f.replace('_IClean', '_firstim'), data = hdu.data[0,0], header=hdu.header, overwrite=True)
