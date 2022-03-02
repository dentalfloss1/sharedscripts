from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from scipy.special import binom
import numpy as np 
import matplotlib.pyplot as plt
import imageio
import os
import glob
import datetime

start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)
# imnames = ['im_'+str(i)+'.fits' for i in range(3,31,2)]
imnames = np.array(glob.glob('*fits'),dtype=str)
dates = np.zeros(imnames.shape, dtype=float)
for i in range(len(imnames)):
    f = imnames[i]
    hdu = fits.open(f)[0]
    dates[i] = (datetime.datetime.strptime(hdu.header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f") - start_epoch).total_seconds()/3600/24 
sortkey = np.argsort(dates)
for f in imnames[sortkey]:
    # filename = get_pkg_data_filename(f) #'/media/sarah/Elements/GRB200219A/1582287955_sdp_l0.GRB200219A_im_3.fits') #'E:\\GRB200219A\\1582287955_sdp_l0.GRB200219A_im_3.fits')
    hdu = fits.open(f)[0]
    wcs = WCS(hdu.header, naxis=2)
    ax = plt.subplot(projection=wcs)
    ax.imshow(hdu.data[0:][0:][0][0], vmin=-2.e-5, vmax=1.e-3, origin='lower')
    # ax.imshow(hdu.data, vmin=-2.e-5, vmax=5.e-4, origin='lower')
    plt.title(f)
    # ax.set_xlim(hdu.data.shape[3]/2 - 810 , hdu.data.shape[3]/2 + 550)
    # ax.set_ylim(hdu.data.shape[2]/2 - 700, hdu.data.shape[2]/2 + 550 )
    plt.savefig(f.rpartition('\\')[-1].replace('.fits','.png'))
with imageio.get_writer('animation.mp4', mode='I') as writer:
    for f in imnames[sortkey]:
        image = imageio.imread(f.rpartition('\\')[-1].replace('.fits','.png'))
        writer.append_data(image)
for f in imnames:
    os.remove(f.rpartition('\\')[-1].replace('.fits','.png'))
