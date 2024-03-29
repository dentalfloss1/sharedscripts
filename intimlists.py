# Used to run 8 second integration images through TraP in batches of adjacent images
# Needs updating for TraP v5.0rc2
# Needs updating for idia
# contact me for advice on how to do this

from astropy.io import fits
import glob
import numpy as np 
import datetime
from tqdm import tqdm
import subprocess
import shutil


start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00) # MJD epoch date
files = glob.glob('/idia/users/sarahchastain/images_to_measure/int/GRB220730A/*/*fits')
images = np.zeros(len(files), dtype={'names': ('imagename', 'date'), 'formats': ('U128','f8')})

for i in tqdm(range(len(files))):
    f = files[i]
    hdu = fits.open(f)
    header = hdu[0].header
    images['imagename'][i] = f
    images['date'][i] = (datetime.datetime.strptime(header['DATE-OBS'], "%Y-%m-%dT%H:%M:%S.%f") - start_epoch).total_seconds()/3600/24

images.sort(order='date')

listgod = [[images['imagename'][0]]]
listindex = 0
for i in range(len(images[1:])):
    if (images['date'][i] - images['date'][i-1]) < (8.01/3600/24):
        listgod[listindex].append(images['imagename'][i])
    else:
        listindex += 1 
        listgod.append([images['imagename'][i]])
for i in range(len(listgod)):
    subprocess.run(['singularity','exec','/idia/software/containers/trap/trapv5.0rc2.simg','trap-manage.py','initjob','chunk'+str(i)])
    shutil.rmtree('chunk'+str(i)+"/images_to_process.py", ignore_errors=True)
    shutil.rmtree('chunk'+str(i)+"/job_params.cfg", ignore_errors=True)
    imagesfilecontent = f"""
###################################################################################
#      List the images for processing by the transient detection pipeline         #
###################################################################################

# This should provide a module-scope iterable named "images" which provides a
# full path to each image to be processed. For example:

images = {listgod[i]}

# Optionally, whatever standard tools are required may be used to generate the
# list:
#
#  import os
#  import glob
#  images = sorted(
#      glob.glob(
#          os.path.expanduser("/home/example/data/*.fits")
#      )
#  )

#Display the list of images to be processed whenever this file is imported:
# (can be used for quick checking via an ipython import)
print "******** IMAGES: ********"
for f in images:
    print f
print "*************************" 
"""
    jobparamscontent = """[persistence]
description = "TRAP dataset"
dataset_id = -1
rms_est_sigma = 4            ; Sigma value used for iterative clipping in RMS estimation
rms_est_fraction = 8         ; Determines size of image subsection used for RMS estimation
bandwidth_max = 0.0          ; if non zero override bandwidth of image, determines which images fall in same band

[quality]
rms_est_history = 100        ; how many images used for calculating rms histogram
rms_est_max = 100            ; global maximum acceptable rms
rms_est_min = 0.0            ; global minimum acceptable rms
rms_rej_sigma = 3            ; threshold for rejecting images using rms histogram
oversampled_x = 30           ; threshold for oversampled check
elliptical_x = 10.0           ; threshold for elliptical check


[quality_lofar]              ; LOFAR only checks for casa images
low_bound = 1                ; multiplied with noise to define lower threshold
high_bound = 80              ; multiplied with noise to define upper threshold
min_separation = 10          ; minimum distance to a bright source (in degrees)

[source_extraction]
detection_threshold = 5      ; extraction threshold (S/N)
analysis_threshold = 3
back_size_x = 50
back_size_y = 50
margin = 10
deblend_nthresh = 0          ; Number of subthresholds for deblending; 0 disables
extraction_radius_pix = 2000
force_beam = True
box_in_beampix = 10
ew_sys_err = 10              ; Systematic errors on ra & decl (units in arcsec)
ns_sys_err = 10
expiration = 10              ; number of forced fits performed after a blind fit

[association]
deruiter_radius = 5.68
beamwidths_limit =  1.0

[transient_search]
new_source_sigma_margin = 1

[pipeline]
mode = 'batch'                            ; batch or stream
; below are the hosts and ports defined. Needs to be a string, if multiple
; hosts split by ,. Lengths need to match.
hosts = ',,,,,'                           ; if stream, the stream server
ports = '6666,6667,6668,6669,6670,6671'   ; the port of the stream
"""

    with open('chunk'+str(i)+'/images_to_process.py', 'w+') as f:
        f.write(imagesfilecontent)
    with open('chunk'+str(i)+'/job_params.cfg', 'w+') as f:
        f.write(jobparamscontent)

for i in range(len(listgod)):
    subprocess.run(['singularity','exec','/idia/software/containers/trap/trapv5.0rc2.simg','trap-manage.py','run','chunk'+str(i)], stderr=subprocess.STDOUT)
