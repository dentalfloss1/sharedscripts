# USED to create images of all the integration times in the dataset
# USAGE:
# in bash:
# for i in {1..5}; do casa -c autointim.py my.ms $i 5 & done


import glob
import datetime
# Used to create unique filenames
starttime = datetime.datetime.now().strftime('%y%m%d%Hh%Mm%Ss')
from casatools import msmetadata, ms 
import sys
import argparse
import os




def main_proc_loop(targetobs):
    """Takes as input a ms, outputs a modified list that will become the obs file"""
   
    def get_integrationtime(s):
        """Takes as input an index to start looking for the integration time. Recursively increases this index until 
        integration time is sucessfully retrieved from the ms without error and returns this time"""

        print(s) # just to let the user know something is happening
        msobj = ms()
        msobj.open(targetobs)
        try:
            return msobj.getscansummary()[str(s)]['0']['IntegrationTime']
        except KeyError:
            return get_integrationtime(s+1)

    try:
        # s here doesn't really matter much it just sets an index that will be incremented if the function fails, so start low
        integration_time = get_integrationtime(1)
        # The epoch commonly used in ms 
        start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)
        def writescansfile(beginepoch,epochduration,scans, scantime,scansfile):
            """Takes as input info about a long observation with gaps, writes a file detailing the scans, and returns the total time in the gaps"""

            gaps = []
            gaptime = 0
            scanlist = []
            for i in range(len(scantime)-1):
                startgap = max(scantime[i])+datetime.timedelta(seconds=round(integration_time)/2.0)
                endgap = min(scantime[i+1])+datetime.timedelta(seconds=round(integration_time)/2.0)
                startscan = min(scantime[i])+datetime.timedelta(seconds=round(integration_time)/2.0)
                endscan = max(scantime[i+1])+datetime.timedelta(seconds=round(integration_time)/2.0)
                gaps.append([startgap,endgap])
                gaptime = gaptime + (min(scantime[i+1]) - max(scantime[i])).total_seconds()
            for s,t in zip(scans,scantime):
                scanstart = (t[0] - datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                scanduration = ((t[-1] - t[0]) + datetime.timedelta(seconds=round(integration_time)/2.0)).total_seconds()
                scanend = (t[0]).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                sensitivity = constant/np.sqrt(scanduration)
                scanlist.append([scanstart, scanend,  rng.normal(sensitivity, 0.08*sensitivity), tmpra, tmpdec, False])
            with open(scansfile, "a+") as f:
                for t in scanlist:
                    f.write("{},{},{},{},{},{},{}\n".format(t[0], t[1], t[2], t[3], t[4], t[5], fov))
            print("wrote to "+scansfile)

            return gaptime
            
        msmd = msmetadata()
        msmd.open(msfile=targetobs)
        # Yes, this just assumes that the intent is TARGET. Maybe that can change with optional user input, but for now this seems fine
        for f in msmd.fieldsforintent('TARGET'):
            intnum = 0
            for tint in msmd.timesforfield(f):
                # YYYY/MM/DD/hh:mm:ss
                intstart = ((start_epoch + datetime.timedelta(seconds=tint)) - datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y/%m/%d/%H:%M:%S.%f')
                intduration = round(integration_time)
                intend = ((start_epoch + datetime.timedelta(seconds=intduration) + datetime.timedelta(seconds=tint)) - datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y/%m/%d/%H:%M:%S.%f')
                timelist.append([intstart+'~'+intend, intnum])
                intnum += 1 
            
    except RuntimeError:
        # Hopefully this doesn't happen, but if it does skipping the entire ms is probably the best option
        print("Issues with "+targetobs+". Skipping this one")

parser=argparse.ArgumentParser()
parser.add_argument("obs", help="ms to make intims of")
parser.add_argument("rank", help="for running multiple at the same time")
parser.add_argument("size", help="for running multiple at the same time")
# parser.add_argument("cutoff", help='cutoff flux for integration image source finding')
timelist= []
args = parser.parse_args()
size = args.size
rank = args.rank 
# cutoff = float(args.cutoff)
import numpy as np 
main_proc_loop(args.obs)  



numfile = 0
import datetime
import shutil
print("My rank is ",rank," of ",size)
# start = datetime.datetime.now()
# for t1, t2 in zip(timelist[(2*int(rank)-1)::2*int(size)],timelist[(1+2*int(rank)-1)::2*int(size)]):
    
#     # print(t)
#     # print(t[0])
#     # start = datetime.datetime.now()
#     tclean(vis=args.obs, timerange=t1[0], pblimit=-1e-12, cell='2arcsec',imsize=5120, imagename='im_alg1'+str(t1[1]), weighting='briggs', robust=0, niter=0, parallel=False)
#     tclean(vis=args.obs, timerange=t2[0], pblimit=-1e-12, cell='2arcsec',imsize=5120, imagename='im_alg1'+str(t2[1]), weighting='briggs', robust=0, niter=0, parallel=False)
#     i1 = 'im_alg1'+str(t1[1]) + '.image'
#     i2 = 'im_alg1'+str(t2[1]) + '.image'
#     immath(imagename=[i1,i2], mode="evalexpr", expr=' ( \"'+i2+'\" - \"'+i1+'\" )', imagemd=i1, outfile=i2.replace('.image','')+'m'+i1.replace('.image','')+'.image')
#     # exportfits('im'+str(t1[1])+'.image', fitsimage='im'+str(t[1])+'.fits')
#     # DO SOURCE FINDING HERE, KEEP FILES IF SOURCES IN FILES
#     for doomedfolder in glob.glob('im_alg1'+str(t1[1])+'.*'):
#         shutil.rmtree(doomedfolder, ignore_errors=True)
#     for doomedfolder in glob.glob('im_alg1'+str(t2[1])+'.*'):
#         shutil.rmtree(doomedfolder, ignore_errors=True)
#     numfile += 1 
#     if numfile == 3:
#         break 

# end = datetime.datetime.now()

# print('alg1: ',end-start)

start = datetime.datetime.now()

interesting_image = []
# for t1, t2 in zip(timelist[(2*int(rank)-1)::2*int(size)],timelist[(1+2*int(rank)-1)::2*int(size)]):
for t1 in timelist[int(rank)-1::int(size)]:
    
    # print(t)
    # print(t[0])
    # start = datetime.datetime.now()
    tclean(vis=args.obs, 
    timerange=t1[0],
    pblimit=-1e-12, 
    cell='2arcsec',
    imsize=5120, 
    imagename='im_'+str(t1[1]), 
    weighting='briggs', 
    robust=0, 
    niter=5000, 
    # gridder='wproject', 
    # wprojplanes=128,
    deconvolver='mtmfs',
    nterms=2,
    scales=[0,5,15],
    parallel=False)
    exportfits('im_'+str(t1[1])+'.image.tt0', fitsimage='im_'+str(t1[1])+'.fits')
    # tclean(vis=args.obs, timerange=t2[0], pblimit=-1e-12, cell='2arcsec',imsize=5120, imagename='im_dirty'+str(t2[1]), weighting='briggs', robust=0, niter=0, parallel=False)
    
#     i1 = 'im_dirty'+str(t1[1]) + '.image'
#     i2 = 'im_dirty'+str(t2[1]) + '.image'
    
#     r1 = imval(i1, box='0,0,5119,5119')
#     r2 = imval(i2, box='0,0,5119,5119')

    for doomedfolder in glob.glob('im_'+str(t1[1])+'.*'):
        if doomedfolder != 'im_'+str(t1[1])+'.fits':
            shutil.rmtree(doomedfolder, ignore_errors=True)
#     for doomedfolder in glob.glob('im_dirty'+str(t2[1])+'.*'):
#         shutil.rmtree(doomedfolder, ignore_errors=True)

#     rdiff = r2['data'] - r1['data']
    
#     print(np.sum(np.abs(rdiff) > cutoff))
#     if np.sum(np.abs(rdiff) > cutoff) > 0:
#         interesting_image.append(t1)
#         interesting_image.append(t2)


# end = datetime.datetime.now()


# print('Number of images that rank ',rank,' will make: ',len(interesting_image))

# for t in interesting_image:
#     tclean(vis=args.obs, 
#         timerange=t[0], 
#         pblimit=-1e-12, 
#         cell='2arcsec',
#         imsize=5120, 
#         imagename='im_cleaned'+str(t[1]),
#         gridder='wproject', 
#         wprojplanes=128,
#         deconvolver='mtmfs',
#         nterms=2,
#         scales=[0,5,15],
#         nsigma=3,
#         gain=0.8,
#         weighting='briggs', 
#         robust=0, 
#         niter=100000, 
#         parallel=False)


