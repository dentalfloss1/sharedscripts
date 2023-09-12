import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.coordinates import SkyCoord
import datetime
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 
parser = argparse.ArgumentParser(prog='slidingetawindow',description='Makes plots and outputs eta values for sources.')
parser.add_argument('csvfile',type=str,help="""Supply csv file of the exact format: dtype=[('date','<U128'),('band','i8'),('freq','f8'),('fpk','f8'),('fpkerr','f8'),('fint','f8'),('finterr','f8'), ('id',int),('det_sigma','f8'),('ra','f8'),('dec','f8'),('distance','f8'),('dataset','i8')])""")
args = parser.parse_args()
dat = np.loadtxt(args.csvfile, delimiter=',', skiprows=1, dtype=[('date','<U128'),('band','i8'),('freq','f8'),('fpk','f8'),('fpkerr','f8'),
    ('fint','f8'),('finterr','f8'), ('id',int),('det_sigma','f8'),('ra','f8'),('dec','f8'),('distance','f8'),('dataset','i8')])
specialtransients = np.array((688400,715724,719331,706803,696463,708691,721255,703942,689430,702639,703580,721648,703513,715880,702355,713717,722680,703195,694407,705714,703706,688586,714302,703241,694226,695801,686091,694120,713623,716372,694324,702999,703664,719717,707066,720097,721079,689811,689071,714807,696899,714671,705049,689546,695276,721006,702852,689539,687700,688132,687657,702679,707416,705887,691897,706436,702643,702324,707219,708405,696983,690365,688613,719978,714170,713985,690736,705197,705126,705482,704780,688164,690696,716711,701445,708685,703726,688222,703223,702092,694245,711246,712134,695623,705414,713705,692317,702565,707852,686547,698797,707094,721883,707144,690557,711209,709448,712910,697775,718738,693964,686640,717041,690829,707028,697401,701994,707287,704660,695871,694542,708781,707005,702209,704100,689549,702086,720496,709415,706011,714730,701632))
# moddat = np.loadtxt('out.csv',delimiter=',', skiprows=1, dtype=[('id','i8'),('ra','f8'),('dec','f8'),('mod','f8'),('moderr','f8')])


# https://stackoverflow.com/questions/18915378/rounding-to-significant-figures-in-numpy
def signif(x, p):
    x = np.asarray(x)
    x_positive = np.where(np.isfinite(x) & (x != 0), np.abs(x), 10**(p-1))
    mags = 10 ** (p - 1 - np.floor(np.log10(x_positive)))
    return np.round(x * mags) / mags


for myfreq in np.unique(dat['band']):
    srcs = np.unique(dat['id'][(dat['band']==myfreq)])
    datsummary = np.zeros(srcs.shape,dtype=[('id','i8'),('ra','f8'),('dec','f8'),('distance','f8')])
    for i, mysrc in enumerate(srcs):
        datsummary['id'][i] = mysrc
        datsummary['ra'][i] = np.average(dat['ra'][dat['id']==mysrc])
        datsummary['dec'][i] = np.average(dat['dec'][dat['id']==mysrc])
        datsummary['distance'][i] = np.average(dat['distance'][dat['id']==mysrc])
    sameoom = (signif(dat['fint'],1) == signif(dat['fpk'],1))
    keepdat = dat[(dat['distance'] < 0.8) & (dat['band']==myfreq) & (sameoom)]
    keepsrcs = np.unique(keepdat['id'][(keepdat['det_sigma'] > 6.4)])
    neweta = np.zeros(len(keepsrcs),dtype=[('id',int),('eta','f8'),('v','f8'),('ra','f8'),('dec','f8'),('distance','f8'),('dataset','i8'),('det_sigma','f8')])
    oldeta = np.zeros(len(keepsrcs),dtype=[('id',int),('eta','f8')])
    for i in range(len(keepsrcs)):
        srcdat = keepdat[keepdat['id']==keepsrcs[i]]
        neweta['ra'][i] = np.average(srcdat['ra'])
        neweta['dec'][i] = np.average(srcdat['dec'])
        neweta['distance'][i] = np.average(srcdat['distance'])
        neweta['dataset'][i] = int(np.unique(srcdat['dataset']))
        neweta['det_sigma'][i] = np.amax(srcdat['det_sigma'])
       #  if np.sum(moddat['id']==keepsrcs[i]) > 0:
       #      mymod = float(moddat['mod'][moddat['id']==keepsrcs[i]])
       #  else:
       #      mymod = 0
        # srcerr = np.sqrt(srcdat['finterr']**2 + (srcdat['fint']*0.1)**2 + (srcdat['fint']*mymod)**2)
        neweta['v'][i] = np.std(srcdat['fint'],ddof=1)/np.average(srcdat['fint'])
        preveta = 0
        datesort =np.argsort( [datetime.datetime.strptime(d, "%Y-%m-%d %H:%M:%S.%f") for d in srcdat['date']])
        
        if len(srcdat['fint']) > 1:
            if len(srcdat['fint']) <= 14:
                srcerr = np.sqrt(srcdat['finterr']**2 + (srcdat['fint']*0.1)**2 )
                xi = np.average(srcdat['fint'],weights=1/srcerr**2)
                neweta['eta'][i] = (1 / (len(srcdat['fint']) - 1)) * np.sum((srcdat['fint'] - xi)**2 / (srcerr**2))

            sorteddata = srcdat[datesort]
            for j in range(len(sorteddata) - 14):
                src = sorteddata[j:(j+14)]
                srcerr = np.sqrt(src['finterr']**2 + (src['fint']*0.1)**2 )
                xi = np.average(src['fint'],weights=1/srcerr**2)
                nexteta = (1 / (len(src['fint']) - 1)) * np.sum((src['fint'] - xi)**2 / (srcerr**2))
                neweta['eta'][i] = max(preveta,nexteta) 
                if src['id'][0]==123366:
                    print(preveta,nexteta)
                    input('presskey')
                preveta = neweta['eta'][i]
            neweta['id'][i] = keepsrcs[i]
        else:
            neweta['eta'][i] = np.inf
            neweta['id'][i] = keepsrcs[i]
        # NOW WE DO IT AGAIN WITHOUT SYSTEMATIC ERROR FOR OLD ETA
        srcerr = np.sqrt(srcdat['finterr']**2 )
        xi = np.average(srcdat['fint'],weights=1/srcerr**2)
        if len(srcdat['fint']) > 1:
            oldeta['eta'][i] = (1 / (len(srcdat['fint']) - 1)) * np.sum((srcdat['fint'] - xi)**2 / (srcerr**2))
            oldeta['id'][i] = keepsrcs[i]
        else:
            oldeta['eta'][i] = np.inf
            oldeta['id'][i] = keepsrcs[i]
    etadat = []
    for src in keepsrcs:
        etadat.append((src,*neweta['eta'][neweta['id']==src], *oldeta['eta'][oldeta['id']==src], *neweta['v'][neweta['id']==src]))  
    etaplot = np.array(etadat, dtype=[('id',int),('neweta','f8'),('oldeta','f8'),('v','f8')])
    etadatarr = np.array(etadat, dtype=[('id',int),('neweta','f8'),('oldeta','f8'),('v','f8')])
    etaplot = etadatarr[(etadatarr['neweta'] > 0) & (etadatarr['neweta'] < np.inf) & (etadatarr['oldeta'] > 0) & (etadatarr['oldeta'] < np.inf)]
    fig, axs=  plt.subplots(1,2, sharex=True, sharey=True, figsize=(16,8))
    for i in range(len(axs)):
        ax = axs[i]
        binmin = np.minimum(etaplot['neweta'].min(), etaplot['oldeta'].min())
        binmax = np.maximum(etaplot['neweta'].max(), etaplot['oldeta'].max())
        print(binmin,binmax)
        bins = np.geomspace(binmin,binmax, num=50)
        if i==0:
            ax.hist(etaplot['neweta'], bins=bins)
        elif i==1:
            ax.hist(etaplot['oldeta'], bins=bins)
        ax.set_xlabel('$\eta$')
        ax.set_ylabel('Number of transient sources')
        ax.set_xscale('log')
        ax.set_title(['$\eta$ with 10% systematic','$\eta$ before'][i])
        ax.axvline(1)
    plt.tight_layout()
    plt.savefig(f'f{myfreq}etahistint.png')
    plt.close()
    
    fig = plt.figure()
    plt.scatter(etaplot['neweta'], etaplot['v'],s=1)
    plt.scatter(etaplot['neweta'][np.isin(etaplot['id'],specialtransients)], np.unique(etaplot['v'][np.isin(etaplot['id'],specialtransients)]),marker='x', color='red')
    ax = plt.gca()
    ax.set_xlabel('new $\eta$')
    ax.set_ylabel('$V$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.axvline(1)
    plt.tight_layout()
    plt.savefig(f'f{myfreq}newetavint.png')
    plt.close()
    
    fig = plt.figure()
    candeta = etaplot['neweta'][etaplot['neweta'] > 1]
    bins = np.geomspace(candeta.min(),candeta.max(), num=25)
    plt.hist(candeta,bins=bins)
    ax = plt.gca()
    ax.set_xlabel('Corrected $\eta$')
    # ax.set_xscale('log')
    plt.tight_layout()
    plt.savefig(f'f{myfreq}newetahist.png')
    plt.close()
    
    for i in range(len(bins)-1):
        lbin = bins[i]
        rbin = bins[i+1]
        # print('-------η > ',lbin,'-------')
        myetas = etadatarr[(etadatarr['neweta'] >=lbin) & (etadatarr['neweta'] <=rbin)]
        for e in np.sort(myetas,order='neweta'):
            print(*e)
    
    # np.savetxt('mysortedsrcs.csv',etadatarr, fmt='%i,%f,%f,%f')
    np.savetxt(f'f{myfreq}mysortedsrcs.csv',neweta[['id','eta','v','ra','dec','distance']],fmt='%i,%f,%f,%f,%f,%f')
    np.savetxt(f'f{myfreq}varmetricdata.csv',neweta[['eta','v','dataset','id','det_sigma']],fmt='%f,%f,%i,%i,%f')
    fig, axs = plt.subplots(2,2, sharex=False, sharey=False, figsize=(25,25))
    binwidth=0.05
    hist, bins = np.histogram(neweta['distance'],bins=int(np.round((neweta['distance'].max() - neweta['distance'].min())/binwidth)))
    myareatotals = np.array([2*np.pi*(1-np.cos(theta*np.pi/180)) for theta in bins])
    myareanorms = np.array([(myareatotals[i+1] - myareatotals[i])*np.pi/180 for i in range(len(myareatotals)-1)])
    axs[0][0].bar((bins[:-1] + (bins[1:] - bins[:-1])/2), hist/myareanorms, bins[1:] - bins[:-1])
    axs[0][0].set_xlabel('degrees from center of field')
    axs[0][0].set_ylabel('Transients per square deg.')
    axs[0][0].set_title('Transients Histogram per square degree')
    axs[0][0].set_xlim([0,1.4])
    axs[1][0].bar((bins[:-1] + (bins[1:] - bins[:-1])/2), hist, bins[1:] - bins[:-1])
    axs[1][0].set_xlabel('degrees from center of field')
    axs[1][0].set_ylabel('Transients sources not normalized')
    axs[1][0].set_title('Transients Histogram')
    axs[1][0].set_xlim([0,1.4])
    fullhist, bins = np.histogram(datsummary['distance'],bins=int(np.round((datsummary['distance'].max() - datsummary['distance'].min())/binwidth)))
    myareatotals = np.array([2*np.pi*(1-np.cos(theta*np.pi/180)) for theta in bins])
    myareanorms = np.array([(myareatotals[i+1] - myareatotals[i])*np.pi/180 for i in range(len(myareatotals)-1)])
    axs[0][1].bar((bins[:-1] + (bins[1:] - bins[:-1])/2), fullhist/myareanorms, bins[1:] - bins[:-1])
    axs[0][1].set_xlabel('degrees from center of field')
    axs[0][1].set_ylabel('All sources per square degree')
    axs[0][1].set_title('All Sources Histogram Normalized by Area')
    axs[0][1].set_xlim([0,1.4])
    axs[1][1].bar((bins[:-1] + (bins[1:] - bins[:-1])/2), fullhist, bins[1:] - bins[:-1])
    axs[1][1].set_xlabel('degrees from center of field')
    axs[1][1].set_ylabel('All sources not normalized')
    axs[1][1].set_title('All Sources Histogram')
    axs[1][1].set_xlim([0,1.4])
    plt.savefig(f'f{myfreq}histdist.png')
    plt.close()
    
    # pointcenter = SkyCoord(186.693958,8.882694,unit='deg')
    # sc = SkyCoord(neweta['ra'],neweta['dec'],unit='deg')
    
    # mysrcs = neweta[(pointcenter.separation(sc).degree < 6)]
    
    # np.savetxt('SN2021smjqvv.csv',mysrcs, fmt='%i,%f,%f,%f,%f,%f')
    # myvarsrcs = mysrcs[mysrcs['eta'] > 2]
    # 
    # fig, axs = plt.subplots(2,2, sharex=False, sharey=False, figsize=(25,25))
    # binwidth=0.05
    # hist, bins = np.histogram(myvarsrcs['distance'],bins=int(np.round((myvarsrcs['distance'].max() - myvarsrcs['distance'].min())/binwidth)))
    # myareatotals = np.array([2*np.pi*(1-np.cos(theta*np.pi/180)) for theta in bins])
    # myareanorms = np.array([(myareatotals[i+1] - myareatotals[i])*np.pi/180 for i in range(len(myareatotals)-1)])
    # axs[0][0].bar((bins[:-1] + (bins[1:] - bins[:-1])/2), hist/myareanorms, bins[1:] - bins[:-1])
    # axs[0][0].set_xlabel('degrees from center of field')
    # axs[0][0].set_ylabel('Transients per square deg.')
    # axs[0][0].set_title('Transients Histogram per square degree')
    # axs[0][0].set_xlim([0,1.4])
    # axs[1][0].bar((bins[:-1] + (bins[1:] - bins[:-1])/2), hist, bins[1:] - bins[:-1])
    # axs[1][0].set_xlabel('degrees from center of field')
    # axs[1][0].set_ylabel('Transients sources not normalized')
    # axs[1][0].set_title('Transients Histogram')
    # axs[1][0].set_xlim([0,1.4])
    # fullhist, bins = np.histogram(mysrcs['distance'],bins=int(np.round((mysrcs['distance'].max() - mysrcs['distance'].min())/binwidth)))
    # myareatotals = np.array([2*np.pi*(1-np.cos(theta*np.pi/180)) for theta in bins])
    # myareanorms = np.array([(myareatotals[i+1] - myareatotals[i])*np.pi/180 for i in range(len(myareatotals)-1)])
    # axs[0][1].bar((bins[:-1] + (bins[1:] - bins[:-1])/2), fullhist/myareanorms, bins[1:] - bins[:-1])
    # axs[0][1].set_xlabel('degrees from center of field')
    # axs[0][1].set_ylabel('All sources per square degree')
    # axs[0][1].set_title('All Sources Histogram Normalized by Area')
    # axs[0][1].set_xlim([0,1.4])
    # axs[1][1].bar((bins[:-1] + (bins[1:] - bins[:-1])/2), fullhist, bins[1:] - bins[:-1])
    # axs[1][1].set_xlabel('degrees from center of field')
    # axs[1][1].set_ylabel('All sources not normalized')
    # axs[1][1].set_title('All Sources Histogram')
    # axs[1][1].set_xlim([0,1.4])
    # plt.savefig('histdist.png')
    
    # fig, axs = plt.subplots(1,2, sharex=False,sharey=False, figsize=(25,15))
    # bins = np.linspace( skyseps['fint'].min(), skyseps['fint'].max(), num =30)
    # axs[0].hist(skyseps['fint'], bins =bins)
    # axs[0].set_ylabel('Transients')
    # axs[0].set_xlabel('Fint (Jy)')
    # # axs[0].set_xscale('log')
    # axs[0].set_title('Histogram of Integrated Fluxes')
    # fullplotdat = fullseps['fint'][(fullseps['fint'] > 0) & (fullseps['fint'] < np.inf)]
    # # Note that we set the bins according to the flux range of the transients. The full
    # # range would be blown out. There is a single value that is extremely high.
    # bins = np.linspace(skyseps['fint'].min(), skyseps['fint'].max(), num=30)
    # axs[1].hist(fullplotdat, bins=bins)
    # axs[1].set_ylabel('Sources')
    # axs[1].set_xlabel('Fint (Jy)')
    # # axs[1].set_xscale('log')
    # axs[1].set_title('Histogram of Integrated Fluxes')
    # plt.tight_layout()
    # plt.savefig('histtransflux.png')
    # plt.close()
    # 
