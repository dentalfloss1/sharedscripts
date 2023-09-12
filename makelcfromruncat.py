import psycopg2
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
import datetime
import argparse
from astropy.coordinates import SkyCoord
parser = argparse.ArgumentParser()
parser.add_argument('csvfile',type=str,help='csvfile from slidingetaint.py')
parser.add_argument('--usefpeak',action='store_true',help='Use peak flux instead of integrated flux')
args = parser.parse_args()



def initdb():
    import psycopg2
    import os
    from psycopg2.extensions import register_adapter, AsIs
    # psycopg2 complains about numpy datatypes, this avoids that error. Got this from stackoverflow but don't recall where
    def addapt_numpy_float64(numpy_float64):
        return AsIs(numpy_float64)
    def addapt_numpy_int64(numpy_int64):
        return AsIs(numpy_int64)
    register_adapter(np.float64, addapt_numpy_float64)
    register_adapter(np.int64, addapt_numpy_int64)
    host = os.getenv('TKP_DBHOST')
    port = os.getenv('TKP_DBPORT')
    user = os.getenv('TKP_DBUSER')
    password = os.getenv('TKP_DBPASSWORD')
    database = os.getenv('TKP_DBNAME')
    if password is None:
        conn = psycopg2.connect("dbname="+database+" user="+user+" host="+host)
    else:
        conn = psycopg2.connect("dbname="+database+" user="+user+" password="+password+" host="+host)

    ################################################################################
    # Typical psycopg2 setup 
    conn = psycopg2.connect("dbname="+database+" user="+user+" password="+password+" host="+host)
    return conn.cursor()


def makelc(dates,flux,fluxerr,runcat):
    from matplotlib import rcParams, rcParamsDefault
    rcParams.update(rcParamsDefault)
    fig = plt.figure()
    plt.scatter(dates,flux)
    plt.errorbar(dates,flux,yerr=np.sqrt(fluxerr**2+(0.1*flux)**2),fmt='none')
    ax=plt.gca()
    # ax.set_xscale('log')
    ax.set_xlabel('Days post trigger')
    ax.set_ylabel('$F_{int}$ (Jy)')
    ax.set_title(f'src {runcat}')
    plt.savefig(f'src{runcat}lc.png')
    plt.close()
    
def make2lc(obsnums,myruncat):
    font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 22}

    matplotlib.rc('font', **font)
    fig,axs = plt.subplots(1,2,figsize=(30,15),sharey=True)
    fig.suptitle(f'src {myruncat}')
    dummyindex = 0 
    for a in axs:
        obs = f'%{obsnums[dummyindex]}%'
        query = """SELECT image.taustart_ts, extractedsource.f_int, extractedsource.f_int_err FROM extractedsource JOIN 
        image ON image.id=extractedsource.image JOIN assocxtrsource ON 
        assocxtrsource.xtrsrc=extractedsource.id WHERE assocxtrsource.runcat=%s AND image.url LIKE %s ;""",
        if args.usefpeak:
            query = query.replace('f_int','f_peak')
        cur.execute(query,
        (myruncat,obs,))
        fetchedfluxerr = cur.fetchall()
        datlist = [ ((d-triggerdate).total_seconds()/3600./24.,f,fe) for d,f,fe in fetchedfluxerr]
        dat = np.array(datlist)
        # print(dat)
        a.scatter(dat[:,0],dat[:,1])
        a.errorbar(dat[:,0],dat[:,1],yerr=np.sqrt(dat[:,2]**2+(dat[:,1]*0.1)**2),fmt='none')
        # a.set_xscale('log')
        a.grid()
        if dummyindex==1:
            a.set_xlabel('Days post trigger')

        a.set_ylabel('$F_{int}$ (Jy)')
        dummyindex +=1 
    plt.tight_layout()
    plt.savefig(f'src{myruncat}lc_twopanel.png')
    plt.close()
    
def make3lc(obsnums,myruncat):
    font = {'family' : 'normal',
    'weight' : 'bold',
    'size'   : 22}

    matplotlib.rc('font', **font)
    fig,axs = plt.subplots(1,3,figsize=(40,15),sharey=True)
    fig.suptitle(f'src {myruncat}')
    dummyindex = 0 
    for a in axs:
        obs = f'%{obsnums[dummyindex]}%'
        query = """SELECT image.taustart_ts, extractedsource.f_int, extractedsource.f_int_err FROM extractedsource JOIN 
        image ON image.id=extractedsource.image JOIN assocxtrsource ON 
        assocxtrsource.xtrsrc=extractedsource.id WHERE assocxtrsource.runcat=%s AND image.url LIKE %s ;""",
        if args.usefpeak:
            query = query.replace('f_int','f_peak')
        cur.execute(query,
        (myruncat,obs,))
        fetchedfluxerr = cur.fetchall()
        datlist = [ ((d-triggerdate).total_seconds()/3600./24.,f,fe) for d,f,fe in fetchedfluxerr]
        dat = np.array(datlist)
        # print(dat)
        a.scatter(dat[:,0],dat[:,1])
        a.errorbar(dat[:,0],dat[:,1],yerr=np.sqrt(dat[:,2]**2+(dat[:,1]*0.1)**2),fmt='none')
        # a.set_xscale('log')
        a.grid()
        if dummyindex==2:
            a.set_xlabel('Days post trigger')

        a.set_ylabel('$F_{int}$ (Jy)')
        dummyindex +=1 
    plt.tight_layout()
    plt.savefig(f'src{myruncat}lc_threepanel.png')
    plt.close()
def make4lc(obsnums,myruncat):
    font = {'family' : 'normal',
    'weight' : 'bold',
    'size'   : 22}

    matplotlib.rc('font', **font)
    title= f'src{myruncat}lc_threepanel.png'
    fig,axs = plt.subplots(2,2,figsize=(30,20),sharey=True)
    fig.suptitle(f'src {title}')
    dummyindex = 0 
    for ax in axs:
        for a in ax:
            obs = f'%{obsnums[dummyindex]}%'
            # print(obs)
            # print(myruncat)
            query = """SELECT image.taustart_ts, extractedsource.f_int, extractedsource.f_int_err FROM extractedsource JOIN 
            image ON image.id=extractedsource.image JOIN assocxtrsource ON 
            assocxtrsource.xtrsrc=extractedsource.id WHERE assocxtrsource.runcat=%s AND image.url LIKE %s ;""",
            if args.usefpeak:
                query = query.replace('f_int','f_peak')
            cur.execute(query,
            (myruncat,obs,))
            fetchedfluxerr = cur.fetchall()
            # print(fetchedfluxerr)
            datlist = [ ((d-triggerdate).total_seconds()/3600./24.,f,fe) for d,f,fe in fetchedfluxerr]
            dat = np.array(datlist)
            # print(dat)
            a.scatter(dat[:,0],dat[:,1])
            a.errorbar(dat[:,0],dat[:,1],yerr=np.sqrt(dat[:,2]**2+(dat[:,1]*0.1)**2),fmt='none')
            # a.set_xscale('log')
            a.grid()
            if dummyindex in [2,3]:
                a.set_xlabel('Days post trigger')
            if dummyindex in [0,2]:
                a.set_ylabel('$F_{int}$ (Jy)')
            dummyindex +=1 
    plt.tight_layout()
    plt.show()
    plt.close()
    
if __name__=='__main__':

    infile = np.loadtxt(args.csvfile,delimiter=',',usecols=(0,1), dtype=[('id','i8'),('eta','f8')])
    srcstoplot = infile['id'][infile['eta'].round(decimals=2) >=2]
    cur = initdb()
    
    query = """SELECT runningcatalog.wm_ra, runningcatalog.wm_decl, runningcatalog.id FROM runningcatalog WHERE runningcatalog.id in %s;"""
    
    cur.execute(query,(tuple(srcstoplot),))
    
    srcinfo = np.array(cur.fetchall(),dtype=[('ra','f8'),('dec','f8'),('id','i8')])   

    
    triggers = np.loadtxt('commensal2triggerdates.csv',skiprows=1,delimiter=',', dtype=[('name','<U32'),('target','<U32'),('trigger','<U128'),('ra','f8'),('dec','f8')])
    triggerarr = np.array([(n,tar,datetime.datetime.strptime(trig,"%Y-%m-%dT%H:%M:%S.%f"),tra,tdec) for n,tar,trig,tra,tdec in triggers],dtype=[('name','<U32'),('target','<U32'),('trigger','O'),('ra','f8'),('dec','f8')])
    
    
    for src in srcinfo:
        sc1 = SkyCoord(src['ra'],src['dec'],unit='deg')
        sc2 = SkyCoord(triggerarr['ra'],triggerarr['dec'],unit='deg')
        seps = sc1.separation(sc2).arcsecond
        nobs = len(triggerarr[seps==seps.min()])
        triggerdate = np.unique(triggerarr['trigger'][seps==seps.min()])[0]
        query = """SELECT image.taustart_ts, extractedsource.f_int, extractedsource.f_int_err FROM extractedsource JOIN 
                image ON image.id=extractedsource.image JOIN assocxtrsource ON 
                assocxtrsource.xtrsrc=extractedsource.id WHERE assocxtrsource.runcat=%s ;"""
        if args.usefpeak:
            query = query.replace('f_int','f_peak')
        cur.execute(query, (src['id'],))
        fetchedfluxerr = cur.fetchall()
        datlist = [ ((d-triggerdate).total_seconds()/3600./24.,f,fe) for d,f,fe in fetchedfluxerr]
        dat = np.array(datlist)
        if nobs==1:
            makelc(dat[:,0],dat[:,1],dat[:,2],src['id'])
        else:
            try:
                if nobs==2:
                    make2lc(list(triggerarr['name'][seps==seps.min()]),src)
                elif nobs==3:
                    make3lc(list(triggerarr['name'][seps==seps.min()]),src)
                elif nobs==4:
                    make4lc(list(triggerarr['name'][seps==seps.min()]),src)                    
            except Exception as e:
                print(e)
                makelc(dat[:,0],dat[:,1],dat[:,2],src['id'])
             
                     
            
        
        
        
        
    
