#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 09:58:04 2022

@author: jason

processing of 1999 Trimble Geoexploror XT data

"""


from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


meta=pd.read_csv('/Users/jason/Dropbox/AWS/GCNET/ancillary/GC-Net_info.csv')
meta=meta.rename(columns={"name": "name_long"})
meta=meta.rename(columns={"name_abv": "name"})
print(meta.columns)
print(meta)

#%%

meta['lat_1999']=np.nan
meta['lon_1999']=np.nan
meta['elev_1999']=np.nan
meta['elev_std_1999']=np.nan
meta['date_1999']=np.nan

inpath='/Users/jason/Dropbox/AWS/GCNET/GCNet_info/GPS/'

dirs = sorted(glob(inpath+'*_gps.dat'))

# for idx in meta.id:

ly='x'

for f in dirs:
    file=f.split('/')[-1]
    idx=file[0:2]
    # print()
    df=pd.read_csv(f,names=['lat','lon','elev','datex','time'], header=None)
    df['date']=pd.to_datetime(df.datex)
    
    # if idx=='11':
    if idx!='null':
        plt.close()
        n_rows=6
        fig, ax = plt.subplots(n_rows,1,figsize=(16,20))
        cc=0
        v=np.where(meta.id==int(idx))
        v=v[0][0]
        ax[cc].set_title(meta.name[meta.id==int(idx)].values[0])
        ax[cc].hist(df.elev)
        ax[cc].set_xlabel('m')
        cc+=1
        ax[cc].plot(df.elev,'.')
        ax[cc].set_ylabel('m')
        cc+=1
        ax[cc].hist(df.lat)
        ax[cc].set_xlabel('deg. N')
        cc+=1
        ax[cc].plot(df.lat,'.')
        ax[cc].set_ylabel('deg. N')        
        cc+=1
        ax[cc].hist(df.lon)
        ax[cc].set_xlabel('deg. N')
        cc+=1
        ax[cc].plot(df.lon,'.')
        ax[cc].set_ylabel('deg. N')
        meta.loc[meta.id==int(idx),'elev_1999']=np.median(df.elev)
        meta.loc[meta.id==int(idx),'elev_1999']=meta.loc[meta.id==int(idx),'elev_1999'].map('{:.0f}'.format)
        meta.loc[meta.id==int(idx),'elev_std_1999']=np.std(df.elev)
        meta.loc[meta.id==int(idx),'elev_std_1999']=meta.loc[meta.id==int(idx),'elev_std_1999'].map('{:.0f}'.format)
        
        meta.loc[meta.id==int(idx),'lat_1999']=np.median(df.lat)
        meta.loc[meta.id==int(idx),'lon_1999']=np.median(df.lon)
        meta.loc[meta.id==int(idx),'date_1999']=df.date[0].strftime('%Y-%m-%d')
        print(idx,np.median(df.lat),np.median(df.lon),np.median(df.elev),np.std(df.elev),v)
        if ly=='p': plt.savefig(inpath+'/figs/'+meta.name[meta.id==int(idx)].values[0]+'.png', dpi=72,bbox_inches = 'tight')
        if ly=='x':plt.show()

#%%
# meta.to_csv('/Users/jason/Dropbox/AWS/GCNET/ancillary/GC-Net_info_incl_1999.csv')

meta.to_csv('/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/output/GC-Net_info_incl_1999.csv')


