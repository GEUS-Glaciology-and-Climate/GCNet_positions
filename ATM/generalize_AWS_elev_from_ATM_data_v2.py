#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 08:51:17 2024

@author: jason
"""


from glob import glob
import pandas as pd
import numpy as np
import os

# ## change to your system's login name to change dir for local work
if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/AWS/GCNET/GCNet_positions.stash/'

os.chdir(base_path)


# read Greenland AWS locations
meta = pd.read_csv('/Users/jason/Dropbox/AWS/GCNET/GC-Net-level-1-data-processing/L1/GC-Net_location.csv')
meta = meta.rename({'Name': 'name'}, axis=1)
meta = meta.rename({'Latitude (°N)': 'lat'}, axis=1)
meta = meta.rename({'Longitude (°E)': 'lon'}, axis=1)
meta = meta.rename({'Elevation (wgs84 m)': 'elev'}, axis=1)
meta = meta.rename({'Crawford Point 1': 'Crawford Pt. 1'}, axis=1)
meta=meta.replace({'Crawford Point 1': 'Crawford Pt. 1'})

nicknames=['SWC', 'CP1', 'CP2', 'JAR', 'JR2', 'JR3',
'NAU', 'GIT', 'HUM', 'SUM', 'TUN', 'DY2',
'SDL', 'SDM', 'NAE', 'NSE', 'NGRP', 'NEM',
'EGP', 'KAR', 'KUL', 'AUR', 'PTG',
'PET']

meta["nickname"]=''

sites=['Swiss Camp', 'Crawford Pt. 1', 'CP2', 'JAR1', 'JAR2', 'JAR3',
'NASA-U', 'GITS', 'Humboldt', 'Summit', 'Tunu-N', 'DYE-2',
'Saddle', 'South Dome', 'NASA-E', 'NASA-SE', 'NGRIP', 'NEEM',
'EastGRIP', 'KAR', 'KULU', 'Aurora', 'Petermann Glacier',
'Petermann ELA']

for ss,site in enumerate(sites):
    meta["nickname"][meta["name"] == sites[ss]]=nicknames[ss]
    
names=['SMS-PET','SMS1','SMS2','SMS3','SMS4','SMS5','LAR1','LAR2','LAR3','Swiss Camp 10m','Roof_GEUS','LYN_T','LYN_L','CEN1','CEN2','KPC_Lv3','KPC_Uv3','THU_L2','WEG_B','ZAK_Uv3']
for name in names:
    meta.drop(meta[meta.name==name].index, inplace=True) # drop original time column

print(meta.name)
# print(meta.columns)


#%%
n_AWS=len(meta)

# meta=meta.sort_values(by='name', ascending=True)
# print(meta.columns)

min_tolerated_dist=20

iyear=1995 ; fyear=1995
iyear=1996 ; fyear=1999 ; year_prefix=1900
iyear=2001 ; fyear=iyear ; 
n_years=fyear-iyear+1
years=np.arange(iyear,fyear+1).astype(str)

years=[
       '1995',
'1996',
'1997',
'1998',
'1999',
'2001',
'2002',
'2003',
'2005',
'2006',
'2008',
]


years=['2010',
       '2011',
       '2012',
       '2013',
       '2014',
       '2015',
       '2016',
       '2017',
       '2019',
       ]

years=[
       '1995',
'1996',
'1997',
'1998',
'1999',
'2001',
'2002',
'2003',
'2005',
'2006',
'2008',

        '2010',
       '2011',
       '2012',
       '2013',
       '2014',
       '2015',
       '2016',
       '2017',
       '2019',
        ]

n_years=len(years)

for k,site in enumerate(meta.name):
# for k in range(1):
    # if k==0: # SWC
    # if k==3: # JAR
    # if k==4: # JAR2
    if k>=11: 
    # if k>=0:

        sentence_list=[]
        for year in years:
            print(year)
            if int(year)>=2000:
                year_prefix=2000
            else:
                year_prefix=1900
            
            path=f'/Users/jason/Dropbox/AWS/GCNET/GCNet_positions.stash/ATM/output/{year}_*'
            
            files = glob(path+"*")
            
            # print(files)
            
            n_files=len(files)
            
            for i,file in enumerate(files):
                # print(i,file)
                df=pd.read_csv(file)
                v=np.where(df.site==site)
                v=v[0]
                df.columns
                if df.dist.values[v]<5:
                    # print(site,v,df.dist.values[v])
                    varsx=[ 'year', 'month', 'day', 'elev_ATM', 'dist','slope_S2N', 'slope_W2E', 'lat', 'lon']
                    n_vars=len(varsx)
                    sentence=np.zeros(n_vars)
                    for vv,var in enumerate(varsx):
                        # print(var,df[var].values[v])
                        sentence[vv]=df[var].values[v]
                    sentence_list.append(sentence)
                # print(df[v[0]])
        if len(sentence_list)>0:
            sentence_list=np.array(sentence_list)
            out=pd.DataFrame({'year':sentence_list[:,0].astype(int),
                  'month':sentence_list[:,1].astype(int),
                  'day':sentence_list[:,2].astype(int),
                  'elev_ATM_m':sentence_list[:,3].astype(float),
                  'distance_km':sentence_list[:,4].astype(float),
                  'slope_S2N_deg':sentence_list[:,5].astype(float),
                  'slope_W2E_deg':sentence_list[:,6].astype(float),
                  'lat':sentence_list[:,7].astype(float),
                  'lon':sentence_list[:,8].astype(float),
                  })
            out['date']=pd.to_datetime(out[['year', 'month', 'day']])
            out = out.sort_values(by='date')

            # out.to_csv(f'./ATM/output/{site}_v2.csv',columns=['date', 'lat', 'lon','elev_ATM_m', 'distance_km','slope_S2N_deg', 'slope_W2E_deg'],index=None)
            out.to_csv(f'./ATM/output/{nicknames[k]}.csv',columns=['date', 'lat', 'lon','elev_ATM_m', 'distance_km','slope_S2N_deg', 'slope_W2E_deg'],index=None)
            print(out)
