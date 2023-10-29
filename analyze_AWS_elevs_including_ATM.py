#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 07:32:22 2023

@author: jason, jeb@geus.dk

compare ATM and other sources of AWS elevations

see ATM folder in this repository

"""

from glob import glob
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
from pathlib import Path

# ## change to your system's login name to change dir for local work
if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/'

os.chdir(base_path)

# read Greenland AWS locations
meta = pd.read_csv('/Users/jason/Dropbox/AWS/GCNET/GC-Net-level-1-data-processing/L1/GC-Net_location.csv')
meta = meta.rename({'Name': 'name'}, axis=1)
meta = meta.rename({'Latitude (°N)': 'lat'}, axis=1)
meta = meta.rename({'Longitude (°E)': 'lon'}, axis=1)
meta = meta.rename({'Elevation (wgs84 m)': 'elev'}, axis=1)

meta["nickname"]=''

sites=['Swiss Camp', 'Crawford Point 1', 'CP2', 'JAR1', 'JAR2', 'JAR3',
'NASA-U', 'GITS', 'Humboldt', 'Summit', 'Tunu-N', 'DYE-2',
'Saddle', 'South Dome', 'NASA-E', 'NASA-SE', 'NGRIP', 'NEEM',
'EastGRIP', 'KAR', 'KULU', 'Aurora', 'Petermann Glacier',
'Petermann ELA']

nicknames=['SWC', 'CP1', 'CP2', 'JR1', 'JR2', 'JR3',
'NAU', 'GIT', 'HUM', 'SUM', 'TUN', 'DY2',
'SDL', 'SDM', 'NAE', 'NSE', 'NGRP', 'NEM',
'EGP', 'KAR', 'KUL', 'AUR', 'PTG',
'PET']

for ss,site in enumerate(sites):
    meta["nickname"][meta["name"] == sites[ss]]=nicknames[ss]

names=['SMS-PET','SMS1','SMS2','SMS3','SMS4','SMS5','LAR1','LAR2','LAR3','Swiss Camp 10m','Roof_GEUS','LYN_T','LYN_L','CEN1','CEN2','KPC_Lv3','KPC_Uv3','THU_L2','WEG_B','ZAK_Uv3']
for name in names:
    meta.drop(meta[meta.name==name].index, inplace=True) # drop original time column

print(meta.name)
print(meta.columns)

#%%
df = pd.read_excel('/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/meta/GC-Net historical positions.xlsx')

df.drop(df[df.date == 'doc_2000'].index, inplace=True)
df.drop(df[df.date == 'ref'].index, inplace=True)
df.drop(df[df.date == 'WMO-DMI_2012'].index, inplace=True)
df.drop(df[df.date == 35596].index, inplace=True)
df.drop(df[df.elev == '-'].index, inplace=True)
df.drop(df[df.elev == np.nan].index, inplace=True)

df.loc[ df["name_long"] == "Crawford Pt. 1","name_long"] = "Crawford Point 1"

df.date=pd.to_datetime(df.date)
df['year']=df.date.dt.year
print(np.array(df.date))
print(len(df))
pos=df.copy()
#%%
n_AWS=len(meta)

# meta=meta.sort_values(by='name', ascending=True)
# print(meta.columns)


# iyear=1995 ; fyear=1998
# n_years=fyear-iyear+1
# years=np.arange(iyear,fyear+1).astype(str)

years=['1995',
'1996',
'1997',
'1998',
'1999',
'2001',
'2002',
'2003',
'2005',
'2006']

years=['2001',
'2002',
'2003',
'2005',
'2006']

years=['1995',
'1996',
'1997',
'1998',
'1999',
'2001',
'2002',
'2005',
'2006',
'2007',
'2008',
'2010',
       '2011',
       '2012',
       '2013',
       '2014',
       '2015',
       '2016',
       '2017',
       '2018',
       '2019',
]

n_years=len(years)

print(n_years)
#%%

th=1 
font_size=12
# plt.rcParams['font.sans-serif'] = ['Georgia']
plt.rcParams["font.size"] = font_size
plt.rcParams['axes.facecolor'] = 'w'
plt.rcParams['axes.edgecolor'] = 'k'
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 1
plt.rcParams['grid.color'] = "#cccccc"
plt.rcParams["legend.facecolor"] ='w'
plt.rcParams["mathtext.default"]='regular'
plt.rcParams['grid.linewidth'] = th
plt.rcParams['axes.linewidth'] = th #set the value globally
plt.rcParams['figure.figsize'] = 17, 10
plt.rcParams["legend.framealpha"] = 0.8
plt.rcParams['figure.figsize'] = 5, 4
ms=12

sites=['Swiss Camp', 'Crawford Point 1', 'CP2', 'JAR1', 'JAR2', 'JAR3',
'NASA-U', 'GITS', 'Humboldt', 'Summit', 'Tunu-N', 'DYE-2',
'Saddle', 'South Dome', 'NASA-E', 'NASA-SE', 'NGRIP', 'NEEM',
'EastGRIP', 'KAR', 'KULU', 'Aurora', 'Petermann Glacier',
'Petermann ELA']

choice_site=0 # SWC
choice_site=1 # CP1
# choice_site=3 # JAR1
# choice_site=4 # JAR2
# choice_site=6 # NAU
# choice_site=7 # GIT
# choice_site=8 # HUM
# choice_site=9 # SUM
# choice_site=10 # TUN
# choice_site=11 # DY2
# choice_site=12 # SDL
# choice_site=13 # SDM
# choice_site=14 # NAE
# choice_site=15 # NSE
# choice_site=23 # PTE

elevs=np.zeros(n_years)
yearx=np.zeros(n_years)
dist=np.zeros(n_years)

min_tolerated_dist=2
# n_AWS=1
for k in range(n_AWS):
# for k in [choice_site]:
    
    plt.close()
    plt.clf()
    fig, ax = plt.subplots(figsize=(9,7))
    yx=np.zeros(n_years)
    for yy,year in enumerate(years):
        df=pd.read_csv(f'./ATM/output/{year}.csv')
        # df.columns
        # if meta.ID.values[k]==choice_site:
        # if df.site.values[k]=='NASA-U':
        # v=np.where(df.site==sites[choice_site])
        v=np.where(df.site==sites[k])
        v=v[0]
        yearx[yy]=int(year)
        elevs[yy]=df.elev_ATM[v]
        dist[yy]=df.dist[v]
        
        # x1=np.sin(np.radians(df.slope_S2N.values[k]))*df.dist.values[k]*1000
        # x2=np.sin(np.radians(df.slope_W2E.values[k]))*df.dist.values[k]*1000
        # x1=df.slope_S2N.values[v][0]*df.dist.values[v][0]*1000
        # x2=df.slope_W2E.values[v][0]*df.dist[v][0]*1000
        x1=df.slope_S2N.values[v]*df.dist.values[v]*1000
        x2=df.slope_W2E.values[v]*df.dist[v]*1000
        yx[yy]=df.elev_ATM.values[v]-x1+x2
        if dist[yy]<min_tolerated_dist:
            # print(yearx[yy],"%.1f"%elevs[yy],"%.1f"%dist[yy],"%.1f"%(df.slope_S2N.values[v][0]*1000),"%.1f"%(df.slope_W2E.values[v][0]*1000),
            #       "%.1f"%x1,"%.1f"%x2,"%.1f"%(x1+x2))
            print(yearx[yy],"%.1f"%elevs[yy],"%.1f"%dist[yy],"%.1f"%(df.slope_S2N.values[v]*1000),"%.1f"%(df.slope_W2E.values[v]*1000),
                  "%.1f"%x1,"%.1f"%x2,"%.1f"%(x1+x2))    
    v=np.where(dist<min_tolerated_dist)
    v=v[0]
    
    x=yearx[v]
    y=elevs[v]
    y2=yx[v]

    xxx=1000
    plt.plot(x,y,'o', fillstyle='none',markersize=ms,label="ATM within %.1f"%min_tolerated_dist+" km: %.0f"%np.mean(y)+"±%.0f"%np.std(y)+' m')
    plt.plot(x,y2,'s', fillstyle='none',markersize=ms,label="ATM with slope cor: %.0f"%np.mean(y2)+"±%.0f"%np.std(y2)+' m')
    
    if len(x)>1:
        b, m = polyfit(x, y, 1)
        xx=[x[0],x[-1]]
        yy=[xx[0]*m+b,xx[1]*m+b]
        plt.plot(xx,yy,c='grey')
        # print(np.mean(y),np.std(y))
    plt.title(sites[k])
    plt.ylabel('altitude m above WGS84 ellipsoid')
    
    year0=1989 ; year1=2023
    plt.xlim(year0,year1)
    
    v=pos.name_long==sites[k]
    y3=pos.elev[v]
    if ~np.isnan(np.std(y3)):
        plt.plot(pos.year[v],y3,'s', fillstyle='none',markersize=ms,c='b',label="GC-Net historical positions.xlsx: %.0f"%np.mean(y3)+"±%.0f"%np.std(y3)+' m')
    
    fn = Path(f'/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/output/Jason/{nicknames[k]}_positions_monthly.csv')
    if fn.is_file():
        print(fn)
        GEUS_AWS_position=pd.read_csv(fn)
        GEUS_AWS_position.elev[GEUS_AWS_position.elev==0]=np.nan
        # GEUS_AWS_position['day']=15
        # GEUS_AWS_position["date"]=pd.to_datetime(GEUS_AWS_position[['year', 'month', 'day']])
        # df['year']=df.date.dt.year
        plt.plot(GEUS_AWS_position.year,GEUS_AWS_position.elev,'s', fillstyle='none',markersize=ms,c='r',
                 label="GEUS GC-Net GPS: %.0f"%np.mean(GEUS_AWS_position.elev)+"±%.0f"%np.std(GEUS_AWS_position.elev)+' m')

    plt.hlines(meta.elev.values[k],year0,year1,color='k',label="Table 4 Vandecrux er al 2023: %.0f"%meta.elev.values[k]+' m')

    plt.legend()
    
    ly='p'
    if ly =='x':plt.show()
    if ly =='p':
        plt.savefig(f'./ATM/Figs/{sites[k]}.png', bbox_inches='tight', dpi=150)