#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 08:40:10 2022

GCN_positions_from_field_notes_kml

@author: jason

"""

from pykml import parser
from glob import glob

site='PET'
site='DY2'

base_dir = '/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/meta/from_fieldbooks/'
files = sorted(glob(base_dir+site+'*.kml'))

print(files)

lats=[]
lons=[]
years=[]
months=[]
days=[]
for file in files:
    # print(file)
    root = parser.fromstring(open(file, 'rb').read())
    temp=str(root.Document.Placemark.Point.coordinates)
    lat=float(temp.split(',')[1])
    lon=float(temp.split(',')[0])
    lons.append(lon)
    lats.append(lat)

    year=file.split(' ')[-3]
    month=file.split(' ')[-2]
    day=file.split(' ')[-1].split('.')[0]
    # print()
    print(file,year,month,day,lat,lon)
    years.append(int(year))
    months.append(month)
    days.append(day)

import numpy as np
import pandas as pd

L = years
from operator import itemgetter
indices, L_sorted = zip(*sorted(enumerate(L), key=itemgetter(1)))

print(indices, L_sorted)

# years=years
#%%
years=np.array(years)
indices=np.array(indices)
months=np.array(months)
days=np.array(days)
lats=np.array(lats)
lons=np.array(lons)
if site=='PET': elevs=[902,915,924,944,907]
if site=='DY2': elevs=[2251,2080,2046]

nams=['year','month','day','lat','lon','elev']
df = pd.DataFrame(columns = nams)
df.index.name = 'index'
df["year"]=pd.Series(years[indices])
df["month"]=pd.Series(months[indices])
df["day"]=pd.Series(days[indices])
df["lat"]=pd.Series(lats[indices])
df['lat'] = df['lat'].apply(lambda x: '%.5f' % x)
df["lon"]=pd.Series(lons[indices])
df['lon'] = df['lon'].apply(lambda x: '%.5f' % x)
df["elev"]=pd.Series(elevs)

df.to_csv('/Users/jason/Dropbox/AWS/GCNET/GCNet_positions/meta/'+site+'.csv',index=None)

