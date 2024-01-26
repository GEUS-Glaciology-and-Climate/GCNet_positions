#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 12:16:39 2024

@author: jason
"""

from glob import glob
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
from pathlib import Path
import calendar
import simplekml

def save_kml(lat,lon,namex,opath,ofile):
    kml = simplekml.Kml(open=1)
    pnt = kml.newpoint(name=namex)
    pnt.coords=[(lon,lat)]
    pnt.style.labelstyle.color = simplekml.Color.red  # Make the text red
    pnt.style.labelstyle.scale = 1.5  # text scaling multiplier
    pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
    pnt.altitudemode = simplekml.AltitudeMode.relativetoground
    
    kml_ofile=opath+ofile+".kml"
    kml.save(kml_ofile)


# ## change to your system's login name to change dir for local work
if os.getlogin() == 'jason':
    base_path = '/Users/jason/Dropbox/AWS/GCNET/GCNet_positions.stash/'


fn='./meta/Stober_et_al_2023/SWC+ST2-Alle Pegel_Koordinaten Geodätisch.xlsx'
# os.system('open '+fn)
stober=pd.read_excel(fn,skiprows=7,names=['latdeg','latmin','latsec','NS','lat','emt','londeg','lonmin','lonsec','EW','lon','elevation'],index_col=0)
# stober = stober.iloc[:-32]
stober_dates=[]
stober_elevs=[]
for i,temp in enumerate(stober.index):
    temp=str(temp)
    temp.replace('106-ex010915','nan')
    temp.replace('ex010915','nan')
    temp.replace('ex010915','nan')
    # print(i,temp)
    if temp!='ST2' and temp!='nan' and temp[4:]!='ex010915':
        lat=stober.lat.values[i]
        lon=-stober.lon.values[i]
        osx=4
        if i>57:
            osx=6
        pegel=temp[0:osx-1]

        datex=pd.to_datetime(temp[osx:],format=('%d%m%y'))
        print(pegel,datex.strftime('%Y-%m-%d'),lat,lon,stober.elevation.values[i])
        # if pegel=='106':
        #     stober_dates.append(datex.strftime('%Y-%m-%d'))
        #     stober_elevs.append(stober.elevation.values[i]-1.08)#-25)
        save_kml(lat,lon,pegel+' '+datex.strftime('%Y-%m-%d'),'./output/kml/Stober_Swiss_Camp/',pegel+' '+datex.strftime('%Y-%m-%d'))
