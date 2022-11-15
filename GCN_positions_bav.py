# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
# import matplotlib.colormap as cm

try:
    url = "https://docs.google.com/spreadsheets/d/1R2SA7rqo9PHfAAGeSVgy7eWVHRugV8Z3nbWga5Xin1U/export?format=csv&gid=0"
    pd.read_csv(url).to_csv("data/GC-Net_yearly_positions.csv", index=None)
except:
    pass

df_pos = pd.read_csv("data/GC-Net_yearly_positions.csv")
plt.close('all')
# col = cm('Spectral',32)
abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
fig, ax=plt.subplots(4,4, figsize=(15,9))
plt.subplots_adjust(left=0.08, right=0.99,
                    bottom=0.08,top=0.95,
                    hspace = 0.35, wspace=0.25)
ax=ax.flatten()
i=0
for j in [1,2,3,4,5,6,7,8,9,10, 11,12,15,22,23,24]:
    i = i+1
    tmp = df_pos.loc[df_pos.id==j,:]

    
    ax[i-1].plot(tmp.lon, tmp.lat, marker = 'o', linestyle='None')
    texts = []
    for k in range(tmp.shape[0]):
        lat, lon, date = tmp.iloc[k,:][['lat','lon','date']].to_list()
        texts.append(ax[i-1].text(lon, lat, date, fontsize=6))
    ax[i-1].set_title(abc[i-1]+'. '+tmp.name_long.unique()[0])
    c_lat = tmp.lat.mean()    
    c_lon = tmp.lon.mean()
    w = 0.07
    if tmp.lon.std() < w:
        ax[i-1].set_xlim(c_lon-w*2,c_lon+w*2)
    if tmp.lat.std() < w:
        ax[i-1].set_ylim(c_lat-w,c_lat+w)
    adjust_text(texts, ax=ax[i-1], arrowprops=dict(arrowstyle='-', color='lightgray'))

        