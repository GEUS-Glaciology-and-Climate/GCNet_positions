# -*- coding: utf-8 -*-
"""
@author: bav@geus.dk

tip list:
    %matplotlib inline
    %matplotlib qt
    import pdb; pdb.set_trace()
"""

import os
import pandas as pd

df_all = pd.DataFrame()
file_list = os.listdir('csv')
for f in file_list:
    print(f)
    tmp = pd.read_csv('csv/'+f)
    
    df_all = df_all.append(tmp,ignore_index=True)

df_all["date"] = df_all.apply(
    lambda x: pd.to_datetime(x["year"], format="%Y")
    + pd.DateOffset(days=x["day_of_year"] - 1, hours=x["decimal_hour"]),
    axis=1,
)

df_all = df_all.set_index('date')
df_all = df_all.resample('D').mean()

df_all.to_csv('gnss_jar_all.csv')