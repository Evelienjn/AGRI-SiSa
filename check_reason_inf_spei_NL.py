#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:09:57 2023

@author: evelien

the goal of this script is to check whether computing the 
spei is leading to issues: infinity
does it come from negative values?

"""

import rpy2

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from tqdm import tqdm
from utils import *
import xarray as xr
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from scipy import stats

#%%

r_time_series = robjects.r('ts')
SPEI_package = importr('SPEI')
r_spei_function = robjects.r['spei']
r_spi_function = robjects.r['spi'] 

tstart = '1990-01'
tfinal='2022-12'


data_path_era5land = '/home/evelien/Agri-SiSa/data/ERA5/ERA5land_tp_ev_soilmoisture_mm_1959_2022_NL.nc' 

area_code = 'NL' 
#select expver=1 (niet erat (prilimilarly), maar era5=validated!)
ds= xr.open_dataset(data_path_era5land).sel(time=(slice(tstart, tfinal)))
#%%

ds.tp.max()
ds.pev.max()
ds.pev.min()

ds['tp_mm'] = ds.tp * 1000
ds['pev_mm'] = ds.pev * 1000

ds['tp_mm'].plot.hist()
ds['pev_mm'].plot.hist()


#%%

ds['precp_plus_potev'] = ds.tp_mm + ds.pev_mm
ds['precp_minus_potev'] = ds.tp_mm - ds.pev_mm

#%%
# the log logistic distribution is uwse for a non-negative random variable.
# in the era5land datasat, negative values of pev occur, this means that the condensation
# would occur given the atmospheric conditions.
ds['precp_plus_potev'].plot.hist()
#ds['precp_minus_potev'].plot.hist()



#%%
# Check positive value of pev
    
df_pev_mm = ds.pev_mm.to_dataframe()
df_pev_mm.query('pev_mm >0')
#no potive pev values
df_tp_mm = ds.tp_mm.to_dataframe()

#%%
df_D2 = ds['precp_plus_potev'].to_dataframe()
df_D2.query('precp_plus_potev < 0')
df_D2['precp_plus_potev'].isna().value_counts()

df_D1 = ds['precp_minus_potev'].to_dataframe()
df_D1.query('precp_minus_potev < 0')
df_D1['precp_minus_potev'].isna().value_counts()

#%%

#inladen era5 gebaseerde SPEI
pad_nl_spei= '/home/evelien/Agri-SiSa/data/output/ERA5land_SPEI136_1990_2022_NL.nc' 
ds_spei = xr.open_dataset(pad_nl_spei)
#%%
df_spei6 = ds_spei['SPEI6'].to_dataframe()
df_spei1 = ds_spei['SPEI1'].to_dataframe()
df_spei3 = ds_spei['SPEI3'].to_dataframe()

#%%
ds_spei['SPEI6'].max()

df_spei6_NL_inf  = df_spei6[df_spei6['SPEI6'] == np.inf].to_csv('/home/evelien/Agri-SiSa/data/output/infinite_spei6_NL_2021.csv', sep = ';')



#join outer ones

df = df_D1.join(other=df_D2,
                         how='left')

df = df.join(other=df_pev_mm,
             how = 'left')

df = df.join(other = df_tp_mm,
             how = 'left')


df = df.join(other = df_spei1,
             how='left')
df = df.join(other = df_spei3,
             how='left')

df = df.join(other = df_spei6,
             how ='left')


#%%

df= df.reset_index(names = ['datetime', 'latitude', 'longitude', 'expver'])
df = df.query('expver == 1')
test = df.query('datetime == "2021-07-01"' )
#%%
#selecteer die met inf SPEI6 en koppel die met 5 maanden ervoor en dezelfde latlonpaar eraan.

test_inf = df.query('SPEI6 == inf' )

#%%
ds_test = ds.sel(expver =1, latitude=50.7, longitude=6.6, method='nearest')
ds_test['precp_plus_potev'].plot()
##ds_test.pev_mm.plot()
ds_test.tp_mm.plot()
df_test = ds_test.to_dataframe()
df_test.reset_index(names=['time','latitude', 'longitude'])
test=df_test.query('time > "2021-01-01"' )[['latitude', 'longitude', 'pev_mm', 'tp_mm', 'precp_plus_potev']]
test.to_csv('/home/evelien/Agri-SiSa/data/output/NL_julys.csv', sep = ';')


julys = df_test.reset_index(names=['time'])
julys = julys.loc[julys['time'].dt.month == 7][['time','latitude', 'longitude', 'pev_mm', 'tp_mm', 'precp_plus_potev']]
julys.to_csv('/home/evelien/Agri-SiSa/data/output/NL_julys.csv', sep = ';')


#%%

dist = getattr(stats, 'fisk')
(c, loc, scale)= dist.fit(julys['precp_plus_potev'])
print(parameters)
ktest_statistic = stats.kstest(rvs = julys['precp_plus_potev'], #random variable sample                               
                               cdf = "fisk", # distribution to be tested against
                               args= [c, loc, scale] #distribution parameters of the fish distribution
                               ).statistic 


 #Mean('m'), variance('v'), skew('s'), and/or kurtosis('k').
alpha = 0.240
n = 33

critical_value_p005 = 0.231

if ktest_statistic  > critical_value_p005:
    print('Reject H0')
else:
    print('pass H0')    
    
#%%
import numpy as np

from scipy.stats import fisk

import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)

rv = fisk(c)
x = np.linspace(-10, 10, 100)
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
ax.set_xlim(-3,6)
ax.hist(julys['precp_plus_potev'])

#%%
dist = getattr(stats, 'gamma')
paramter = (c, loc, scale)= dist.fit(julys['precp_plus_potev'])
print(parameters)
ktest_statistic = stats.kstest(rvs = julys['precp_plus_potev'], #random variable sample                               
                               cdf = "gamma", # distribution to be tested against
                               args= parameters #distribution parameters of the fish distribution
                               ).statistic 


 #Mean('m'), variance('v'), skew('s'), and/or kurtosis('k').
alpha = 0.240
n = 33

critical_value_p005 = 0.231

if ktest_statistic  > critical_value_p005:
    print('Reject H0')
else:
    print('pass H0')    
    
#%%
import numpy as np

from scipy.stats import gamma
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)

rv = genlogistic(c)
x = np.linspace(-10, 10, 100)
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
ax.hist(julys['precp_plus_potev'])


#%%
dist = getattr(stats, 'genhyperbolic')
parameters= (a, b, loc, scale, size)= dist.fit(julys['precp_plus_potev'])
print(parameters)
ktest_statistic = stats.kstest(rvs = julys['precp_plus_potev'], #random variable sample                               
                               cdf = "genhyperbolic" , # distribution to be tested against
                               args= parameters #distribution parameters of the fish distribution
                               ).statistic 


 #Mean('m'), variance('v'), skew('s'), and/or kurtosis('k').
alpha = 0.240
n = 33

critical_value_p005 = 0.231

if ktest_statistic  > critical_value_p005:
    print('Reject H0')
else:
    print('pass H0')    
    
#%%
import numpy as np

from scipy.stats import genhyperbolic
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)

rv = genhyperbolic(a=a, b=b, p=-0.3)
x = np.linspace(-10, 10, 100)
ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
#ax.hist(julys['precp_plus_potev'])



    