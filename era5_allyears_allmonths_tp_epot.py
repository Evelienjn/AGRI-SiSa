#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:30:21 2023

@author: evelien

data download
"""
import numpy as np
import cdsapi
import os
#%%

start_year = 1959
end_year = 2022

years = [str(year) for year in np.arange(start_year, end_year)]
proj = '/home/evelien/Agri-SiSa' 
proj_era_data = os.path.join(proj, 'data', 'ERA5')
filename= 'ERA5_tp_pev_'+str(start_year)+'_' + str(end_year) + '_NL.nc'

if os.path.isdir(proj_era_data):
    print(proj_era_data, ' exists')
    filename_path = os.path.join(proj_era_data, filename)
    #%%

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels-monthly-means',
    {
        'product_type': 'monthly_averaged_reanalysis',
        'variable': [
            'potential_evaporation', 'total_precipitation',
        ],
        'year': years,
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        'area': [
            54, 3.2, 50.5,
            7.2,
        ],
        'format': 'netcdf',
    },
    filename_path)

