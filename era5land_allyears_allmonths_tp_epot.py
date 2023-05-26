#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 11:16:55 2023

@author: evelien
"""

import numpy as np
import cdsapi
import os

lat_min_NL = 50.5
lat_max_NL = 54

lon_min_NL = 3.2
lon_max_NL = 7.2


lat_min_FIN =59.661825
lat_max_FIN = 65.195894
lon_min_FIN =  21.079833
lon_max_FIN =  32.090713


#%%


lon_min = lon_min_FIN
lon_max = lon_max_FIN

lat_min = lat_min_FIN
lat_max = lat_max_FIN

area_selection_shortname = 'FIN'
start_year = 1959
end_year = 2023

years = [str(year) for year in np.arange(start_year, end_year)]
proj = '/home/evelien/Agri-SiSa' 
proj_era_data = os.path.join(proj, 'data', 'ERA5')
filename= 'ERA5land_tp_ev_soilmoisture_mm_'+str(start_year)+'_' + str(end_year-1) + '_moda_'+ area_selection_shortname + '.zip'

if os.path.isdir(proj_era_data):
    print(proj_era_data, ' exists')
    filename_path = os.path.join(proj_era_data, filename)


download_variables = [
    'total_precipitation',
     'potential_evaporation'
     #'evaporation_from_bare_soil', 'evaporation_from_open_water_surfaces_excluding_oceans', 'evaporation_from_the_top_of_canopy',
    #'evaporation_from_vegetation_transpiration', , 'skin_reservoir_content',
    #'surface_latent_heat_flux', 'surface_net_solar_radiation', 'surface_solar_radiation_downwards',
    #'total_evaporation', 'total_precipitation', 'volumetric_soil_water_layer_1',
    #'volumetric_soil_water_layer_2', 'volumetric_soil_water_layer_3', 'volumetric_soil_water_layer_4',
]
#%%


c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-land-monthly-means',
    {
        'product_type': 'monthly_averaged_reanalysis',
        'stream':'moda',
        'variable': download_variables,
        'year': years,
        
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'time': '00:00',
        'format': 'netcdf.zip',
        'area': [
            lat_max, lon_min,
            lat_min, lon_max

        ],
    },
    filename)