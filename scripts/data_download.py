#!/usr/bin/env python
# coding: utf-8

# # Setup
# 
# ## Note:
# - Download of ERA5 data below will require around 80 GB of storage (plus another 330 GB if global analysis is desired), and may take several days due to queueing in the ECMWF Climate Data Store (CDS)
# - Download of GLASS data will require around 30 GB of storage
# - So a total of around 100 GB of storage will be required (plus another 330 GB if global analysis is desired)

# In[ ]:


# Import libraries
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import cdsapi
import wget


# # Download ERA5 data
# 
# If haven't already, first set up ECMWF CDS API using instructions from here: https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5#HowtodownloadERA5-4-DownloadERA5familydatathroughtheCDSAPI

# In[ ]:


# Open CDS API client for ERA5 downloads
c = cdsapi.Client()


# In[ ]:


# Areas for each region in NWSE format to retrieve data for
area = {
    "ca": [17, -91, 7, -81],
    "sa": [0, -65, -15, -30],
    "wa": [-26, 114, -36, 124],
    "global": [90, -180, -90, 180]
}


# In[ ]:


# Months to retrieve data for
months = [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
]

# Days to retrieve data for
days = [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
]

# Times to retrieve data for
times = [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
]


# In[ ]:


# Variable names in ECMWF CDS, available here:
# https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation

vars_era5_static = [
            'angle_of_sub_gridscale_orography', 'anisotropy_of_sub_gridscale_orography', 'geopotential',
            'high_vegetation_cover', 'lake_cover', 'lake_depth',
            'land_sea_mask', 'low_vegetation_cover', 'slope_of_sub_gridscale_orography',
            'soil_type', 'standard_deviation_of_filtered_subgrid_orography', 'standard_deviation_of_orography',
            'type_of_high_vegetation', 'type_of_low_vegetation',
]

vars_era5 = {
    "sfc": [
            '100m_u_component_of_wind', '100m_v_component_of_wind', '10m_u_component_of_wind',
            '10m_v_component_of_wind', '2m_temperature', 'mean_sea_level_pressure',
            'surface_latent_heat_flux', 'surface_sensible_heat_flux',
    ],
    "atm": [
            'evaporation', 'total_column_cloud_liquid_water', 'total_column_water_vapour',
            'vertical_integral_of_divergence_of_moisture_flux', 'vertical_integral_of_energy_conversion', 
            'vertical_integral_of_kinetic_energy', 'vertical_integral_of_potential_internal_and_latent_energy',
    ],
    "cld": [
            'boundary_layer_height', 'cloud_base_height', 'convective_available_potential_energy',
            'convective_inhibition', 'forecast_albedo', 'total_cloud_cover',
    ]
}


# ## Static data (global)

# ### Request all static variables

# In[ ]:


def retrieve_era5_slv_static_all():
    # Create global_era5-slv_static folder to download data into (if it doesn't already exist)
    Path("../data_raw/global_era5-slv_static").mkdir(parents=True, exist_ok=True)
    file_name = "../data_raw/global_era5-slv_static/global_era5-slv_static_all.nc"
    if Path(file_name).exists():
        print(file_name, "already exists")
    else:
        try:
            c.retrieve(
                'reanalysis-era5-single-levels-monthly-means',
                {
                    'product_type': 'monthly_averaged_reanalysis',
                    'variable': vars_era5_static,
                    'year': '2022',
                    'month': '01',
                    'time': '00:00',
                    'format': 'netcdf',
                },
                file_name)
            print("Retrieved", file_name)
        except:
            print("Failed to retrieve " + file_name)


# In[ ]:


# Retrieve all static data
retrieve_era5_slv_static_all()


# ## Monthly averaged reanalysis by hour of day, on single levels

# In[ ]:


# Define request to retrieve monthly reanalysis (to be used with retrieve function)
def request_era5_slv_month_hour(region, year, vars_type):
    file_name = '../data_raw/{region}_era5-slv-{vars_type}_month-hour/{region}_era5-slv-{vars_type}_month-hour_{year}.nc'.format(
        region=region, vars_type=vars_type, year=year)
    if Path(file_name).exists():
        print(file_name + " already exists")
    else:
        try:
            c.retrieve(
                'reanalysis-era5-single-levels-monthly-means',
                {
                    'product_type': 'monthly_averaged_reanalysis_by_hour_of_day',
                    'variable': vars_era5[vars_type],
                    'year': year,
                    'month': months,
                    'time': times,
                    'format': 'netcdf',
                    'area': area[region],
                },
                file_name)
            print("Retrieved " + file_name)
        except:
            print("Failed to retrieve " + file_name)
        
# Define function to retrieve monthly reanalysis    
def retrieve_era5_slv_month_hour(region, vars_type, year_start, year_end):
    # Assert region and vars_type is valid so we don't unneccessarily create a folder in the next part
    # And so we don't unnecessarily trigger the exception message
    assert region in area.keys(), f"region not one of: {*[*area],}"
    assert vars_type in [*vars_era5], f"vars_type not one of {[*vars_era5]}"
    # Create {region}_era5-slv-{vars_type}_month-hour folder to download data into
    # (if it doesn't already exist)
    Path("../data_raw/{region}_era5-slv-{vars_type}_month-hour".format(region=region, vars_type=vars_type)).mkdir(parents=True, exist_ok=True)
    # Run up to 10 parallel retrieve requests (ECMWF CDS only allows 1 year per request for hourly data)
    with ThreadPoolExecutor(max_workers=10) as executor:
        for year in range(year_start, year_end+1):
            executor.submit(request_era5_slv_month_hour, region, year, vars_type)


# In[ ]:


# Retrieve monthly averaged by hour reanalysis data for surface analysis in Western Australia
retrieve_era5_slv_month_hour("wa", "sfc", 1980, 2021)


# In[ ]:


# Retrieve monthly averaged by hour reanalysis data for atmospheric analysis in Western Australia
retrieve_era5_slv_month_hour("wa", "atm", 1980, 2021)


# In[ ]:


# Retrieve monthly averaged by hour reanalysis data for cloud analysis in Western Australia
retrieve_era5_slv_month_hour("wa", "cld", 1980, 2021)


# In[ ]:


# Retrieve monthly averaged by hour reanalysis data for surface analysis in Central America
retrieve_era5_slv_month_hour("ca", "sfc", 1980, 2021)


# In[ ]:


# Retrieve monthly averaged by hour reanalysis data for atmospheric analysis in Central America
retrieve_era5_slv_month_hour("ca", "atm", 1980, 2021)


# In[ ]:


# Retrieve monthly averaged by hour reanalysis data for cloud analysis in Central America
retrieve_era5_slv_month_hour("ca", "cld", 1980, 2021)


# In[ ]:


# Retrieve monthly averaged by hour reanalysis data for surface analysis in South America
retrieve_era5_slv_month_hour("sa", "sfc", 1980, 2021)


# In[ ]:


# Retrieve monthly averaged by hour reanalysis data for atmospheric analysis in South America
retrieve_era5_slv_month_hour("sa", "atm", 1980, 2021)


# In[ ]:


# Retrieve monthly averaged by hour reanalysis data for cloud analysis in South America
retrieve_era5_slv_month_hour("sa", "cld", 1980, 2021)


# ## Hourly reanalysis, on single levels

# In[ ]:


# Define request to retrieve hourly reanalysis (to be used with retrieve function)
def request_era5_slv_hour(region, year, vars_type):
    file_name = '../data_raw/{region}_era5-slv-{vars_type}_hour/{region}_era5-slv-{vars_type}_hour_{year}.nc'.format(
        region=region, vars_type=vars_type, year=year)
    if Path(file_name).exists():
        print(file_name + " already exists")
    else:
        try:
            c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': vars_era5[vars_type],
                'year': year,
                'month': months,
                'day': days,
                'time': times,
                'area': area[region],
                'format': 'netcdf',
            },
            file_name)
            print("Retrieved " + file_name)
        except:
            print("Failed to retrieve " + file_name)

# Define function to retrieve hourly reanalysis    
def retrieve_era5_slv_hour(region, vars_type, year_start, year_end):
    # Assert region and vars_type is valid so we don't unneccessarily create a folder in the next part
    # And so we don't unnecessarily trigger the exception message
    assert region in area.keys(), f"region not one of: {*[*area],}"
    assert vars_type in [*vars_era5], f"vars_type not one of {[*vars_era5]}"
    # Create {region}_era5-slv-{vars_type}_hour folder to download data into
    # (if it doesn't already exist)
    Path("../data_raw/{region}_era5-slv-{vars_type}_hour".format(region=region, vars_type=vars_type)).mkdir(parents=True, exist_ok=True)
    # Run up to 10 parallel retrieve requests (ECMWF CDS only allows 1 year per request for hourly data)
    with ThreadPoolExecutor(max_workers=10) as executor:
        for year in range(year_start, year_end+1):
            executor.submit(request_era5_slv_hour, region, year, vars_type)


# In[ ]:


# Retrieve hourly reanalysis data for surface analysis in Central America
retrieve_era5_slv_hour("ca", "sfc", 1980, 2021)


# In[ ]:


# Retrieve hourly reanalysis data for surface analysis in South America
retrieve_era5_slv_hour("sa", "sfc", 1980, 2021)


# In[ ]:


# Retrieve hourly reanalysis data for surface analysis in Western Australia
retrieve_era5_slv_hour("wa", "sfc", 1980, 2021)


# # Download GLASS data
# 
# First check that the naming conventions for the file urls are still up to date here: http://www.glass.umd.edu/Overview.html
# 
# Of the different components in the naming convention, the product_version and production_date is most likely to change so check these by browsing the downloads page here: http://www.glass.umd.edu/Download.html
# 
# Take care in that the production_date may vary for different years within the same dataset

# In[ ]:


# Define function to retrieve GLASS data
# Check on server that the url components are up to date (especially product_version)
def retrieve_glass_8_day(variable, data_source, year_start, year_end, production_date):
    # Define dictionaries according to GLASS dataset names (user may need to update these)
    variable_number = {"lai": "01", "fapar": "09"}
    data_source_number = {"modis": "01", "avhrr": "02"}
    product_version = {"modis": "V60", "avhrr": "V40"}
    dl_dir_path = {"modis": "MODIS/0.05D", "avhrr": "AVHRR"}
    # Assert arguments are valid so we don't unneccessarily create a folder in the next part
    # And so we don't unnecessarily trigger the exception message
    assert variable in variable_number.keys(), f"variable not one of: {*[*variable_number],}"
    assert data_source in data_source_number.keys(), f"data_source not one of: {*[*data_source_number],}"
    # Create global_glass-{variable}-{data_source}_8-day folder to download data into
    # (if it doesn't already exist)
    file_path = "../data_raw/global_glass-{variable}-{data_source}_8-day".format(
        variable=variable, data_source=data_source)
    Path(file_path).mkdir(parents=True, exist_ok=True)
    # Download data from the correct url, and to the correct file name
    for year in range(year_start, year_end+1):
        for day in range(1, 361+1, 8):
            file_name = file_path + "/global_glass-{variable}-{data_source}_8-day_{year}-{day:03}.hdf".format(
                variable=variable, data_source=data_source, year=year, day=day)
            dl_url = ("http://www.glass.umd.edu/{variable}/{dl_dir_path}/{year}/GLASS{variable_number}" + 
                      "B{data_source_number}.{product_version}.A{year}{day:03}.{production_date}.hdf").format(
                variable=variable.upper(), dl_dir_path=dl_dir_path[data_source], year=year, day=day,
                variable_number=variable_number[variable], data_source_number = data_source_number[data_source],
                product_version=product_version[data_source], production_date=production_date)
            if Path(file_name).exists():
                print(file_name + " already exists")
            else:
                try:
                    wget.download(dl_url, file_name)
                    print("Retrieved " + file_name)
                except:
                    print("Failed to retrieve " + file_name + "; check on server if there is missing data " +
                          "for the given date, and/or if the product_version and production_date for that " +
                          "file is correct")


# In[ ]:


# Retrieve LAI data derived from MODIS
# Check on server for the correct production_date to enter for each year
retrieve_glass_8_day("lai", "modis", 2000, 2019, "2022010")
retrieve_glass_8_day("lai", "modis", 2020, 2021, "2022138")


# In[ ]:


# Retrieve LAI data derived from AVHRR
# Check on server for the correct production_date to enter for each year
retrieve_glass_8_day("lai", "avhrr", 1981, 2017, "2019353")
retrieve_glass_8_day("lai", "avhrr", 2018, 2018, "2019358")


# In[ ]:


# Retrieve FAPAR data derived from MODIS
# Check on server for the correct production_date to enter for each year
retrieve_glass_8_day("fapar", "modis", 2000, 2020, "2022092")


# In[ ]:


# Retrieve FAPAR data derived from AVHRR
# Check on server for the correct production_date to enter for each year
retrieve_glass_8_day("fapar", "avhrr", 1982, 2015, "2019353")
retrieve_glass_8_day("fapar", "avhrr", 2016, 2018, "2019358")


# # Download other data

# ## Definition for the State Barrier Fence of Western Australia
# 
# - Download from https://data-downloads.slip.wa.gov.au/ExtractDownload/DownloadFile/187882
# - Information on dataset available here: https://catalogue.data.wa.gov.au/dataset/state-barrier-fence-dafwa-030/resource/2b867d37-26d3-4be0-9cbd-6969dc30df4d?inner_span=True

# In[ ]:


# Create wa_sbf_static folder to download data into
# (if it doesn't already exist)
Path("../data_raw/wa_sbfwa_static").mkdir(parents=True, exist_ok=True)


# ## NOAA Climate Indices
# 
# Catalog of indices available here: https://psl.noaa.gov/data/climateindices/list/
# 
# Indices used in this analysis are:
# - AMO: NOAA/PSL AMO Index (AMOI)
# - PDO: NOAA/PSL PDO Index (PDOI)
# - ENSO: NOAA/CPC Oceanic Nino Index (ONI)
# - IOD: NOAA/PSL Dipole Mode Index (DMI)
# - AAO/SAM: NOAA/CPC AAO Index (AAOI)
# - AO/NAM: NOAA/CPC AO Index (AOI)
# - NAO: NOAA/CPC NAO Index (NAOI)
# - EPO: NOAA/CPC EPO Index (EPOI)
# 
# Abbreviations:
# - NOAA: National Oceanic and Atmospheric Administration
# - PSL: Physical Sciences Laboratory
# - CPC: Climate Prediction Center
# - AMO: Atlantic Multidecadal Oscillation
# - PDO: Pacific Decadal Oscillation
# - ENSO: El Nino-Southern Oscillation
# - IOD: Indian Ocean Dipole
# - AAO/SAM: Antarctic Oscillation/Southern Annular Mode
# - AO/NAM: Arctic Oscillation/Northern Annular Mode
# - NAO: North Atlantic Oscillation
# - EPO: Eastern Pacific Oscillation

# In[ ]:


indices = {
    "amoi": {"output": "noaa-psl_amoi", "url": "https://psl.noaa.gov/data/correlation/amon.us.data"},
    "pdoi": {"output": "noaa-psl_pdoi", "url": "https://psl.noaa.gov/data/correlation/pdo.data"},
    "oni": {"output": "noaa-cpc_oni", "url": "https://psl.noaa.gov/data/correlation/oni.data"},
    "dmi": {"output": "noaa-psl_dmi", "url": "https://psl.noaa.gov/gcos_wgsp/Timeseries/Data/dmi.had.long.data"},
    "aaoi": {"output": "noaa-cpc_aaoi", "url": "https://psl.noaa.gov/data/correlation/aao.data"},
    "aoi": {"output": "noaa-cpc_aoi", "url": "https://psl.noaa.gov/data/correlation/ao.data"},
    "naoi": {"output": "noaa-cpc_naoi", "url": "https://psl.noaa.gov/data/correlation/nao.data"},
    "epoi": {"output": "noaa-cpc_epoi", "url": "https://psl.noaa.gov/data/correlation/epo.data"},
}


# In[ ]:


Path("../data_raw/global_noaa-climate-indices").mkdir(parents=True, exist_ok=True)
for index in [*indices]:
    file_name = "../data_raw/global_noaa-climate-indices/" + indices[index]["output"]
    if Path(file_name).exists():
        print(file_name + " already exists")
    else:
        try:
            wget.download(indices[index]["url"], file_name)
            print("Retrieved " + file_name)
        except:
            print("Failed to retrieve " + file_name + "; check on catalog if url has changed.")


# # (Deprecated)
# 
# Data downloads which were executed but were not used within the thesis project due to time constraints, change of scope, or finding that it was unnecessary.

# ## (Deprecated) Global monthly averaged reanalysis by hour of day

# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global surface analysis
# retrieve_era5_slv_month_hour("global", "sfc", 1981, 1985)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global surface analysis
# retrieve_era5_slv_month_hour("global", "sfc", 1992, 1996)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global surface analysis
# retrieve_era5_slv_month_hour("global", "sfc", 2007, 2011)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global surface analysis
# retrieve_era5_slv_month_hour("global", "sfc", 2017, 2021)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global atmospheric analysis
# retrieve_era5_slv_month_hour("global", "atm", 1981, 1985)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global atmospheric analysis
# retrieve_era5_slv_month_hour("global", "atm", 1992, 1996)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global atmospheric analysis
# retrieve_era5_slv_month_hour("global", "atm", 2007, 2011)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global atmospheric analysis
# retrieve_era5_slv_month_hour("global", "atm", 2017, 2021)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global cloud analysis
# retrieve_era5_slv_month_hour("global", "cld", 1981, 1985)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global cloud analysis
# retrieve_era5_slv_month_hour("global", "cld", 1992, 1996)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global cloud analysis
# retrieve_era5_slv_month_hour("global", "cld", 2007, 2011)


# In[ ]:


# # Retrieve monthly averaged by hour reanalysis data for global cloud analysis
# retrieve_era5_slv_month_hour("global", "cld", 2017, 2021)


# ## (Deprecated) BoM hourly observation data
# 
# Request from http://www.bom.gov.au/catalogue/data-feeds.shtml

# In[ ]:


# # Create wa_bom_hour folder to download data into
# # (if it doesn't already exist)
# Path("../data_raw/wa_bom_hour").mkdir(parents=True, exist_ok=True)


# ## (Deprecated) BoM minutely observation data
# 
# Request from http://www.bom.gov.au/catalogue/data-feeds.shtml

# In[ ]:


# # Create wa_bom_minute folder to download data into
# # (if it doesn't already exist)
# Path("../data_raw/wa_bom_minute").mkdir(parents=True, exist_ok=True)


# ## (Deprecated) Bunny Fence Experiment (2005-2007) data
# 
# Request from https://www.eol.ucar.edu/field_projects/bufex

# In[ ]:


# # Create wa_bufex folder to download data into
# # (if it doesn't already exist)
# Path("../data_raw/wa_bufex").mkdir(parents=True, exist_ok=True)


# ## (Deprecated) Original set of ERA5 variables which were retrieved
# 
# In the original set of data retrievals, the atmospheric variables included 'vertical_integral_of_divergence_of_cloud_frozen_water_flux' and 'vertical_integral_of_divergence_of_cloud_liquid_water_flux'. This is because it was thought that the vertical integral of divergence of *water vapour* flux needed to be derived by subtracting these two variables from 'vertical_integral_of_divergence_of_moisture_flux'. But correspondence with ECMWF specialist support confirmed that 'vertical_integral_of_divergence_of_moisture_flux' includes the divergence for water vapour only (i.e. does not include cloud liquid or frozen water fluxes).

# In[ ]:


# # Variable names in ECMWF CDS, available here:
# # https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation

# vars_era5_static = [
#             'angle_of_sub_gridscale_orography', 'anisotropy_of_sub_gridscale_orography', 'geopotential',
#             'high_vegetation_cover', 'lake_cover', 'lake_depth',
#             'land_sea_mask', 'low_vegetation_cover', 'slope_of_sub_gridscale_orography',
#             'soil_type', 'standard_deviation_of_filtered_subgrid_orography', 'standard_deviation_of_orography',
#             'type_of_high_vegetation', 'type_of_low_vegetation',
# ]

# vars_era5 = {
#     "sfc": [
#             '100m_u_component_of_wind', '100m_v_component_of_wind', '10m_u_component_of_wind',
#             '10m_v_component_of_wind', '2m_temperature', 'mean_sea_level_pressure',
#             'surface_latent_heat_flux', 'surface_sensible_heat_flux',
#     ],
#     "atm": [
#             'evaporation', 'total_column_cloud_liquid_water', 'total_column_water_vapour',
#             'vertical_integral_of_divergence_of_cloud_frozen_water_flux', 'vertical_integral_of_divergence_of_cloud_liquid_water_flux', 'vertical_integral_of_divergence_of_moisture_flux',
#             'vertical_integral_of_energy_conversion', 'vertical_integral_of_kinetic_energy', 'vertical_integral_of_potential_internal_and_latent_energy',
#     ]
# }


# ## (Deprecated) Request one static ERA5 variable at a time
# 
# Abandoned in favour of retrieving all static variables simultaneously (since these don't take up much storage anyway).

# In[ ]:


# # Define request to retrieve static data (to be used with retrieve function)
# def request_era5_slv_static(variable):
#     file_name = "../data_raw/global_era5-slv_static/global_era5-slv_static_{output_name}.nc".format(
#         output_name=variable.replace("_", "-"))
#     if Path(file_name).exists():
#         print(file_name, "already exists")
#     else:
#         try:
#             c.retrieve(
#                 'reanalysis-era5-single-levels-monthly-means',
#                 {
#                     'product_type': 'monthly_averaged_reanalysis',
#                     'variable': variable,
#                     'year': '2022',
#                     'month': '01',
#                     'time': '00:00',
#                     'format': 'netcdf',
#                 },
#                 file_name)
#             print("Retrieved", file_name)
#         except:
#             print("Failed to retrieve " + file_name)
            
# # Define function to retrieve static data 
# def retrieve_era5_slv_static(variables):
#     # Assert variables are valid so we don't unneccessarily create folders in the next part
#     # And so we don't unnecessarily trigger the exception message
#     assert all(variable in vars_era5_static for variable in variables), \
#         "variables not subset of: {}".format(vars_era5["static"])
#     # Create global_era5-slv_static folder to download data into
#     # (if it doesn't already exist)
#     Path("../data_raw/global_era5-slv_static").mkdir(parents=True, exist_ok=True)
#     # Run up to 10 parallel retrieve requests (to queue and download data for all variables simultaneously)
#     with ThreadPoolExecutor(max_workers=min(len(variables), 10)) as executor:
#         for variable in variables:
#             executor.submit(request_era5_slv_static, variable)


# In[ ]:


# # Retrieve static data for slope of sub-gridscale orography, geopotential (to plot elevation) and land-sea mask
# retrieve_era5_slv_static(["slope_of_sub_gridscale_orography", "geopotential", "land_sea_mask"])


# In[ ]:




