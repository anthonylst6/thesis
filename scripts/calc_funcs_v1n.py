#!/usr/bin/env python
# coding: utf-8

# In[ ]:


## Import libraries for calculations

import logging
import inspect
import copy
import math
import pandas as pd
import numpy as np
import xarray as xr
import metpy.calc as mpcalc
from glob import glob
from pathlib import Path
from datetime import datetime
from dateutil.relativedelta import relativedelta
from scipy.special import gamma
from scipy.interpolate import interp1d


# In[ ]:


## Settings and global variables for calculations

# This is to get the version number of the calc_funcs script being used so that it
# can be appended to the file name of any outputs. The reason this is done is because
# the calculation functions below and the plot functions in the plot_funcs script both
# output intermediate files one at a time from low level to high level, and that each
# file name is used in recognising whether there is a need to re-run a script (if the
# file already exists then the script is not run so as to save on computation).
# However, this method can propagate errors from low level through to high level
# if there has been a change to the code and/or output at the lower levels. By
# appending the version number of the calc_funcs script being used, it forces all
# intermediate files to be recreated from scratch rather than reuse intermediate files
# which was outputted by outdated code. "v00" is used as a placeholder version number
# if there is an error: it is used mostly for scripting purposes within an
# interactive python notebook where the file name cannot be directly extracted
# using the __file__ python variable.
try:
    calc_funcs_ver = "cf" + Path(__file__).stem[-3:]
except:
    calc_funcs_ver = "cfv00"

# Set level of logging out of (in decreasing order of detail): 
# [logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL]
calc_log_level = logging.INFO
log_levels = [logging.DEBUG, logging.INFO, logging.WARNING, 
              logging.ERROR, logging.CRITICAL]
assert calc_log_level in log_levels, \
    f"calc_log_level (global variable in settings section) must be one of: {log_levels}"

# Set priority of code to one of: ["speed", "memory"]. "speed" sacrifices additional 
# memory for faster processing while "memory" sacrifices faster processing for lower 
# RAM consumption. 
priority = "memory"
priorities = ["speed", "memory"]
assert priority in priorities, \
    f"priority (global variable in settings section) must be one of: {priorities}"

# Valid regions and their mapping to axis extents in [W, E, S, N] format 
# as well as timezones in hour +- GMT
regions = {
    # Central America (mostly Honduras-Nicaragua-Costa Rica)
    "ca": {"extent": [-91, -81, 7, 17], "tz": -6},
    # South America (mostly central and eastern Brazil)
    "sa": {"extent": [-65, -30, -15, 0], "tz": -3},
    # Western Australia (mostly near the west coast)
    "wa": {"extent": [113, 123, -35, -30], "tz": +8},
    # Global (requires a lot of memory)
    "global": {"extent": [-180, 180, -90, 90], "tz": +0}
}

# Size of chunks
chunksize = "500MB"

# Valid subset strings to use as argument in climatologies and
# their mapping to month numbers for use in xarray time slicing
months_subsets = {
    "all": list(range(1, 12+1)),
    "djf": [12,1,2], "mam": [3,4,5], "jja": [6,7,8], "son": [9,10,11]
}

# File number check to make sure data_download notebook was run correctly
number_of_glass_files = {"lai": {"avhrr": 1748, "modis": 1005},
                         "fapar": {"avhrr": 1702, "modis": 960}
                        }
number_of_era5_month_hour_files = 42
number_of_era5_hour_files = 42

# GLASS data sources
glass_sources_all = ["avhrr", "modis"]

# Earliest and latest entries in each GLASS dataset
avhrr_earliest = "Jan-1981"
modis_earliest = "Mar-2000"
avhrr_latest = "Dec-2018"
modis_latest = "Dec-2021"
fapar_earliest = "Jan-1982"
fapar_latest = "Dec-2020"

assert (datetime.strptime(fapar_earliest, "%b-%Y") > 
        datetime.strptime(avhrr_earliest, "%b-%Y")), \
    "fapar_earliest must be later than avhrr_earliest since this was a design assumption"
assert (datetime.strptime(fapar_latest, "%b-%Y") < 
        datetime.strptime(modis_latest, "%b-%Y")), \
    "fapar_latest must be earlier than modis_latest since this was a design assumption"
assert (datetime.strptime(modis_earliest, "%b-%Y") > 
        datetime.strptime(avhrr_earliest, "%b-%Y")), \
    "modis_earliest must be later than avhrr_earliest since this was a design assumption"
assert (datetime.strptime(modis_latest, "%b-%Y") > 
        datetime.strptime(avhrr_latest, "%b-%Y")), \
    "modis_latest must be later than avhrr_latest since this was a design assumption"

avhrr_earliest_year = int(avhrr_earliest[-4:])
avhrr_latest_year = int(avhrr_latest[-4:])
modis_earliest_year = int(modis_earliest[-4:])
modis_latest_year = int(modis_latest[-4:])
fapar_earliest_year = int(fapar_earliest[-4:])
fapar_latest_year = int(fapar_latest[-4:])

# Resolution of ERA5 dataset in degrees (used for regridding)
res_era5 = 0.25

# Valid ERA5/ERA5-derived parameters for use in analysis
hours_all = list(range(0, 23+1))
vars_and_dvars_era5 = {
    "vars": {
        "sfc": ["u10", "v10", "ws10", "wv10", "u100", "v100", "ws100", "wv100", 
                "mslp", "t2", "slhf", "sshf"],
        "atm": ["nse", "vidmf", "viec", "vipile", "vike", "tcclw", "tcwv", "nac"],
        "cld": ["blh", "fa", "cbh", "tcc", "cape", "ci"]
    },
    "dvars": {
        "sfc": ["du10", "dv10", "dws10", "dwv10", "du100", "dv100", "dws100", "dwv100", 
                "dmslp", "dt2", "dslhf", "dsshf"],
        "atm": ["dnse", "dvidmf", "dviec", "dvipile", "dvike", "dtcclw", "dtcwv", "dnac"],
        "cld": ["dblh", "dfa", "dcbh", "dtcc", "dcape", "dci"]
    }
}

vars_era5_all = []
for _, var_list in vars_and_dvars_era5["vars"].items():
    for var in var_list:
        vars_era5_all.append(var)

dvars_era5_all = []
for _, dvar_list in vars_and_dvars_era5["dvars"].items():
    for dvar in dvar_list:
        dvars_era5_all.append(dvar)

vars_and_dvars_era5_all = vars_era5_all + dvars_era5_all

# Valid time strings to use as argument in plot_funcs script,
times = {
    "0-5": list(range(0, 5+1)), "6-11": list(range(6, 11+1)),
    "12-17": list(range(12, 17+1)), "18-23": list(range(18, 23+1)),
    "night": list(range(0, 5+1)), "morning": list(range(6, 11+1)),
    "afternoon": list(range(12, 17+1)), "evening": list(range(18, 23+1))
}

# Parameters which are vectors
params_vector = ["wv10", "wv100", "dwv10", "dwv100"]

# Output parameters for calc_glass_mean_clim
params_glass_mean = ["mlai", "mfapar"]

# Output parameters for calc_era5_mdp_clim_stats_given_var_or_dvar
params_stat = ["hour_max", "hour_min", "max", "max_u", "max_v", 
               "min", "min_u", "min_v", "mean", "mean_u", "mean_v", "range"]

# Output parameters for calc_era5_wsd_clim
params_wsd = ["ws10_mean", "ws10_std", "c10", "k10", "ws100_mean", 
              "ws100_std", "c100", "k100", "eroe100", "tgcf100"]

# Output parameters for static orographic calculations
params_orog = ["lse", "ssgo"]

# Values which arg_extra in plot_funcs script can take on.
arg_extra_all = params_glass_mean + hours_all + params_stat + params_wsd

# Mapping from ERA5 dataset variable names to own desired names, while
# also accounting for any var_or_dvar variable dependencies
vars_deps_and_rename = {"u10": {"u10": "u10"},
                        "v10": {"v10": "v10"},
                        "ws10": {"u10": "u10", "v10": "v10"},
                        "wv10": {"u10": "u10", "v10": "v10"},
                        "u100": {"u100": "u100"},
                        "v100": {"v100": "v100"},
                        "ws100": {"u100": "u100", "v100": "v100"},
                        "wv100": {"u100": "u100", "v100": "v100"},
                        "mslp": {"msl": "mslp"},
                        "t2": {"t2m": "t2"},
                        "slhf": {"slhf": "slhf"},
                        "sshf": {"sshf": "sshf"},
                        "nse": {"e": "nse"},
                        "vidmf": {"p84.162": "vidmf"},
                        "viec": {"p64.162": "viec"},
                        "vipile": {"p62.162": "vipile"},
                        "vike": {"p59.162": "vike"},
                        "tcclw": {"tclw": "tcclw"},
                        "tcwv": {"tcwv": "tcwv"},
                        "nac": {"e": "nse", "p84.162": "vidmf", "tcwv": "tcwv"},
                        "blh": {"blh": "blh"},
                        "fa": {"fal": "fa"},
                        "cbh": {"cbh": "cbh"},
                        "tcc": {"tcc": "tcc"},
                        "cape": {"cape": "cape"},
                        "ci": {"cin": "ci"}
                       }

# Speed (in m/s) for expected rate of exceedance analysis at 100 m
speed_eroe = 42.5

# Typical power curve for a 100 m turbine with 100 m rotor diameter and 
# a nameplate rating of around 2500 kW (used to compute gross capacity factor)
# Speeds are in m/s, powers are in kW, data from https://www.thewindpower.net
speeds_common = np.append(np.linspace(0, 25.5, 52), 999)
power_nameplate = 2500
# Vestas V100/2600
powers_vestas = np.array([0, 0, 0, 0, 0, 0, 21, 63, 115, 172, 239, 318, 405, 550,
                          706, 890, 1080, 1283, 1485, 1641, 1796, 1944, 2092, 2225,
                          2351, 2440, 2502, 2560, 2584, 2597, 2600, 2600, 2600,
                          2600, 2600, 2600, 2600, 2600, 2600, 2600, 2600, 2600,
                          2600, 2600, 2600, 2600, 2600, 2600, 2600, 2600, 2600, 0, 0])
# Goldwind GW100/2500
powers_gw = np.array([0, 0, 0, 0, 0, 6, 34, 65, 101, 165, 235, 320, 409, 530, 655,
                      826, 997, 1196, 1394, 1669, 1943, 2170, 2313, 2415, 2458,
                      2485, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500,
                      2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500,
                      2500, 2500, 2500, 2500, 2500, 2500, 0, 0])
# GE Energy 2.5-100
powers_ge = np.array([0, 0, 0, 0, 0, 0, 10, 80, 160, 250, 340, 460, 590, 770, 952,
                      1170, 1389, 1650, 1869, 2100, 2260, 2400, 2487, 2500, 2500,
                      2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500,
                      2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500,
                      2500, 2500, 2500, 2500, 2500, 2500, 0, 0])
# Average
powers_avg = (powers_vestas/26*25 + powers_gw + powers_ge) / 3

# Names of main functions
calc_func_names = ["calc_glass_mean_clim",
                   "calc_era5_mdp_clim_given_var_or_dvar",
                   "calc_era5_mdp_clim_stats_given_var_or_dvar",
                   "calc_era5_wsd_clim"]
extra_func_names = ["calc_diff",
                    "calc_era5_orog",
                    "calc_glass_rolling_avg_of_annual_diff"]
main_func_names = calc_func_names + extra_func_names

# Maximum percentage of annual difference files used in computing glass rolling avg
# of annual difference which can be NaN DataArray's
per_diff_nan_max = 25

# Plot-specific arguments for use in the plot_funcs script which can take on None
# values (these args are excepted from the check_args_for_none function below)
args_plot = ["mask_period1", "mask_period2", "extents", "vmin", "vmax", 
             "vmin_periods", "vmax_periods", "vmin_diff", "vmax_diff", 
             "ax", "ax_period1", "ax_period2", "ax_diff", "cfv_data"]


# In[ ]:


# Attributes for DataArray's within output datasets.
attrs_da = {
    
    # For calc_era5_glass_mean
    "mlai": {"abbreviation": "MLAI",
             "full_name": "Mean Leaf Area Index",
             "units": "dimensionless", 
             "source": ""},
    "mfapar": {"abbreviation": "MFAPAR",
               "full_name": ("Mean Fraction of Absorbed Photosynthetically " +
                             "Active Radiation"),
               "units": "dimensionless",
               "source": ""},
    
    # For calc_era5_mdp_clim_given_var_or_dvar
    "u10": {"abbreviation": "U10",
            "full_name": "Zonal Component of Wind Velocity at 10 m Above Surface",
            "units": "$m s^{-1}$"},
    "v10": {"abbreviation": "V10",
            "full_name": "Meridional Component of Wind Velocity at 10 m Above Surface",
            "units": "$m s^{-1}$"},
    "ws10": {"abbreviation": "WS10",
             "full_name": "Wind Speed at 10 m Above Surface",
             "units": "$m s^{-1}$"},
    "wv10": {"abbreviation": "WV10",
             "full_name": "Wind Velocity at 10 m Above Surface",
             "units": "$m s^{-1}$"},
    "u100": {"abbreviation": "U100",
             "full_name": "Zonal Component of Wind Velocity at 100 m Above Surface",
             "units": "$m s^{-1}$"},
    "v100": {"abbreviation": "V100",
             "full_name": "Meridional Component of Wind Velocity at 100 m Above Surface",
             "units": "$m s^{-1}$"},
    "ws100": {"abbreviation": "WS100",
              "full_name": "Wind Speed at 100 m Above Surface",
              "units": "$m s^{-1}$"},
    "wv100": {"abbreviation": "WV100",
              "full_name": "Wind Velocity at 100 m Above Surface",
              "units": "$m s^{-1}$"},
    "mslp": {"abbreviation": "MSLP",
             "full_name": "Mean Sea Level Pressure",
             "units": "$Pa$"},
    "t2": {"abbreviation": "T2",
           "full_name": "Temperature at 2 m Above Surface",
           "units": "$K$"},
    "slhf": {"abbreviation": "SLHF",
             "full_name": "Surface Latent Heat Flux",
             "units": "$W m^{-2}$"},
    "sshf": {"abbreviation": "SSHF",
             "full_name": "Surface Sensible Heat Flux",
             "units": "$W m^{-2}$"},
    "viec": {"abbreviation": "VIEC",
             "full_name": "Vertical Integral of Energy Conversion",
             "units": "$W m^{-2}$"},
    "vipile": {"abbreviation": "VIPILE",
               "full_name": "Vertical Integral of Potential, Internal and Latent Energy",
               "units": "$J m^{-2}$"},
    "vike": {"abbreviation": "VIKE",
             "full_name": "Vertical Integral of Kinetic Energy",
             "units": "$J m^{-2}$"},
    "tcclw": {"abbreviation": "TCCLW",
              "full_name": "Total Column Cloud Liquid Water",
              "units": "$kg m^{-2}$"},
    "tcwv": {"abbreviation": "TCWV",
             "full_name": "Total Column Water Vapour",
             "units": "$kg m^{-2}$"},
    "nse": {"abbreviation": "NSE",
            "full_name": "Net Surface Evaporation",
            "units": "$kg m^{-2} s^{-1}$"},
    "vidmf": {"abbreviation": "VIDMF",
              "full_name": "Vertical Integral of Divergence of Moisture Flux",
              "units": "$kg m^{-2} s^{-1}$"},
    "vidcfwf": {"abbreviation": "VIDCFWF",
                "full_name": "Vertical Integral of Divergence of Cloud Frozen Water Flux",
                "units": "$kg m^{-2} s^{-1}$"},
    "vidclwf": {"abbreviation": "VIDCLWF",
                "full_name": "Vertical Integral of Divergence of Cloud Liquid Water Flux",
                "units": "$kg m^{-2} s^{-1}$"},
    "nac": {"abbreviation": "NAC",
            "full_name": "Net Atmospheric Condensation",
            "units": "$kg m^{-2} s^{-1}$"},
    "blh": {"abbreviation": "BLH",
            "full_name": "Boundary Layer Height",
            "units": "$m$"},
    "fa": {"abbreviation": "FA",
           "full_name": "Forecast Albedo",
           "units": "dimensionless"},
    "cbh": {"abbreviation": "CBH",
            "full_name": "Cloud Base Height",
            "units": "$m$"},
    "tcc": {"abbreviation": "TCC",
            "full_name": "Total Cloud Cover",
            "units": "dimensionless"},
    "cape": {"abbreviation": "CAPE",
             "full_name": "Convective Available Potential Energy",
             "units": "$J kg^{-1}$"},
    "ci": {"abbreviation": "CI",
           "full_name": "Convective Inhibition",
           "units": "$J kg^{-1}$"},
    
    # For calc_era5_mdp_clim_stats_given_var_or_dvar
    "hour_max": {"abbreviation": "$hour_{{max}}$({})",
                 "full_name": "Hour of Maximum for {}",
                 "units": "UTC {}"},
    "hour_min": {"abbreviation": "$hour_{{min}}$({})",
                 "full_name": "Hour of Minimum for {}",
                 "units": "UTC {}"},
    "max": {"abbreviation": "$max$({})",
            "full_name": "Maximum for {}",
            "units": "{}"},
    "max_u": {"abbreviation": "$max_u$({})",
              "full_name": "Zonal Component of Maximum for {}",
              "units": "{}"},
    "max_v": {"abbreviation": "$max_v$({})",
              "full_name": "Meridional Component of Maximum for {}",
              "units": "{}"},
    "min": {"abbreviation": "$min$({})",
            "full_name": "Minimum for {}",
            "units": "{}"},
    "min_u": {"abbreviation": "$min_u$({})",
              "full_name": "Zonal Component of Minimum for {}",
              "units": "{}"},
    "min_v": {"abbreviation": "$min_v$({})",
              "full_name": "Meridional Component of Minimum for {}",
              "units": "{}"},
    "mean": {"abbreviation": "$mean$({})",
             "full_name": "Mean for {}",
             "units": "{}"},
    "mean_u": {"abbreviation": "$mean_u$({})",
               "full_name": "Zonal Component of Mean for {}",
               "units": "{}"},
    "mean_v": {"abbreviation": "$mean_v$({})",
               "full_name": "Meridional Component of Mean for {}",
               "units": "{}"},
    "range": {"abbreviation": "$range$({})",
              "full_name": "Range for {}",
              "units": "{}"},
    
    # For calc_era5_wsd_clim
    "ws10_mean": {"abbreviation": "$mean$(WS10)",
                  "full_name": "Mean of Wind Speed at 10 m Above Surface",
                  "units": "$m s^{-1}$"},
    "ws10_std": {"abbreviation": "$std$(WS10)",
                 "full_name": "Standard Deviation of Wind Speed at 10 m Above Surface",
                 "units": "$m s^{-1}$"},
    "c10": {"abbreviation": "C10",
            "full_name": ("Scale Parameter of Wind Speed Weibull Distribution at " +
                          "10 m Above Surface"),
            "units": "$m s^{-1}$"},
    "k10": {"abbreviation": "K10",
            "full_name": ("Shape Parameter of Wind Speed Weibull Distribution at " +
                          "$10 m Above Surface$"),
            "units": "dimensionless"},
    "ws100_mean": {"abbreviation": "$mean$(WS100)",
                   "full_name": "Mean of Wind Speed at 100 m Above Surface",
                   "units": "$m s^{-1}$"},
    "ws100_std": {"abbreviation": "$std$(WS100)",
                  "full_name": "Standard Deviation of Wind Speed at 100 m Above Surface",
                  "units": "$m s^{-1}$"},
    "c100": {"abbreviation": "C100",
             "full_name": ("Scale Parameter of Wind Speed Weibull Distribution at " +
                           "100 m Above Surface"),
             "units": "$m s^{-1}$"},
    "k100": {"abbreviation": "K100",
             "full_name": ("Shape Parameter of Wind Speed Weibull Distribution at " +
                           "100 m Above Surface"),
             "units": "dimensionless"},
    "eroe100": {"abbreviation": "EROE100",
                "full_name": (f"Expected Rate of Exceeding {speed_eroe} m/s for Wind " +
                              "Speed Weibull Distribution at 100 m Above Surface"),
                "units": "dimensionless"},
    "tgcf100": {"abbreviation": "TGCF100",
                "full_name": ("Gross Capacity Factor for Typical Turbine at 100 m " +
                              "Above Surface given the Wind Speed Weibull Distribution " +
                              "at 100 m Above Surface"),
                "units": "dimensionless"},
    
    # For calc_era5_orog
    "lse": {"abbreviation": "LSE",
            "full_name": "Land Surface Elevation",
            "units": "$m$"},
    "ssgo": {"abbreviation": "SSGO",
             "full_name": "Slope of Sub-Gridscale Orography",
             "units": "dimensionless"}

}

# Attributes for coordinates within output datasets.
coord_attrs = {
    "longitude": {"abbreviation": "lon",
                  "full_name": "Longitude",
                  "units": "$^{\circ} E$"},
    "latitude": {"abbreviation": "lat",
                 "full_name": "Latitude",
                 "units": "$^{\circ} N$"},
    "hour": {"abbreviation": "h",
             "full_name": "Hour",
             "units": "UTC {}"},
    "year": {"abbreviation": "y",
             "full_name": "Year",
             "units": "CE"}
}


# In[ ]:


## Supplementary functions for calculations

def create_log_if_directly_executed(time_exec_1up, func_1up=None, func_2up=None, 
                                    args_1up=None, args_1up_values=None):
    
    """
    Creates a debugging log for the function which was directly executed.
    Designed to be called within the start of other functions.
    
    Arguments:
        time_exec_1up (datetime.datetime): Time of when func_1up was executed.
        func_1up (str): Name of function calling this function.
        func_2up (str): Name of function calling the function calling this function.
        args_1up (list): List of argument names for func_1up.
        args_1up_values (dict): Mapping of argument names for func_1up to their
            input values.
    
    Returns:
        ../logs/{func_1up}/{calc_funcs_ver}_{func_1up}({args})_{time_str}.txt:
            Output log file in logs folder. {calc_funcs_ver} is the 
            version of the calc_funcs script being used. {args} is the set
            of arguments input into {func_1up}. {time_str} is a string giving
            the time of when {func_1up} was executed, in "%Y-%m-%d-%H-%M-%S" format.
    """
    
    assert str(type(time_exec_1up)) == "<class 'datetime.datetime'>", \
        "time_exec_1up must be a datetime.datetime object"
    
    if time_exec_1up == None:
        time_exec_1up = datetime.today()
    
    if func_1up == None:
        func_1up = inspect.stack()[1][3]
        
    if func_2up == None:
        func_2up = inspect.stack()[2][3]
        
    if (func_2up == "<cell line: 1>") | (func_2up == "<module>"):
        
        if (args_1up == None) | (args_1up_values == None):
            frame_1up = inspect.currentframe().f_back
            args_1up, _, _, args_1up_values = inspect.getargvalues(frame_1up)
        
        time_str = time_exec_1up.strftime("%Y-%m-%d-%H-%M-%S")
        
        args_1up_list = []
        
        for arg in args_1up:
            arg_value = args_1up_values[arg]
            arg_value_type = str(type(arg_value))
            
            if ((arg_value_type == "<class 'xarray.core.dataset.Dataset'>") | 
                (arg_value_type == "<class 'xarray.core.dataarray.DataArray'>")):
                arg_str = arg
            else:
                arg_str = str(arg_value)
                
            if arg_value_type == "<class 'str'>":
                arg_str = arg_str.replace(arg_value, f"'{arg_value}'")
                
            if arg_value_type == "<class 'function'>":
                arg_str = arg_str.split(" ")[1]
                
            args_1up_list.append(arg_str)
            
        args_1up_str = ", ".join(arg_input for arg_input in args_1up_list)
        path_log = (f"../logs/{func_1up}/({args_1up_str})_" +
                    f"{calc_funcs_ver}_{time_str}")
        Path(f"../logs/{func_1up}").mkdir(parents=True, exist_ok=True)
        # File names can only have maximum 255 characters.
        logging.basicConfig(level=calc_log_level, filename=path_log[:255], force=True)
        
        msg_log = f"CREATED: log file for {func_1up}: {path_log}."
        logging.info(msg_log)
        print(msg_log)


# In[ ]:


def remove_handlers_if_directly_executed(func_2up=None):
    
    """
    Remove all handlers associated with the root logger object if the function was 
    directly executed. Designed to be called within the end of other functions or 
    just before each function return where create_log_if_directly_executed 
    was called within the start of the function. This is used to avoid accidentally 
    appending additional entries to the created log.
    
    Arguments:
        func_2up (str): Name of function calling the function calling this function.
    """
    
    if func_2up == None:
        func_2up = inspect.stack()[2][3]
        
    if (func_2up == "<cell line: 1>") | (func_2up == "<module>"):
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)


# In[ ]:


def terminate_if_file_exists(path_output, func_1up=None, func_2up=None):
    
    """
    Terminates the code if the given path_output already exists. This is to avoid
    overwriting existing files or creating duplicates. This function is designed to 
    be called within one of the calc_func functions.
    
    Arguments:
        path_output (str): Output path for func_1up.
        func_1up (str): Name of function calling this function.
        func_2up (str): Name of function calling the function calling this function.
    """
    
    assert isinstance(path_output, str), \
        "path_output must have data type str"
    
    if func_1up == None:
        func_1up = inspect.stack()[1][3]
        
    if func_2up == None:
        func_2up = inspect.stack()[2][3]
        
    if Path(path_output).exists():
        msg_exist = f"TERMINATED: {func_1up} because file already exists: {path_output}."
        logging.error(msg_exist)
        remove_handlers_if_directly_executed(func_2up)
        raise Exception(msg_exist)


# In[ ]:


def create_output_file(ds, path_output, func_2up=None):
    
    """
    Output the dataset ds as a netcdf4 file into the given path_output.
    
    Arguments:
        ds (xarray.Dataset): Dataset to create output file for.
        path_output (str): Output path for the output file of ds.
        func_2up (str): Name of function calling the function calling this function.
    """
    
    assert str(type(ds)) == "<class 'xarray.core.dataset.Dataset'>", \
        "ds must be an xarray.Dataset"
    assert isinstance(path_output, str), \
        "path_output must have data type str"
    
    if func_2up == None:
        func_2up = inspect.stack()[2][3]
    
    logging.info(f"Creating: file: {path_output}.")
    path_output_dir = "/".join(path_output.split("/")[:-1])
    Path(path_output_dir).mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(path_output)
    
    if (func_2up == "<cell line: 1>") | (func_2up == "<module>"):
        msg_cre_1up = f"CREATED: file: {path_output}."
        logging.info(msg_cre_1up)
        print(msg_cre_1up)
    else:
        msg_cre_2up = f"CREATED: file for use in {func_2up}: {path_output}."
        logging.info(msg_cre_2up)
        print(msg_cre_2up)


# In[ ]:


def check_args_for_none(func_name, args_1up=None, args_1up_values=None):
    
    """
    Function to check whether input arguments are correctly or incorrectly
    set to None. Used to avoid erroneous passes in check_args function.
    
    Arguments:
        func_name (str): Name of function to check Nones for.
        args_1up (list): List of argument names for func_1up.
        args_1up_values (dict): Mapping of argument names for func_1up to their
            input values.
            
    Returns:
        AssertionError if any of the input arguments are incorrectly set to None
            or incorrectly set to something other than None.
    """
    
    if (args_1up == None) | (args_1up_values == None):
        frame_1up = inspect.currentframe().f_back
        args_1up, _, _, args_1up_values = inspect.getargvalues(frame_1up)
        
    if func_name in ["calc_era5_mdp_clim_given_var_or_dvar", 
                     "calc_era5_mdp_clim_stats_given_var_or_dvar",
                     "calc_era5_wsd_clim"]:
        assert args_1up_values["glass_source_pref"] == None, \
            f"glass_source_pref must be None if calc_func = {func_name}"
        args_1up.remove("glass_source_pref")
    
    if func_name in ["calc_glass_mean_clim", "calc_era5_wsd_clim"]:
        assert args_1up_values["var_or_dvar"] == None, \
            f"var_or_dvar must be None if calc_func = {func_name}"
        args_1up.remove("var_or_dvar")
        
    if func_name == "create_orog_static_plot":
        args_1up.remove("region")
    
    # Make exceptions for args in plot_funcs script which can be None.
    for arg_plot in args_plot:
        try:
            args_1up.remove(arg_plot)
        except:
            pass
        
    for arg in args_1up:
        assert args_1up_values[arg] != None, \
            f"{arg} cannot be None"


# In[ ]:


def check_args(
    calc_func=None, region=None, period_start=None, period_end=None, period1_start=None, 
    period1_end=None, period2_start=None, period2_end=None, months_subset=None, 
    glass_source_pref=None, var_or_dvar=None, year_start=None, year_end=None, 
    window_size=None, arg_extra=None, hour=None, time=None, param_orog=None,
    param_glass_mean=None, var_or_dvar_layer=None, var_or_dvar_type=None, perc=None, 
    mask_perc_quantile=None, mask_period1=None, mask_period2=None, extents=None, 
    vmin=None, vmax=None, vmin_periods=None, vmax_periods=None, vmin_diff=None, 
    vmax_diff=None, ax=None, ax_period1=None, ax_period2=None, ax_diff=None, 
    cfv_data=None, output=None
):
    
    """
    Function to check whether input arguments are valid.
    
    Arguments:
        calc_func (function): Calculation function to use in analysis. Must be one of: 
            [calc_glass_mean_clim,
            calc_era5_mdp_clim_given_var_or_dvar,
            calc_era5_mdp_clim_stats_given_var_or_dvar,
            calc_era5_wsd_clim].
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period1_start (str): Start of first period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period1_end (str): End of first period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period2_start (str): Start of second period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period2_end (str): End of second period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        var_or_dvar (str): Variable or value of change in variable to perform
            calculation over. Must be one of: ['u10', 'v10', 'ws10', 'wv10', 'u100', 
            'v100', 'ws100', 'wv100', 'mslp', 't2', 'slhf', 'sshf', 'nse', 'vidmf', 
            'viec', 'vipile', 'vike', 'tcclw', 'tcwv', 'nac', 'blh', 'fa', 'cbh', 'tcc', 
            'cape', 'ci', 'du10', 'dv10', 'dws10', 'dwv10', 'du100', 'dv100', 'dws100', 
            'dwv100', 'dmslp', 'dt2', 'dslhf', 'dsshf', 'dnse', 'dvidmf', 'dviec', 
            'dvipile', 'dvike', 'dtcclw', 'dtcwv', 'dnac', 'dblh', 'dfa', 'dcbh', 
            'dtcc', 'dcape', 'dci'].
        year_start (int): Earliest year to compute the rolling average for.
        year_end (int): Latest year to compute the rolling average for.
        window_size (int): Rolling window size (in years) to compute average for.
            Must be an odd number and greater than or equal to 3.
        arg_extra (str or int): Extra plotting argument used to specify which GLASS 
            parameter to plot, which hour for the mean diurnal profile of an ERA5 
            parameter to plot, which statistic of the mean diurnal profile to plot, 
            or which parameter of the wind speed distribution to plot. Must be one of:
            ["mlai", "mfapar", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
            16, 17, 18, 19, 20, 21, 22, 23, "hour_max", "hour_min", "max", "max_u", 
            "max_v", "min", "min_u", "min_v", "mean", "mean_u", "mean_v", "range", 
            "ws10_mean", "ws10_std", "c10", "k10", "ws100_mean", "ws100_std", "c100", 
            "k100", "eroe100", "tgcf100"].
        hour (int): Hour of mean diurnal profile to plot values for. This is used for
            the plot_funcs script. Must be one of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23].
        time (str): Which times of the day to display MDP values for. This
            is used for the plot_funcs script. Must be one of: 
            ["0-5"/"night", "6-11"/"morning", "12-17"/"afternoon", "18-23"/"evening"].
        param_orog (str): 
        param_glass_mean (str): 
        var_or_dvar_layer (str): Spatial layer from which to draw ERA5 parameters for 
            analysis. This is used for the plot_funcs script. Must be one of: 
            ["sfc", "atm", "cld"].
        var_or_dvar_type (str): Whether to analyse the variables themselves or the 
            change in their mean diurnal profile values as compared with their values
            in the previous hour. This is used for the plot_funcs script.
            Must be one of: ["vars", "dvars"].
        perc (bool): Whether to plot the difference in values as a percentage of the
            value (magnitude if negative) in period1. This is used for the comp plots
            in the plot_funcs script. Must be one of: [True, False].
        mask_perc_quantile (int): If perc is True, specify the quantile of values
            (magnitude if negative) from period1 to mask for the difference plot.
            This is used because percentage differences may be particularly high
            for values which had a low magnitude as a base in period 1.
        mask_period1 (str): Whether to mask grid cells in a comp plot depending on
            whether the value in period 1 was positive or negative. Must be one of:
            ["pos", "neg"].
        mask_period2 (str): Whether to mask grid cells in a comp plot depending on
            whether the value in period 2 was positive or negative. Must be one of:
            ["pos", "neg"].
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes -180
            to 180 and latitudes -90 to 90.
        vmin (float or int): Minimum of colourbar extents for a calc plot.
        vmax (float or int): Maximum of colourbar extents for a calc plot.
        vmin_periods (float or int): Minimum of colourbar extents for a calc plot.
            Used to set common colourbar extents for both periods in a comp plot.
        vmax_periods (float or int): Maximum of colourbar extents for a calc plot.
            Used to set common colourbar extents for both periods in a comp plot.
        vmin_diff (float or int): Minimum of colourbar extents for a diff plot.
            Used to set colourbar extents for the diff plot within a comp plot.
        vmax_diff (float or int): Maximum of colourbar extents for a diff plot.
            Used to set colourbar extents for the diff plot within a comp plot.
        ax (cartopy.GeoAxesSubplot): Figure axis to create plot on.
        ax_period1 (cartopy.GeoAxesSubplot): Figure axis to create calc plot on.
            Used for period 1 plot within a comp plot.
        ax_period2 (cartopy.GeoAxesSubplot): Figure axis to create calc plot on.
            Used for period 2 plot within a comp plot.
        ax_diff (cartopy.GeoAxesSubplot): Figure axis to create diff plot on.
            Used for period2 - period1 diff plot within a comp plot.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns:
        AssertionError if any of the input arguments are invalid.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to check whether input arguments are" + 
                      "valid.")
    else:
        logging.debug(f"Executing: {func_cur} to check whether input arguments into " +
                      f"{func_1up} are valid.")
    
    if calc_func:
        assert callable(calc_func), \
            f"calc_func must be a function and one of: {calc_func_names}"
        calc_func_name = calc_func.__name__
        assert calc_func_name in calc_func_names, \
            f"calc_func must be a function and one of: {calc_func_names}"
        
    if region:
        assert region in [*regions], \
            f"region must be one of: {[*regions]}"
        
    if period_start:
        period_start = datetime.strptime(period_start, "%b-%Y")
        assert period_start >= datetime.strptime(avhrr_earliest, "%b-%Y"), \
            f"period_start must be equal to or later than {avhrr_earliest}"
        
    if period_end:
        period_end = datetime.strptime(period_end, "%b-%Y")
        assert period_end <= datetime.strptime(modis_latest, "%b-%Y"), \
            f"period_end must be equal to or earlier than {modis_latest}"
        
    if (period_start is not None) & (period_end is not None):
        assert period_end >= period_start, \
            "period_end must be equal to or later than period_start"
    
    if period1_start:
        period1_start = datetime.strptime(period1_start, "%b-%Y")
        assert period1_start >= datetime.strptime(avhrr_earliest, "%b-%Y"), \
            f"period1_start must be equal to or later than {avhrr_earliest}"
        
    if period1_end:
        period1_end = datetime.strptime(period1_end, "%b-%Y")
        assert period1_end <= datetime.strptime(modis_latest, "%b-%Y"), \
            f"period1_end must be equal to or earlier than {modis_latest}"
        
    if (period1_start is not None) & (period1_end is not None):
        assert period1_end >= period1_start, \
            "period1_end must be equal to or later than period1_start"
        
    if period2_start:
        period2_start = datetime.strptime(period2_start, "%b-%Y")
        assert period2_start >= datetime.strptime(avhrr_earliest, "%b-%Y"), \
            f"period2_start must be equal to or later than {avhrr_earliest}"
        
    if period2_end:
        period2_end = datetime.strptime(period2_end, "%b-%Y")
        assert period2_end <= datetime.strptime(modis_latest, "%b-%Y"), \
            f"period2_end must be equal to or earlier than {modis_latest}"
        
    if (period2_start is not None) & (period2_end is not None):
        assert period2_end >= period2_start, \
            "period2_end must be equal to or later than period2_start"
    
    if months_subset:
        assert (isinstance(months_subset, list) & (months_subset != []) & 
                all(month in months_subsets["all"] for month in months_subset)
               ) | (months_subset in [*months_subsets]), \
            (f"months_subset must be type str and one of: {[*months_subsets]}, or type " +
             "list and subset of: {} with at least one item".format(months_subsets["all"]))
        
    if (period_start is not None) & (period_end is not None) & (months_subset is not None):
        dates_in_period = pd.date_range(period_start, period_end, freq = "MS")
        months_in_period = set(map(int, dates_in_period.strftime("%-m")))
        if isinstance(months_subset, str):
            months_subset = months_subsets[months_subset]
        assert any(month in months_subset for month in months_in_period), \
            "period must contain at least one month within the given months_subset"
    
    if (period1_start is not None) & (period1_end is not None) & (months_subset is not None):
        dates_in_period1 = pd.date_range(period1_start, period1_end, freq = "MS")
        months_in_period1 = set(map(int, dates_in_period1.strftime("%-m")))
        if isinstance(months_subset, str):
            months_subset = months_subsets[months_subset]
        assert any(month in months_subset for month in months_in_period1), \
            "period1 must contain at least one month within the given months_subset"
        
    if (period2_start is not None) & (period2_end is not None) & (months_subset is not None):
        dates_in_period2 = pd.date_range(period2_start, period2_end, freq = "MS")
        months_in_period2 = set(map(int, dates_in_period2.strftime("%-m")))
        if isinstance(months_subset, str):
            months_subset = months_subsets[months_subset]
        assert any(month in months_subset for month in months_in_period2), \
            "period2 must contain at least one month within the given months_subset"
    
    if glass_source_pref:
        assert glass_source_pref in glass_sources_all, \
            f"glass_source_pref must be one of: {glass_sources_all}"
    
    if var_or_dvar:
        assert var_or_dvar in vars_and_dvars_era5_all, \
            f"var_or_dvar must be one of: {vars_and_dvars_era5_all}"
        
    if year_start:
        year_earliest = int(avhrr_earliest[-4:])
        assert isinstance(year_start, int) & (year_start >= year_earliest), \
            f"year_start must be an integer year equal to or later than {year_earliest}"
        
    if year_end:
        year_latest = int(modis_latest[-4:])
        assert isinstance(year_end, int) & (year_end <= year_latest), \
            f"year_end must be an integer year equal to or earlier than {year_latest}"
        
    if (year_start is not None) & (year_end is not None):
        assert year_end >= year_start, \
            "year_end must be equal to or later than year_start"
        
    if window_size:
        assert (window_size % 2 == 1) & (window_size >= 3), \
            "window_size must be an odd integer greater than or equal to 3"
        
    if (year_start is not None) & (window_size is not None):
        year_earliest_roll = year_earliest + int((window_size-1)/2)
        assert year_start >= year_earliest_roll, \
            (f"year_start must be equal to or later than {year_earliest_roll} " + 
             f"if window_size is {window_size}")
        
    if (year_end is not None) & (window_size is not None):
        year_latest_roll = year_latest - int((window_size-1)/2)
        assert year_end <= year_latest_roll, \
            (f"year_end must be equal to or earlier than {year_latest_roll} " + 
             f"if window_size is {window_size}")
        
    if arg_extra:
        assert arg_extra in arg_extra_all, \
            f"arg_extra must be one of: {arg_extra_all}"
        
    if (arg_extra is not None) & (calc_func is not None):
        if calc_func_name == "calc_glass_mean_clim":
            assert arg_extra in params_glass_mean, \
                (f"arg_extra must be one of: {params_glass_mean} " +
                 f"for calc_func = {calc_func_name}")
        if calc_func_name == "calc_era5_mdp_clim_given_var_or_dvar":
            assert arg_extra in hours_all, \
                (f"arg_extra must be one of: {hours_all} " +
                 f"for calc_func = {calc_func_name}")
            if (arg_extra in ["max_u", "max_v", "min_u", "min_v", 
                              "mean_u", "mean_v"]) & (var_or_dvar is not None):
                assert var_or_dvar in params_vector, \
                    (f"var_or_dvar must be one of: {params_vector} " +
                     f"for arg_extra = {arg_extra}")
        if calc_func_name == "calc_era5_mdp_clim_stats_given_var_or_dvar":
            assert arg_extra in params_stat, \
                (f"arg_extra must be one of: {params_stat} " +
                 f"for calc_func = {calc_func_name}")
        if calc_func_name == "calc_era5_wsd_clim":
            assert arg_extra in params_wsd, \
                (f"arg_extra must be one of: {params_wsd} " +
                 f"for calc_func = {calc_func_name}")
        
    if hour:
        assert hour in hours_all, \
            f"hour must be one of: {hours_all}"
    
    if time:
        assert time in [*times], \
            f"time must be one of: {[*times]}"
    
    if param_orog:
        assert param_orog in params_orog, \
            f"param_orog must be one of {params_orog}"
        
    if param_glass_mean:
        assert param_glass_mean in params_glass_mean, \
            f"param_glass must be one of {params_glass_mean}"
    
    if var_or_dvar_layer:
        assert var_or_dvar_layer in [*[*vars_and_dvars_era5.values()][0].keys()], \
            ("var_or_dvar_layer must be one of: "+
             f"{[*[*vars_and_dvars_era5.values()][0].keys()]}")
        
    if var_or_dvar_type:
        assert var_or_dvar_type in [*vars_and_dvars_era5], \
            f"var_or_dvar_type must be one of: {[*vars_and_dvars_era5]}"
        
    if perc:
        assert perc in [True, False], \
            "perc must be one of: [True, False]"
        
    if mask_perc_quantile:
        assert mask_perc_quantile in range(0, 100+1), \
            "mask_perc_quantile must be an integer between 0 and 100 (inclusive)"
        
    if mask_period1:
        assert mask_period1 in ["pos", "neg"], \
            "mask_period1 must be one of: ['pos', 'neg']"
        
    if mask_period2:
        assert mask_period2 in ["pos", "neg"], \
            "mask_period2 must be one of: ['pos', 'neg']"
        
    if extents:
        assert (isinstance(extents, list) & (len(extents) == 4) & 
                (extents[0] >= -180) & (extents[1] <= 180) &
                (extents[2] >= -90) & (extents[3] <= 90) &
                (extents[1] > extents[0]) & (extents[3] > extents[2])), \
            ("extents must a 4 element list in [W, E, S, N] format " + 
             "with longitudes -180 to 180 and latitudes -90 to 90")
        if region:
            extents_default = regions[region]["extent"]
        else:
            extents_default = [-180, 180, -90, 90]
        assert ((extents[0] >= extents_default[0]) & 
                (extents[1] <= extents_default[1]) & 
                (extents[2] >= extents_default[2]) & 
                (extents[3] <= extents_default[3])), \
            ("extents must be completely contained within " +
             f"{extents_default} for region = {region}")
        
    if vmin:
        assert isinstance(vmin, float) | isinstance(vmin, int), \
            "vmin must have data type float or int"
    
    if vmax:
        assert isinstance(vmax, float) | isinstance(vmax, int), \
            "vmax must have data type float or int"
    
    if (vmin is not None) & (vmax is not None):
        assert vmax >= vmin, \
            "vmax must be equal to or greater than vmin"
        
    if vmin_periods:
        assert isinstance(vmin_periods, float) | isinstance(vmin_periods, int), \
            "vmin_periods must have data type float or int"
    
    if vmax_periods:
        assert isinstance(vmax_periods, float) | isinstance(vmax_periods, int), \
            "vmax_periods must have data type float or int"
    
    if (vmin_periods is not None) & (vmax_periods is not None):
        assert vmax_periods >= vmin_periods, \
            "vmax_periods must be equal to or greater than vmin_periods"
        
    if vmin_diff:
        assert isinstance(vmin_diff, float) | isinstance(vmin_diff, int), \
            "vmin_diff must have data type float or int"
    
    if vmax_diff:
        assert isinstance(vmax_diff, float) | isinstance(vmax_diff, int), \
            "vmax_diff must have data type float or int"
    
    if (vmin_diff is not None) & (vmax_diff is not None):
        assert vmax_diff >= vmin_diff, \
            "vmax_diff must be equal to or greater than vmin_diff"
        
    if ax:
        assert str(type(ax)) == "<class 'cartopy.mpl.geoaxes.GeoAxesSubplot'>", \
            "ax must be a cartopy.GeoAxesSubplot"
        
    if ax_period1:
        assert str(type(ax_period1)) == "<class 'cartopy.mpl.geoaxes.GeoAxesSubplot'>", \
            "ax_period1 must be a cartopy.GeoAxesSubplot"
        
    if ax_period2:
        assert str(type(ax_period2)) == "<class 'cartopy.mpl.geoaxes.GeoAxesSubplot'>", \
            "ax_period2 must be a cartopy.GeoAxesSubplot"
        
    if ax_diff:
        assert str(type(ax_diff)) == "<class 'cartopy.mpl.geoaxes.GeoAxesSubplot'>", \
            "ax_diff must be a cartopy.GeoAxesSubplot"
        
    if cfv_data:
        assert (isinstance(cfv_data, str) & (len(cfv_data) == 5) & 
                (cfv_data[:3] == "cfv") & cfv_data[3].isnumeric() & 
                cfv_data[4].isalpha() & cfv_data[4].islower()) | (cfv_data == "cfv00"), \
            ("cfv_data must be 'cfv00' or of form 'cfvXY' where X is a single digit " +
             "number and Y is a lowercase alphabet character. eg. cfv1n")
    
    if output:
        assert output in [True, False], \
            "output must be one of: [True, False]"
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info("Passed: validity check for input arguments.")
    else:
        logging.info(f"Passed: validity check for input arguments into {func_1up}.")
        
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def select_glass_source(period_start, period_end, glass_source_pref):
    
    """
    Select which GLASS data source (AVHRR or MODIS) to use.
    
    Arguments:
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        
    Returns:
        glass_source (str): String indicating whether to use
            "avhrr" or "modis" data for the given period.
    
    For the given period, select the most appropriate GLASS data source to use 
    (out of AVHRR and MODIS). Where a period is completely contained within
    the time ranges of both AVHRR and MODIS data, glass_source_pref 
    is selected as the data source for use. Otherwise, AVHRR data 
    is used where the given period is completely contained only within the time 
    range of AVHRR data, and conversely for MODIS data. Periods which simultaneously 
    cover both an AVHRR-only period (i.e. before Mar-2000) and a MODIS-only period 
    (i.e. after Dec-2018) are prevented from selection since summary statistics 
    over this range are subject to artefacts from the change in instruments.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(period_start=period_start, period_end=period_end, 
               glass_source_pref=glass_source_pref)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to select the appropriate glass " +
                      f"data source for use between {period_start} and {period_end}.")
    else:
        logging.debug(f"Executing: {func_cur} to select the appropriate glass " + 
                      f"data source for use between {period_start} and " +
                      f"{period_end} in {func_1up}.")
    
    period_start = datetime.strptime(period_start, "%b-%Y")
    period_end = datetime.strptime(period_end, "%b-%Y")
    
    if ((period_start >= datetime.strptime(modis_earliest, "%b-%Y")) &
        (period_end <= datetime.strptime(avhrr_latest, "%b-%Y"))
       ):
        glass_source = glass_source_pref
    elif ((period_start >= datetime.strptime(avhrr_earliest, "%b-%Y")) &
          (period_end <= datetime.strptime(avhrr_latest, "%b-%Y"))
         ):
        glass_source = "avhrr"
    elif ((period_start >= datetime.strptime(modis_earliest, "%b-%Y")) &
          (period_end <= datetime.strptime(modis_latest, "%b-%Y"))
         ):
        glass_source = "modis"
    else:
        raise Exception(f"If period_start is before {modis_earliest}, " +
                        f"period_end cannot be after {avhrr_latest} " +
                        "(since this would cover both an " +
                        "AVHRR-only and a MODIS-only period)")
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info("Selected: {} glass data source for use between {} and {}."
                     .format(glass_source, period_start.strftime("%b-%Y"),
                             period_end.strftime("%b-%Y"))
                    )
    else:
        logging.info("Selected: {} glass data source for use between {} and {} in {}."
                     .format(glass_source, period_start.strftime("%b-%Y"),
                             period_end.strftime("%b-%Y"), func_1up)
                    )
    
    remove_handlers_if_directly_executed(func_1up)    
    return glass_source


# In[ ]:


def get_var_or_dvar_layer_and_type(var_or_dvar):
    
    """
    Obtain the var_or_dvar_layer and var_or_dvar_type classification for the 
    given var_or_dvar.
    
    Arguments:
        var_or_dvar (str): Variable or value of change in variable to perform
            calculation over. Must be one of: ['u10', 'v10', 'ws10', 'wv10', 'u100', 
            'v100', 'ws100', 'wv100', 'mslp', 't2', 'slhf', 'sshf', 'nse', 'vidmf', 
            'viec', 'vipile', 'vike', 'tcclw', 'tcwv', 'nac', 'blh', 'fa', 'cbh', 'tcc', 
            'cape', 'ci', 'du10', 'dv10', 'dws10', 'dwv10', 'du100', 'dv100', 'dws100', 
            'dwv100', 'dmslp', 'dt2', 'dslhf', 'dsshf', 'dnse', 'dvidmf', 'dviec', 
            'dvipile', 'dvike', 'dtcclw', 'dtcwv', 'dnac', 'dblh', 'dfa', 'dcbh', 
            'dtcc', 'dcape', 'dci'].
    
    Returns:
        var_or_dvar_layer (str): var_or_dvar_layer classification indicating which 
            spatial layer var_or_dvar primarily sits in. Can be one of: 
            ["sfc", "atm", "cld"].
        var_or_dvar_type (str): var_or_dvar_type classification indicating whether 
            var_or_dvar is the variable itself ("vars") or the change in the value of a
            variable as compared with the previous hour ("dvars").
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(var_or_dvar=var_or_dvar)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to obtain the var_or_dvar_layer and " +
                      f"var_or_dvar_type classification of {var_or_dvar}.")
    else:
        logging.debug(f"Executing: {func_cur} to obtain the var_or_dvar_layer and " +
                      f"var_or_dvar_type classification of {var_or_dvar} for use in " +
                      f"{func_1up}.")
    
    for var_or_dvar_type_ans, dict_ans in vars_and_dvars_era5.items():
        for var_or_dvar_layer_ans, list_ans in dict_ans.items():
            if var_or_dvar in list_ans:
                var_or_dvar_layer = var_or_dvar_layer_ans
                var_or_dvar_type = var_or_dvar_type_ans
                
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Obtained: {var_or_dvar} classifications of " +
                     f"{var_or_dvar_layer, var_or_dvar_type}.")
    else:
        logging.info(f"Obtained: {var_or_dvar} classifications of " +
                     f"{var_or_dvar_layer, var_or_dvar_type} for use in {func_1up}.")
        
    remove_handlers_if_directly_executed(func_1up)    
    return var_or_dvar_layer, var_or_dvar_type


# In[ ]:


def add_ds_attrs(ds, time_exec_1up=None, func_1up=None, func_2up=None,
                 args_1up=None, args_1up_values=None):
    
    """
    Add attributes to the output dataset from one of the main functions.
    Designed to be called from within one of the main functions.
    
    Arguments:
        ds (xarray.Dataset): Dataset to add attributes for.
        time_exec_1up (datetime.datetime): Time of when func_1up was executed.
        func_1up (str): Name of function calling this function.
            Must be one of: [calc_glass_mean_clim,
            calc_era5_mdp_clim_given_var_or_dvar,
            calc_era5_mdp_clim_stats_given_var_or_dvar,
            calc_era5_wsd_clim,
            calc_diff,
            calc_era5_orog,
            calc_glass_rolling_avg_of_annual_diff].
        func_2up (str): Name of function calling the function calling this function.
        args_1up (list): List of argument names for func_1up.
        args_1up_values (dict): Mapping of argument names for func_1up to their
            input values.
    
    Add attributes for the abbreviation, full name and units for each coordinate in 
    the Dataset. Add attributes for the function executed and time that it was 
    executed to the input dataset itself.
    """
    
    # Obtain arguments if they were not explicitly input.
    
    if time_exec_1up == None:
        time_exec_1up = datetime.today()
        
    if func_1up == None:
        func_1up = inspect.stack()[1][3]
        
    if func_2up == None:
        func_2up = inspect.stack()[2][3]
        
    if (args_1up == None) | (args_1up_values == None):
        frame_1up = inspect.currentframe().f_back
        args_1up, _, _, args_1up_values = inspect.getargvalues(frame_1up)
        
    func_cur = inspect.stack()[0][3]
    
    # The if statement here is redundant since the create_log function also invokes this 
    # same if statement. But by doing this earlier we avoid having to inspect the frame
    # stack an additional time were it not necessary, thus saving some time.
    
    if (func_2up == "<cell line: 1>") | (func_2up == "<module>"):
        frame_cur = inspect.currentframe()
        args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
        create_log_if_directly_executed(time_exec_1up, func_cur, func_1up, 
                                        args_cur, args_cur_values)
    
    # Assert ds is an xarray.Dataset and is the output from one of the main functions.
    
    assert str(type(ds)) == "<class 'xarray.core.dataset.Dataset'>", \
        "ds must be an xarray.Dataset"
    assert func_1up in main_func_names, \
        f"func_1up must be one of: {main_func_names}"
    
    # Logging to indicate this function is being executed.
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to add appropriate attributes to dataset.")
    else:
        logging.debug(f"Executing: {func_cur} to add appropriate attributes to dataset " +
                      f"output from {func_1up}.")
    
    # Add attributes for the coordinates of the input dataset.
    
    for coord in ds.coords:
        ds[coord].attrs = coord_attrs[coord]
        if coord == "hour":
            # Obtain string for timezone relative to UTC (used for hour_max and hour_min).
            region = args_1up_values["region"]
            tz = regions[region]["tz"]
            tz_str = str(tz) if tz < 0 else "+" + str(tz)
            ds[coord].attrs["units"] = ds[coord].attrs["units"].format(tz_str)
            
    # Add attributes for input dataset indicating the function executed and time executed.
    
    time_str = time_exec_1up.strftime("%Y-%m-%d-%H-%M-%S")
    
    args_1up_list = []
        
    for arg in args_1up:
        arg_value = args_1up_values[arg]
        arg_value_type = str(type(arg_value))
        
        if ((arg_value_type == "<class 'xarray.core.dataset.Dataset'>") | 
            (arg_value_type == "<class 'xarray.core.dataarray.DataArray'>")):
            arg_str = arg
        else:
            arg_str = str(arg_value)
            
        if arg_value_type == "<class 'str'>":
            arg_str = arg_str.replace(arg_value, f"'{arg_value}'")
              
        if arg_value_type == "<class 'function'>":
            arg_str = arg_str.split(" ")[1]
                
        args_1up_list.append(arg_str)
            
    args_1up_str = ", ".join(arg_input for arg_input in args_1up_list)
        
    ds.attrs = {"func_executed": f"{func_1up}({args_1up_str})",
                "time_executed": time_str}
    
    # Logging to indicate this function was executed successfully.
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info("Added: attributes to dataset.")
    else:
        logging.info(f"Added: attributes to dataset output from {func_1up}.")
    
    # Remove handlers so if log was created, entries won't be accidentally appended later.
    
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def regrid_era5(ds):
    
    """
    Regrid ERA5 xarray dataset.
    
    Arguments:
        ds (xarray.Dataset): Dataset containing ERA5 data, loaded in
            with xarray using the netcdf4 engine.
                                
    Returns:
        ds_rg (xarray.Dataset): Dataset with regridded coordinates.
        
    Shifts each latitude coordinate south and longitude coordinate east by
    half a grid cell. This reflects the fact that coordinates in the
    original ERA5 dataset defines the north-western corner of each grid
    cell, whereas xarray plots assuming the coordinates refer to the centre.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    assert str(type(ds)) == "<class 'xarray.core.dataset.Dataset'>", \
        "ds must be an xarray.Dataset"
    
    file_name = ds.encoding["source"]
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} on file to make ERA5 coordinates " +
                      f"consistent with xarray plotting: {file_name}.")
    else:
        logging.debug(f"Executing: {func_cur} on file to make ERA5 coordinates " +
                      f"consistent with xarray plotting for use in {func_1up}: " +
                      f"{file_name}.")
    
    ds_rg = (ds
            .assign_coords({"latitude": ds.latitude - res_era5/2,
                            "longitude": (ds.longitude + 180 + res_era5/2) % 360 - 180})
            # Redundant measure just in case longitudes exceed 180 degrees.
            .sortby("longitude")
            )
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Regridded: ERA5 coordinates in file: {file_name}.")
    else:
        logging.debug("Regridded: ERA5 coordinates in file for use in " +
                      f"{func_1up}: {file_name}.")
    
    remove_handlers_if_directly_executed(func_1up)    
    return ds_rg


# In[ ]:


def get_magnitude(da_x, da_y):
    
    """
    Calculate magnitude of 2D vectors given their components.
    
    Arguments:
        da_x (xarray.DataArray): x-component of vectors [m s-1].
        da_y (xarray.DataArray): y-component of vectors [m s-1].
        
    Returns:
        da_r (xarray.DataArray): Magnitude of vectors [m s-1].
        
    Performs a vectorised computation on two different data arrays
    containing the x and y component of some vectors and returns
    the magnitude of the vectors. Dask is allowed.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    assert str(type(da_x)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_x must be an xarray.DataArray"
    assert str(type(da_y)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_y must be an xarray.DataArray"
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug("Executing: {} to calculate vector magnitudes from {} and {}."
                      .format(func_cur, da_x.name, da_y.name)
                     )
    else:
        logging.debug("Executing: {} to calculate vector magnitudes from {} and {} "
                      .format(func_cur, da_x.name, da_y.name) + f"for use in {func_1up}."
                     )
    
    r = lambda x, y: np.sqrt(x**2 + y**2)
    da_r = xr.apply_ufunc(r, da_x, da_y, dask = "allowed")
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info("Obtained: vector magnitudes from {} and {}."
                     .format(da_x.name, da_y.name)
                    )
    else:
        logging.info("Obtained: vector magnitudes from {} and {} for use in {}."
                     .format(da_x.name, da_y.name, func_1up)
                    )
        
    remove_handlers_if_directly_executed(func_1up)
    return da_r


# In[ ]:


def get_da_mdp_for_and_lag(da_mdp):
    
    """
    Obtain DataArray's with values shifted forward (for) or back (lag) 1 hour from input.
    
    Arguments:
        da_mdp (xarray.DataArray): Mean diurnal profile of a variable.
        
    Returns:
        da_for (xarray.DataArray): DataArray with values shifted forward 1 hour.
        da_lag (xarray.DataArray): DataArray with values shifted back 1 hour.
        
    Obtain da_for and da_lag by reassigning hour coordinates and exploiting the fact
    that fixing the mean diurnal profile (MDP) values then moving the corresponding 
    hour back by one position is equivalent to fixing the hour values then moving the
    corresponding MDP values forward by one position (and conversely with moving hours
    forward / values back). The coordinates are then sorted for subtraction or addition,
    for use in handling the fact that some ERA5 variables are instantaneous while others
    are accumulated with hour representing end of accumulation period.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    assert str(type(da_mdp)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_mdp must be an xarray.DataArray"
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to obtain forwarded and lagged " +
                      f"{da_mdp.name} DataArray's.")
    else:
        logging.debug(f"Executing: {func_cur} to obtain forwarded and lagged " +
                      f"{da_mdp.name} DataArray's for use in {func_1up}.")
        
    da_mdp_for = (da_mdp
                  .assign_coords({"hour": (da_mdp.hour - 1) % 24})
                  .sortby("hour")
                 )
    da_mdp_lag = (da_mdp
                  .assign_coords({"hour": (da_mdp.hour + 1) % 24})
                  .sortby("hour")
                 )
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Obtained: forwarded and lagged {da_mdp.name} DataArray's.")
    else:
        logging.info(f"Obtained: forwarded and lagged {da_mdp.name} DataArray's " +
                     f"for use in {func_1up}.")
    
    remove_handlers_if_directly_executed(func_1up)
    return da_mdp_for, da_mdp_lag


# In[ ]:


def get_nac(da_nse, da_vidmf, da_tcwv):
    
    """
    Calculate net atmospheric condensation from ERA5 atmospheric variables.
    
    Arguments:
        da_nse (xarray.DataArray): Net surface evaporation [kg m-2 s-1].
        da_vidmf (xarray.DataArray): Vertical integral of
            divergence of moisture flux [kg m-2 s-1].
        da_tcwv (xarray.DataArray): Total column water vapour [kg m-2].
        
    Returns:
        da_nac (xarray.DataArray): Net atmospheric condensation [kg m-2 s-1].
        
    Performs a vectorised computation on data arrays containing ERA5 atmospheric
    variables and returns the net atmospheric condensation. Dask is allowed.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    assert str(type(da_nse)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_nse must be an xarray.DataArray"
    assert str(type(da_vidmf)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_vidmf must be an xarray.DataArray"
    assert str(type(da_tcwv)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_tcwv must be an xarray.DataArray"
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to calculate net atmospheric condensation " +
                      "from ERA5 atmospheric variables.")
    else:
        logging.debug(f"Executing: {func_cur} to calculate net atmospheric condensation " +
                      f"from ERA5 atmospheric variables for use in {func_1up}.")
            
    logging.debug("Computing: average rate of change in tcwv over consecutive hours to " +
                  "estimate instantaneous rate of change at each hour for use in " +
                  f"{func_cur}.")
    da_tcwv_for, da_tcwv_lag = get_da_mdp_for_and_lag(da_tcwv)
    # Estimate midpoint of rate of change using mean rate (result is in kg m-2 s-1).
    da_tcwv_change = (da_tcwv_for - da_tcwv_lag)/(2*3600)
    
    def nac(nse, vidmf, tcwv_change):
        return nse - vidmf - tcwv_change
    da_nse = xr.apply_ufunc(nac, da_nse, da_vidmf, da_tcwv_change, dask = "allowed")
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info("Obtained: net atmospheric condensation from ERA5 atmospheric " +
                     "variables.")
    else:
        logging.info("Obtained: net atmospheric condensation from ERA5 atmospheric " +
                     f"variables for use in {func_1up}.")
    
    remove_handlers_if_directly_executed(func_1up)
    return da_nse


# In[ ]:


def get_da_range_for_vector_mdp_values(ds_era5_mdp, var_or_dvar):
    
    """
    Calculate the mean diurnal profile (MDP) climatology statistics for a particular
    variable or change in variable (compared to previous hour) using ERA5 data.
    
    Arguments:
        ds_era5_mdp (xarray.Dataset): xarray Dataset containing MDP values for
            the given var_or_dvar.
        var_or_dvar (str): Vector variable or value of change in vector variable to 
            perform calculation over. Must be one of: 
            ['wv10', 'wv100', 'dwv10', 'dwv100'].
                        
    Returns:
        da_range (xarray.DataArray): Range of vector MDP values.
    
    Range for vector values here refers to the largest magnitude of all possible 
    subtractions (corresponding to all possible hourly combinations) between the 24 
    vector MDP values. The calculation creates a dataset for each possible subtraction, 
    merges them into a single dataset along the dimension "iteration", computes the
    magnitude of the vectors from its components, then finds the maximum magnitude
    along the "iteration" dimension.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    assert str(type(ds_era5_mdp)) == "<class 'xarray.core.dataset.Dataset'>", \
        "ds_era5_mdp must be an xarray.Dataset"
    assert var_or_dvar in params_vector, \
        f"var_or_dvar must be one of: {params_vector}"
    assert ((var_or_dvar.replace("wv", "u") in [*ds_era5_mdp.keys()]) & 
            (var_or_dvar.replace("wv", "v") in [*ds_era5_mdp.keys()])
           ), \
        f"ds_era5_mdp must be the MDP dataset for {var_or_dvar}"
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to compute {var_or_dvar} MDP range " +
                      "DataArray.")
    else:
        logging.debug(f"Executing: {func_cur} to compute {var_or_dvar} MDP range " +
                      f"DataArray for use in {func_1up}.")
    
    if priority == "speed":
        # Append every possible subtraction to a list, concatenate into a dataset by 
        # the "iteration" dimension, compute magnitudes for each iteration, then find 
        # maximum magnitude along the iteration dimension.
        ds_subtract_list = []
        for i in range(0, 23+1):
            for j in range(0, 23+1):
                ds_subtract_ij = (ds_era5_mdp.isel(hour = i, drop = True) - 
                                  ds_era5_mdp.isel(hour = j, drop = True))
                ds_subtract_list.append(ds_subtract_ij)
        ds_subtract = xr.concat(ds_subtract_list, "iteration")
        da_range = (get_magnitude(ds_subtract[var_or_dvar.replace("wv", "u")],
                                  ds_subtract[var_or_dvar.replace("wv", "v")])
                    .max("iteration")
                   )
        da_range.name = "range"
        
    elif priority == "memory":
        # For each possible subtraction, first compute the magnitude of that iteration
        # and append to a list, then only update the result corresponding to a
        # (latitude, longitude) coordinate if a following iteration has a larger
        # magnitude at that (latitude, longitude) coordinate.
        da_sub_mag_list = []
        for i in range(0, 23+1):
            for j in range(0, 23+1):
                ds_sub_ij = (ds_era5_mdp.isel(hour = i, drop = True) - 
                             ds_era5_mdp.isel(hour = j, drop = True))
                da_sub_mag_ij = get_magnitude(ds_sub_ij[var_or_dvar.replace("wv", "u")], 
                                              ds_sub_ij[var_or_dvar.replace("wv", "v")])
                da_sub_mag_list.append(da_sub_mag_ij)
                da_sub_mag_list = [xr.concat(da_sub_mag_list, "iteration")
                                   .max("iteration")]
        da_range = xr.DataArray(da_sub_mag_list[0], name = "range")
        
    else:
        raise Exception("priority (global variable in settings section) must be one of: " +
                        "['speed', 'memory'].")
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Obtained: {var_or_dvar} MDP range DataArray.")
    else:
        logging.info(f"Obtained: {var_or_dvar} MDP range DataArray " +
                     f"for use in {func_1up}.")
    
    remove_handlers_if_directly_executed(func_1up)
    return da_range


# In[ ]:


def get_weibull_params(da_mean, da_std):
    
    """
    Obtain Weibull parameters for wind speed distribution from
    mean and standard deviation.
    
    Arguments:
        da_mean (xarray.DataArray): Mean wind speed [m s-1].
        da_std (xarray.DataArray): Standard deviation of wind speed [m s-1].
    
    Returns:
        da_c (xarray.DataArray): Scale parameter for empirical Weibull fit [m s-1].
        da_k (xarray.DataArray): Shape parameter for empirical Weibull fit.
    
    Performs a vectorised computation on two different data arrays
    containing the mean and standard deviation of wind speed and returns
    the Weibull scale and shape parameters for the fit. Dask is allowed.
    This method uses equations (15) and (16) from an article by Justus et al.
    (1977) titled "Methods for Estimating Wind Speed Frequency Distributions".
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    assert str(type(da_mean)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_mean must be an xarray.DataArray"
    assert str(type(da_std)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_std must be an xarray.DataArray"
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug("Executing: {} to obtain Weibull parameters from {} and {}."
                      .format(func_cur, da_mean.name, da_std.name)
                     )
    else:
        logging.debug("Executing: {} to obtain Weibull parameters from {} and {} "
                      .format(func_cur, da_mean.name, da_std.name) +
                      f"for use in {func_1up}."
                     )
    
    k = lambda mean, std: (std / mean)**(-1.086)
    da_k = xr.apply_ufunc(k, da_mean, da_std, dask = "allowed")
    c = lambda mean, k: mean / gamma(1 + 1/k)
    da_c = xr.apply_ufunc(c, da_mean, da_k, dask = "allowed")
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info("Obtained: Weibull parameters from {} and {}."
                     .format(da_mean.name, da_std.name)
                    )
    else:
        logging.info("Obtained: Weibull parameters from {} and {} for use in {}."
                     .format(da_mean.name, da_std.name, func_1up)
                    )
    
    remove_handlers_if_directly_executed(func_1up)
    return da_c, da_k


# In[ ]:


def get_weibull_eroe(da_c, da_k, ws_exc):
    
    """
    Obtain the expected rate of exceedance for a particular wind speed
    from fitted Weibull parameters.
    
    Arguments:
        da_c (xarray.DataArray): Scale parameter for empirical Weibull fit [m s-1].
        da_k (xarray.DataArray): Shape parameter for empirical Weibull fit.
        ws_exc (float or int): Particular wind speed on which to conduct the
            expected rate of exceedance analysis [m s-1].
        
    Returns:
        da_eroe (xarray.DataArray): Expected rate of exceedance from Weibull.
        
    For the given Weibull parameters, the expected rate of exceedance for a 
    particular wind speed is computed as 1 minus the cumulative probability
    distribution for the Weibull fit.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    assert str(type(da_c)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_c must be an xarray.DataArray"
    assert str(type(da_k)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_k must be an xarray.DataArray"
    assert (isinstance(ws_exc, float) | isinstance(ws_exc, int)), \
        "ws_exc must have data type float or int"
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(("Executing: {} to obtain {} m/s expected rate of exceedance " +
                       "from {} and {}.").format(func_cur, ws_exc, da_c.name, da_k.name)
                     )
    else:
        logging.debug(("Executing: {} to obtain {} m/s expected rate of exceedance " +
                       "from {} and {} for use in {}.")
                      .format(func_cur, ws_exc, da_c.name, da_k.name, func_1up)
                     )
    
    eroe = lambda c, k: np.exp(-(ws_exc / c)**k)
    da_eroe = xr.apply_ufunc(eroe, da_c, da_k, dask = "allowed")
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info("Obtained: {} m/s expected rate of exceedance from {} and {}."
                     .format(ws_exc, da_c.name, da_k.name)
                    )
    else:
        logging.info(("Obtained: {} m/s expected rate of exceedance from {} and {} " +
                      "for use in {}.").format(ws_exc, da_c.name, da_k.name, func_1up)
                    )
    
    remove_handlers_if_directly_executed(func_1up)
    return da_eroe


# In[ ]:


def get_gcf(da_ws, speeds, powers, power_max):
    
    """
    Compute the gross capacity factor for a typical wind turbine.
    
    Arguments:
        da_ws (xarray.DataArray): Wind speed data over a period [m s-1].
        speeds (numpy.ndarray): Speed bins for the turbine's power curve [m s-1].
        powers (numpy.ndarray): Powers for each speed bin according to manufacturer
            power curve data [kW].
        power_max (float or int): The maximum power which the turbine can produce [kW].
        
    Returns:
        da_gcf (xarray.DataArray): Gross capacity factor over the period.
        
    First uses the speeds and powers arguments to produce an interpolation
    function for the power curve. Then this interpolation function is applied
    to obtain the power at each wind speed data point in da_ws. Finally, this
    is divided over the maximum power and averaged over the period to obtain
    the gross capacity factor.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    assert str(type(da_ws)) == "<class 'xarray.core.dataarray.DataArray'>", \
        "da_ws must be an xarray.DataArray"
    assert speeds.ndim == 1, \
        (f"speeds argument must be a 1D numpy array with " +
         "data type 'float64' or 'int64'")
    assert (speeds.dtype == "float64") | (speeds.dtype == "int64"), \
        f"speeds must be a 1D numpy array with data type 'float64' or 'int64'"
    assert powers.ndim == 1, \
        f"powers must be a 1D numpy array with data type 'float64' or 'int64'"
    assert (powers.dtype == "float64") | (powers.dtype == "int64"), \
        f"powers must be a 1D numpy array with data type 'float64' or 'int64'"
    assert isinstance(power_max, float) | isinstance(power_max, int), \
        f"power_max must have data type 'float' or 'int'"
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug("Executing: {} to obtain {} gross capacity factor."
                      .format(func_cur, da_ws.name)
                     )
    else:
        logging.debug("Executing: {} to obtain {} gross capacity factor for use in {}."
                      .format(func_cur, da_ws.name, func_1up)
                     )
    
    power_curve = interp1d(speeds, powers, kind = "nearest")
    gcf_instant = lambda ws: power_curve(ws) / power_max
    da_gcf = (xr.apply_ufunc(gcf_instant, da_ws, dask = "allowed")
              .mean("time")
             )
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info("Obtained: {} gross capacity factor.".format(da_ws.name))
    else:
        logging.info("Obtained: {} gross capacity factor for use in {}."
                     .format(da_ws.name, func_1up)
                    )              
    
    remove_handlers_if_directly_executed(func_1up)
    return da_gcf


# In[ ]:


def get_months_subset_str(months_subset):
    
    """
    Obtain the string representation of months_subset for use in output paths.
    
    Arguments:
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
            
    Returns:
        months_subset_str (str): String representation of months_subset.
        
    If months_subset is one of the string inputs then this function does nothing. But
    if months_subset is a list input, this function concatenates the elements of the
    list into a string joined by "-" and in ascending numerical order.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to obtain string representation of " +
                      f"{months_subset}.")
    else:
        logging.debug(f"Executing: {func_cur} to obtain string representation of " +
                      f"{months_subset} for use in {func_1up}.")
    
    if isinstance(months_subset, str):
        months_subset_str = months_subset
    
    if isinstance(months_subset, list):
        months_subset.sort()
        months_subset_str = "-".join(str(month) for month in months_subset)
        
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Obtained: string representation of {months_subset}.")
    else:
        logging.info(f"Obtained: string representation of {months_subset} " +
                     f"for use in {func_1up}.")
    
    remove_handlers_if_directly_executed(func_1up)
    return months_subset_str


# In[ ]:


def convert_period_data_types(period_start, period_end, months_subset):
    
    """
    Convert period_start, period_end and months_subset to appropriate data
    types for use within calculation functions.
    
    Arguments:
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
            
    Returns:
        period_start (datetime.datetime): Start of period to perform calculation over.
        period_end (datetime.datetime): End of period to perform calculation over.
        months_subset (list): Subset of period to perform calculation over.
        
    Converts period_start to period_end to datetime.datetime objects. If months_subset
    is a list then there is no conversion but if it is a string specifying "all" or a 
    season then it will be converted to a list with the months in that subset.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(period_start=period_start, period_end=period_end, 
               months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to convert period data types.")
    else:
        logging.debug(f"Executing: {func_cur} to convert period data types " +
                      f"for use in {func_1up}.")
    
    period_start = datetime.strptime(period_start, "%b-%Y")
    period_end = datetime.strptime(period_end, "%b-%Y")
    if isinstance(months_subset, str):
        months_subset = months_subsets[months_subset]
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Obtained: converted period data types.")
    else:
        logging.info(f"Obtained: converted period data types for use in {func_1up}.")
    
    remove_handlers_if_directly_executed(func_1up)
    return period_start, period_end, months_subset


# In[ ]:


def get_path_for_calc_func(calc_func_name, region, period_start, period_end, 
                           months_subset, glass_source_pref=None,
                           var_or_dvar=None):
    
    """
    Obtain output path for calc_func function.
    
    Arguments:
        calc_func_name (str): Name of calculation function to obtain path for.
            Must be one of: ["calc_glass_mean_clim",
            "calc_era5_mdp_clim_given_var_or_dvar",
            "calc_era5_mdp_clim_stats_given_var_or_dvar",
            "calc_era5_wsd_clim"].
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        var_or_dvar (str): Variable or value of change in variable to perform
            calculation over. Must be one of: ['u10', 'v10', 'ws10', 'wv10', 'u100', 
            'v100', 'ws100', 'wv100', 'mslp', 't2', 'slhf', 'sshf', 'nse', 'vidmf', 
            'viec', 'vipile', 'vike', 'tcclw', 'tcwv', 'nac', 'blh', 'fa', 'cbh', 'tcc', 
            'cape', 'ci', 'du10', 'dv10', 'dws10', 'dwv10', 'du100', 'dv100', 'dws100', 
            'dwv100', 'dmslp', 'dt2', 'dslhf', 'dsshf', 'dnse', 'dvidmf', 'dviec', 
            'dvipile', 'dvike', 'dtcclw', 'dtcwv', 'dnac', 'dblh', 'dfa', 'dcbh', 
            'dtcc', 'dcape', 'dci'].
            
    Returns:
        path_output_calc_func (str): Output path for results from calc_func.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    # Assert that input arguments are valid, and create path for months_subset.
    
    assert calc_func_name in calc_func_names, \
        f"calc_func_name must be one of: {calc_func_names}"
    check_args_for_none(calc_func_name, args_cur, args_cur_values)
    check_args(region=region, period_start=period_start, period_end=period_end, 
               months_subset=months_subset, glass_source_pref=glass_source_pref,
               var_or_dvar=var_or_dvar)
    
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to obtain {calc_func_name} output path.")
    else:
        logging.debug(f"Executing: {func_cur} to obtain {calc_func_name} output path " +
                      f"for use in {func_1up}.")

    # Define path stem.
    
    path_output_calc_func = (f"../data_processed/{calc_func_name[5:None]}/" +
                             f"{calc_funcs_ver}_calc_{region}_{period_start}_" + 
                             f"{period_end}_{months_subset_str}_")
    
    # Append path endings.
    
    if calc_func_name == "calc_glass_mean_clim":
        glass_source = select_glass_source(period_start=period_start, 
                                           period_end=period_end, 
                                           glass_source_pref=glass_source_pref)
        path_output_calc_func += f"glass-mean_{glass_source}.nc"
        
    if calc_func_name == "calc_era5_mdp_clim_given_var_or_dvar":
        path_output_calc_func += f"era5-mdp_{var_or_dvar}.nc"
        
    if calc_func_name == "calc_era5_mdp_clim_stats_given_var_or_dvar":
        path_output_calc_func += f"era5-mdp_{var_or_dvar}_stats.nc"
        
    if calc_func_name == "calc_era5_wsd_clim":
        path_output_calc_func += f"era5-wsd.nc"
        
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Obtained: {calc_func_name} output path.")
    else:
        logging.info(f"Obtained: {calc_func_name} output path for use in {func_1up}.")
    
    remove_handlers_if_directly_executed(func_1up)
    return path_output_calc_func


# In[ ]:


def get_path_for_calc_diff(calc_func_name, region, period1_start, period1_end,
                           period2_start, period2_end, months_subset,
                           glass_source_pref=None, var_or_dvar=None):
    
    """
    Obtain output path for calc_diff function.
    
    Arguments:
        calc_func_name (str): Name of calculation function to compute difference in
            results for. Must be one of: ["calc_glass_mean_clim",
            "calc_era5_mdp_clim_given_var_or_dvar",
            "calc_era5_mdp_clim_stats_given_var_or_dvar",
            "calc_era5_wsd_clim"].
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period1_start (str): Start of first period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period1_end (str): End of first period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period2_start (str): Start of second period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period2_end (str): End of second period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        var_or_dvar (str): Variable or value of change in variable to perform
            calculation over. Must be one of: ['u10', 'v10', 'ws10', 'wv10', 'u100', 
            'v100', 'ws100', 'wv100', 'mslp', 't2', 'slhf', 'sshf', 'nse', 'vidmf', 
            'viec', 'vipile', 'vike', 'tcclw', 'tcwv', 'nac', 'blh', 'fa', 'cbh', 'tcc', 
            'cape', 'ci', 'du10', 'dv10', 'dws10', 'dwv10', 'du100', 'dv100', 'dws100', 
            'dwv100', 'dmslp', 'dt2', 'dslhf', 'dsshf', 'dnse', 'dvidmf', 'dviec', 
            'dvipile', 'dvike', 'dtcclw', 'dtcwv', 'dnac', 'dblh', 'dfa', 'dcbh', 
            'dtcc', 'dcape', 'dci'].
            
    Returns:
        path_output_calc_diff (str): Output path for results from calc_diff.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    # Assert that input arguments are valid, and create path for months_subset.
    
    assert calc_func_name in calc_func_names, \
        f"calc_func_name must be one of: {calc_func_names}"
    check_args_for_none(calc_func_name, args_cur, args_cur_values)
    check_args(region=region, period1_start=period1_start, period1_end=period1_end,
               period2_start=period2_start, period2_end=period2_end,
               months_subset=months_subset, var_or_dvar=var_or_dvar, 
               glass_source_pref=glass_source_pref)
    
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to obtain calc_diff output path " +
                      f"corresponding to {calc_func_name}.")
    else:
        logging.debug(f"Executing: {func_cur} to obtain calc_diff output path " +
                      f"corresponding to {calc_func_name} for use in {func_1up}.")
    
    # Define path stem.
    
    path_output_calc_diff = (f"../data_processed/{calc_func_name[5:None]}/" + 
                             f"{calc_funcs_ver}_diff_{region}_{period1_start}_" +
                             f"{period1_end}_{period2_start}_{period2_end}_" +
                             f"{months_subset_str}_")
    
    # Append path endings.
    
    if calc_func_name == "calc_glass_mean_clim":
        glass_source_period1 = select_glass_source(period_start=period1_start, 
                                                   period_end=period1_end, 
                                                   glass_source_pref=glass_source_pref)
        glass_source_period2 = select_glass_source(period_start=period2_start, 
                                                   period_end=period2_end, 
                                                   glass_source_pref=glass_source_pref)
        
        if glass_source_period1 == glass_source_period2:
            glass_source = glass_source_period1
        else:
            glass_source = "mixed"
            
        path_output_calc_diff += f"glass-mean_{glass_source}.nc"
        
    if calc_func_name == "calc_era5_mdp_clim_given_var_or_dvar":
        path_output_calc_diff += f"era5-mdp_{var_or_dvar}.nc"
        
    if calc_func_name == "calc_era5_mdp_clim_stats_given_var_or_dvar":
        path_output_calc_diff += f"era5-mdp_{var_or_dvar}_stats.nc"
        
    if calc_func_name == "calc_era5_wsd_clim":
        path_output_calc_diff += f"era5-wsd.nc"
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Obtained: calc_diff output path corresponding to " +
                     f"{calc_func_name}.")
    else:
        logging.info(f"Obtained: calc_diff output path corresponding to " +
                     f"{calc_func_name} for use in {func_1up}.")
    
    remove_handlers_if_directly_executed(func_1up)
    return path_output_calc_diff


# In[ ]:


def get_path_for_era5_orog():
    
    """
    Obtain output path for calc_era5_orog function.
            
    Returns:
        path_output_orog (str): Output path for results from calc_era5_orog.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to obtain calc_era5_orog output path.")
    else:
        logging.debug(f"Executing: {func_cur} to obtain calc_era5_orog output path " +
                      f"for use in {func_1up}.")
        
    # Obtain output path.
    
    path_output_orog = (f"../data_processed/era5_orog/{calc_funcs_ver}_" +
                        "calc_global_era5-orog.nc")
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Obtained: calc_era5_orog output path.")
    else:
        logging.info(f"Obtained: calc_era5_orog output path for use in {func_1up}.")
    
    remove_handlers_if_directly_executed(func_1up)
    return path_output_orog


# In[ ]:


def get_path_for_calc_glass_rolling(region, year_start, year_end, months_subset, 
                                    window_size, glass_source_pref):
    
    """
    Obtain output path for calc_glass_rolling_avg_of_annual_diff function.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of ["ca", "sa", "wa"].
        year_start (int): Earliest year to compute the rolling average for.
        year_end (int): Latest year to compute the rolling average for.
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        window_size (int): Rolling window size (in years) to compute average for.
            Must be an odd number and greater than or equal to 3.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
            
    Returns:
        path_output_glass_roll (str): Output path for results from 
            calc_glass_rolling_avg_of_annual_diff.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    # Assert that input arguments are valid, and create path for months_subset.
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(region=region, year_start=year_start, year_end=year_end, 
               months_subset=months_subset, window_size=window_size, 
               glass_source_pref=glass_source_pref)
    
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.debug(f"Executing: {func_cur} to obtain " 
                      + "calc_glass_rolling_avg_of_annual_diff output path.")
    else:
        logging.debug(f"Executing: {func_cur} to obtain " +
                      "calc_glass_rolling_avg_of_annual_diff output path for" +
                      f"use in {func_1up}.")
        
    # Obtain output path.
    
    path_output_glass_roll = ("../data_processed/glass_rolling_avg_of_annual_diff/" +
                              f"{calc_funcs_ver}_calc_{region}_{year_start}_" +
                              f"{year_end}_{months_subset_str}_{window_size}-year_" +
                              f"glass-rolling-diff_pref-{glass_source_pref}.nc")
        
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Obtained: calc_glass_rolling_avg_of_annual_diff output path.")
    else:
        logging.info("Obtained: calc_glass_rolling_avg_of_annual_diff output path " +
                     f"for use in {func_1up}.")
    
    remove_handlers_if_directly_executed(func_1up)
    return path_output_glass_roll


# In[ ]:


## Main calculation functions

def calc_glass_mean_clim(region, period_start, period_end, months_subset, 
                         glass_source_pref, var_or_dvar=None):
    
    """
    Calculate mean leaf area index (MLAI) and mean fraction of absorbed
    photosynthetically active radiation (MFAPAR) climatology using GLASS data.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        var_or_dvar (None): This argument is not used for this analysis. It is used 
            for applying the calc_diff function over an arbitrary calc_func.
                        
    Returns:
        ../data_processed/glass_mean_clim/{calc_funcs_ver}_calc_{region}_{period_start}_
        {period_end}_{months_subset_str}_glass-mean_{glass_source}.nc:
            Output netcdf4 file in data_processed folder containing both MLAI and
            MFAPAR. {calc_funcs_ver} is the version of the calc_funcs script being
            used. {months_subset_str} is a string representing the list of selected 
            months to use as a subset. {glass_source} is automatically selected 
            between ["avhrr", "modis"] based on the selected period. 
    
    For each grid cell, calculate the mean glass climatology (MLAI and MFAPAR). These
    values are computed over the period from period_start to period_end (inclusive),
    and only using a subset of data within this period (if months_subset not "all" is 
    specified). The calculation uses 8-day satellite HDF data from the data_raw
    folder as input, then outputs the result as a netcdf4 file into the
    data_processed folder.
    
    Where a period is completely contained within the time ranges of both AVHRR and 
    MODIS data, glass_source_pref is selected as the data source for use. Otherwise, 
    AVHRR data is used where the given period is completely contained only within 
    the time range of AVHRR data, and conversely for MODIS data. Periods which 
    simultaneously cover both an AVHRR-only period (i.e. before Mar-2000) and 
    a MODIS-only period (i.e. after Dec-2018) are prevented from selection since 
    summary statistics over this range are subject to artefacts from the change in 
    instruments.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    # Assert that input arguments are valid, select the appropriate data source 
    # (AVHRR or MODIS) to use depending on period, and create path for months_subset.
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(region=region, period_start=period_start, period_end=period_end,
               months_subset=months_subset, glass_source_pref=glass_source_pref)
    
    glass_source = select_glass_source(period_start=period_start, 
                                       period_end=period_end,
                                       glass_source_pref=glass_source_pref)
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Executing: {func_cur} to obtain {months_subset_str} " +
                     f"climatology of MLAI and MFAPAR using {glass_source} data " + 
                     f"between {period_start} and {period_end}.")
    else:
        logging.info(f"Executing: {func_cur} to obtain {months_subset_str} " +
                     f"climatology of MLAI and MFAPAR using {glass_source} data " + 
                     f"between {period_start} and {period_end} for use in {func_1up}.")
    
    # Define the output path, convert period_start and period_end to
    # datetime.datetime objects, and months_subset to list if a str was used as input.
    
    path_output_glass_mean = get_path_for_calc_func(
        calc_func_name=func_cur, region=region, period_start=period_start, 
        period_end=period_end, months_subset=months_subset, 
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    terminate_if_file_exists(path_output_glass_mean, func_cur, func_1up)
    
    period_start, period_end, months_subset = convert_period_data_types(
        period_start=period_start, period_end=period_end, months_subset=months_subset)

    # The two functions below are used with xarray's open_mfdataset for parallel
    # computing using dask. The region and times (period with months_subset) are
    # selected within separate functions and uses different logic as compared with
    # filtering in the ERA5 datasets. This is because each GLASS file contains
    # global data whereas the ERA5 datasets were downloaded for each local region.
    
    def filter_glass_files(file_name):
        # This function is used as a mask in conjunction with the default python
        # filter function later, in order to select out the raw data files within
        # the input period and climatological subset within the original function 
        # arguments (by using dates contained within each file name). This is done 
        # (as opposed to using open_mfdataset then filtering) for scalability reasons 
        # since we may need to persist the data in RAM to speed up certain computations.
        time = file_name[-12:-4]
        time = datetime.strptime(time, "%Y-%j")
        if ((time.month in months_subset) &
            # We add an extra month to period_end here because period_end was
            # specified as a month, and conversion into a datetime object
            # defaults to the first (rather than last) day of that month.
            (period_start <= time < period_end + relativedelta(months=1))
           ):
            return True
        else:
            return False
    
    def preprocess_glass(ds):
        # This function is used for the preprocess argument in open_mfdataset.
        # It uses the dates in each raw data file name to assign a time dimension
        # and coordinate for the corresponding dataset. This then forms the
        # dimension along which the files are combined into a single dataset
        # and rechunked. This function also selects out the input region within
        # the original function arguments before the files are concatenated using
        # open_mfdataset (this is again done for persist scalability).
        file_name = ds.encoding["source"]
        logging.debug(f"Preprocessing: file for use in {func_cur}: {file_name}.")
        time = datetime.strptime(file_name[-12:-4], "%Y-%j")
        ds = (ds
              .expand_dims({"time": [time]})
              # Redundant measure just in case longitudes exceed 180 degrees.
              .assign_coords({"x": (ds.x + 180) % 360 - 180})
              .sortby("x")
              .rename({"x": "longitude", "y": "latitude"})
              .drop_vars("spatial_ref")
              .squeeze("band", drop=True)
              )
        ds = ds.sel(longitude=slice(regions[region]["extent"][0],
                                    regions[region]["extent"][1]),
                    latitude=slice(regions[region]["extent"][3],
                                   regions[region]["extent"][2])
                   )
        return ds
    
    # The following code creates the mean climatology datasets for each GLASS
    # variable, by using the previous functions along with open_mfdataset.
    # An initally empty dataset is iteratively appended then merged so that
    # future scalability is possible in case one wishes to add more GLASS
    # variables to the params_glass_mean global python variable.
    
    datasets = []
    
    for param_glass_mean in params_glass_mean:
        param_glass = param_glass_mean[1:]
        files_glass_all = glob(f"../data_raw/global_glass-{param_glass}-" +
                               f"{glass_source}_8-day/global_glass-" +
                               f"{param_glass}-{glass_source}*")
        
        if len(files_glass_all) != number_of_glass_files[param_glass][glass_source]:
            msg_files = (
                f"WARNING: Expected " +
                f"{number_of_glass_files[param_glass][glass_source]} files in " +
                f"../data_raw/global_glass-{param_glass}-{glass_source}_8-day/ " +
                f"but got {len(files_glass_all)}. This could be because the " +
                "data_download.ipynb notebook was not run properly. Or it could be " +
                "that the number of GLASS files on the server from which the data " +
                "was downloaded has changed. Or it may be that the user has changed " +
                "the period coverage from the original values in the " +
                "data_download.ipynb notebook. Alternatively, the user may have " +
                "changed some files in this folder."
            )
            logging.warning(msg_files)
            print(msg_files)
            
        logging.debug(f"Filtering: {param_glass} files from data_raw folder " +
                      f"for use in {func_cur}.")
        files_glass_filtered = list(filter(filter_glass_files, files_glass_all))
        files_glass_filtered.sort()
        
        # This if statement is to ensure an array full of NaNs is returned for MFAPAR
        # when the input period includes 1981 or 2021. At the time of writing, GLASS
        # FAPAR data is not available for these years.
        
        if (param_glass == "fapar") & (
            (period_start < datetime.strptime(fapar_earliest, "%b-%Y")) |
            (period_end > datetime.strptime(fapar_latest, "%b-%Y"))
        ):
            # This line exploits the fact that the for loop runs in sequence and 
            # will have computed MLAI before it attempts to compute MFAPAR. 
            # Therefore an MLAI array with appropriate coordinates already exists 
            # in the datasets list and can be used to create an array of NaNs.
            ds_mean = (datasets[0]
                       .where(np.isnan(datasets[0]["mlai"]))
                       .rename({"mlai": "mfapar"})
                      )
            msg_avail = ("WARNING: GLASS FAPAR data is not available for dates " +
                         f"on or before {fapar_earliest}, and dates on or after " +
                         f"{fapar_latest}. A data array with NaNs was returned " +
                         "for MFAPAR instead.")
            logging.warning(msg_avail)
            print(msg_avail)
            
        else:
            logging.debug(f"Opening: {param_glass} files from data_raw folder for use " +
                          f"in {func_cur}.")
            ds_mean = (xr.open_mfdataset(files_glass_filtered, engine = "rasterio",
                                         preprocess=preprocess_glass, parallel = True)
                       # Rechunking after open_mfdataset here is actually bad practice
                       # since it requires extra computation, but the chunks argument
                       # for open_mfdataset doesn't seem to work here for some reason.
                       .chunk(chunks = {"time": chunksize})
                      )
            if priority == "speed":
                ds_mean = ds_mean.persist()
            
            logging.debug(f"Computing: {param_glass_mean} values for use in {func_cur}.")
            ds_mean = (ds_mean
                       .mean("time")
                       .rename({"band_data": f"{param_glass_mean}"})
                      )
            
        datasets.append(ds_mean)
        
    logging.debug("Merging: glass mean datasets.")
    ds_glass_mean = xr.merge(datasets)
    
    # Add attributes to each DataArray in Dataset.
    
    logging.info("Adding: attributes for each DataArray within output "+
                 f"Dataset from {func_cur}.")
    for da_name in [*ds_glass_mean.keys()]:
            ds_glass_mean[da_name].attrs = copy.deepcopy(attrs_da[da_name])
            ds_glass_mean[da_name].attrs["source"] = glass_source
            
    # Add attributes to Dataset.
    
    add_ds_attrs(ds_glass_mean, time_exec, func_cur, func_1up, args_cur, args_cur_values)
    
    # Create output file in data_processed folder.
    
    create_output_file(ds_glass_mean, path_output_glass_mean, func_1up)
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def calc_era5_mdp_clim_given_var_or_dvar(region, period_start, period_end, 
                                         months_subset, var_or_dvar, 
                                         glass_source_pref=None):
    
    """
    Calculate the mean diurnal profile (MDP) climatology for a particular
    variable or change in variable (compared to previous hour) using ERA5 data.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        var_or_dvar (str): Variable or value of change in variable to perform
            calculation over. Must be one of: ['u10', 'v10', 'ws10', 'wv10', 'u100', 
            'v100', 'ws100', 'wv100', 'mslp', 't2', 'slhf', 'sshf', 'nse', 'vidmf', 
            'viec', 'vipile', 'vike', 'tcclw', 'tcwv', 'nac', 'blh', 'fa', 'cbh', 'tcc', 
            'cape', 'ci', 'du10', 'dv10', 'dws10', 'dwv10', 'du100', 'dv100', 'dws100', 
            'dwv100', 'dmslp', 'dt2', 'dslhf', 'dsshf', 'dnse', 'dvidmf', 'dviec', 
            'dvipile', 'dvike', 'dtcclw', 'dtcwv', 'dnac', 'dblh', 'dfa', 'dcbh', 
            'dtcc', 'dcape', 'dci'].
        glass_source_pref (None): This argument is not used for this analysis. It is used 
            for applying the calc_diff function over an arbitrary calc_func.
                        
    Returns:
        ../data_processed/era5_mdp_clim_given_var_or_dvar/{calc_funcs_ver}_calc_{region}_
        {period_start}_{period_end}_{months_subset_str}_era5-mdp_{var_or_dvar}.nc:
            Output netcdf4 file in data_processed folder containing the MDP for
            var_or_dvar. {calc_funcs_ver} is the version of the calc_funcs script
            being used. {months_subset_str} is a string representing the list of 
            selected months to use as a subset.
    
    For each grid cell, calculate the MDP for the selected var_or_dvar. The MDP
    values are computed over the period from period_start to period_end (inclusive),
    and only using a subset of data within this period (if months_subset not "all" is
    specified). The calculation uses monthly averaged reanalysis by hour of day netcdf4
    data from the data_raw folder as input, then outputs the result as a netcdf4 file
    into the data_processed folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    # Assert that input arguments are valid, and create path for months_subset.
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(region=region, period_start=period_start, period_end=period_end,
               months_subset=months_subset, var_or_dvar=var_or_dvar)
    
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Executing: {func_cur} to obtain {months_subset_str} climatology " +
                     f"of {var_or_dvar} MDP between {period_start} and {period_end}.")
    else:
        logging.info(f"Executing: {func_cur} to obtain {months_subset_str} climatology " +
                     f"of {var_or_dvar} MDP between {period_start} and {period_end} " +
                     f"for use in {func_1up}.")
    
    # Define the output path, convert period_start and period_end to
    # datetime.datetime objects, and months_subset to list if a str was used as input.
    
    path_output_mdp_clim = get_path_for_calc_func(
        calc_func_name=func_cur, region=region, period_start=period_start, 
        period_end=period_end, months_subset=months_subset, var_or_dvar=var_or_dvar)
    terminate_if_file_exists(path_output_mdp_clim, func_cur, func_1up)
    
    period_start, period_end, months_subset = convert_period_data_types(
        period_start=period_start, period_end=period_end, months_subset=months_subset)
    
    # Obtain the var_or_dvar_layer and var_or_dvar_type classifications for the given 
    # var_or_dvar. This is used later to identify which folder to open files from, 
    # and whether to compute changes since the previous hour.
    
    var_or_dvar_layer, var_or_dvar_type = get_var_or_dvar_layer_and_type(
        var_or_dvar=var_or_dvar)
    if var_or_dvar_type == "vars":
        var = var_or_dvar
    if var_or_dvar_type == "dvars":
        var = var_or_dvar[1:]
    
    # The two functions below are used with xarray's open_mfdataset for parallel
    # computing using dask. Together they select out the relevant files to read
    # and persist in memory only the data which is necessary for the computation.
    
    def filter_era5_month_hour_files(file_name):
        # This function is used as a mask in conjunction with the default python
        # filter function later, in order to select out the raw data files with
        # years within the input period. The following preprocess function also
        # selects out the relevant years (as well as months) but by applying a
        # filter on the list of file names first we can avoid preprocessing a 
        # lot of files and hence save on memory.
        year = int(file_name[-7:-3])
        if period_start.year <= year <= period_end.year:
            return True
        else:
            return False
    
    def preprocess_era5_month_hour(ds):
        # This function is used for the preprocess argument in open_mfdataset.
        # It selects out only the subset months for persist scalability,
        # renames variables and sorts data in time order.
        file_name = ds.encoding["source"]
        logging.debug(f"Preprocessing: file for use in {func_cur}: {file_name}.")
        ds = (regrid_era5(ds=ds)[[*vars_deps_and_rename[var].keys()]]
              .rename(vars_deps_and_rename[var])
              .sel(time = ds.time.dt.month.isin(months_subset))
              # The downloaded ERA5 dataset is not sorted in time order (which
              # is a necessity for open_mfdataset, so we sort first over here.
              .sortby("time")
             )
        return ds
    
    # The following code loads in an existing MDP file for the var corresponding to
    # var_or_dvar, or creates one if it doesn't already exist. The True component of
    # the if statement will only ever be run for dvars since if it were a var then
    # path_output_var == path_output_mdp_clim and the code would have terminated earlier.
    
    path_output_var = path_output_mdp_clim.replace(var_or_dvar, var)
    if Path(path_output_var).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_output_var}."
        logging.info(msg_open)
        print(msg_open)
        ds_era5_mdp = xr.open_dataset(path_output_var, engine = "netcdf4")
    else:      
        # The following code opens the relevant monthy ERA5 files then computes
        # mean over each hour of the day.
        files_era5_month_hour = glob(
            f"../data_raw/{region}_era5-slv-{var_or_dvar_layer}_month-hour/*.nc")
        files_era5_month_hour.sort()
        if len(files_era5_month_hour) != number_of_era5_month_hour_files:
            msg_files = (
                f"WARNING: Expected {number_of_era5_month_hour_files} files in " +
                f"../data_raw/{region}_era5-slv-{var_or_dvar_layer}_month-hour/ but " +
                f"got {len(files_era5_month_hour)}. This could be because the " + 
                "data_download.ipynb notebook was not run properly. Or it could " +
                "be that the user has selected a different number of years to " +
                "retrieve data for in the data_download.ipynb notebook as " +
                "compared with the original analysis. Or it may be that the " +
                "user has changed some files in this folder."
            )
            logging.warning(msg_files)
            print(msg_files)
            
        logging.debug(f"Filtering: ERA5 {var_or_dvar_layer} files from data_raw " +
                      f"folder for use in {func_cur}.")
        files_era5_month_hour_filtered = list(filter(filter_era5_month_hour_files, 
                                                     files_era5_month_hour))
        files_era5_month_hour_filtered.sort()
        
        logging.debug(f"Opening: ERA5 {var_or_dvar_layer} files from data_raw " +
                      f"folder for use in {func_cur}.")
        ds_era5_mdp = (xr.open_mfdataset(files_era5_month_hour_filtered,
                                         preprocess=preprocess_era5_month_hour,
                                         engine = "netcdf4", parallel = True)
                       # We add an extra month to period_end here because period_end was
                       # specified as a month, and conversion into a datetime object
                       # defaults to the first (rather than last) day of that month. The
                       # -1 hr is to avoid selecting first hour of the following month.
                       .sel(time = slice(period_start, period_end +
                                         relativedelta(months=1, hours = -1)))
                       # Rechunking after open_mfdataset here is actually bad practice
                       # since it requires extra computation, but the chunks argument
                       # for open_mfdataset doesn't seem to work here for some reason.
                       .chunk(chunks = {"time": chunksize})
                      )
        if priority == "speed":
            ds_era5_mdp = ds_era5_mdp.persist()
        
        if var == "ws10":
            ds_era5_mdp = (get_magnitude(ds_era5_mdp["u10"], ds_era5_mdp["v10"])
                           .to_dataset(name = "ws10")
                          )
            
        if var == "ws100":
            ds_era5_mdp = (get_magnitude(ds_era5_mdp["u100"], ds_era5_mdp["v100"])
                           .to_dataset(name = "ws100")
                          )
        
        logging.debug(f"Computing: MDPs of {[*ds_era5_mdp.keys()]} for use in {func_cur}.")
        ds_era5_mdp = (ds_era5_mdp
                       # The time coordinates for the downloaded ERA5 dataset do not
                       # have perfect alignment across all variables, and sometimes
                       # the hour components for a particular month are distributed
                       # across different days (01 or 02) in that month. But this
                       # shouldn't affect the results of an hour-wise average.
                       .groupby("time.hour")
                       .mean("time")
                      )
        
        # For slhf, sshf and nse: average values from hour before and hour after.
        # This is because these variables are hourly accumulations ending at each hour.
        # The average then provides an estimate of the instantaneous value at the 
        # midpoint between these two accumulation periods. Afterwards, calculate net 
        # atmospheric condensation from atmospheric variables.
        
        if var == "slhf":
            logging.debug(f"Computing: averages for consecutive {var} accumulation " +
                          "periods to estimate instantaneous hourly values for use " +
                          f"in {func_cur}.")
            da_slhf_for, da_slhf_lag = get_da_mdp_for_and_lag(ds_era5_mdp["slhf"])
            # Estimate midpoint and convert from J m-2 to W m-2.
            ds_era5_mdp["slhf"] = (da_slhf_for + da_slhf_lag)/(2*3600)
            
        if var == "sshf":
            logging.debug(f"Computing: averages for consecutive {var} accumulation " +
                          "periods to estimate instantaneous hourly values for use " +
                          f"in {func_cur}.")
            da_sshf_for, da_sshf_lag = get_da_mdp_for_and_lag(ds_era5_mdp["sshf"])
            # Estimate midpoint and convert from J m-2 to W m-2.
            ds_era5_mdp["sshf"] = (da_sshf_for + da_sshf_lag)/(2*3600)
            
        if var == "nac":
            logging.debug("Computing: averages for consecutive nse accumulation " +
                          f"periods to estimate instantaneous {var} hourly values " +
                          f"for use in {func_cur}.")
            da_nse_for, da_nse_lag = get_da_mdp_for_and_lag(ds_era5_mdp["nse"])
            # Estimate midpoint and convert from m of water equivalent to kg m-2 s-1.
            ds_era5_mdp["nse"] = (da_nse_for + da_nse_lag)/(2*3.6)
            ds_era5_mdp["nac"] = get_nac(ds_era5_mdp["nse"], ds_era5_mdp["vidmf"],  
                                         ds_era5_mdp["tcwv"])
        
        # Add attributes to each DataArray within Dataset.
        
        logging.info("Adding: attributes for each DataArray within output "+
                     f"Dataset from {func_cur}.")
        for da_name in [*ds_era5_mdp.keys()]:
            ds_era5_mdp[da_name].attrs = copy.deepcopy(attrs_da[da_name])
        
        # Add attributes to Dataset.
    
        add_ds_attrs(ds_era5_mdp, time_exec, func_cur, func_1up, 
                     args_cur, args_cur_values)
        
        # Modify attributes in Dataset and output if var_or_dvar_type == "dvars" 
        # (i.e. the MDP for var is being computed as an intermediate output 
        # for the dvar MDP).
        
        if var_or_dvar_type == "dvars":
            ds_era5_mdp.attrs["func_executed"] = (
                ds_era5_mdp.attrs["func_executed"]
                .replace(var_or_dvar, var)
            )
            create_output_file(ds_era5_mdp, path_output_var, func_1up)
    
    # If a dvar was specified for var_or_dvar, calculate the change in the value
    # of the variable as compared with its value in the previous hour.
    
    if var_or_dvar_type == "dvars":
        for da_name in [*ds_era5_mdp.keys()]:
            logging.info(f"Obtaining d{da_name} MDP from {da_name} MDP for " +
                         f"use in {func_cur}.")
            _, da_var_lag = get_da_mdp_for_and_lag(ds_era5_mdp[da_name])
            ds_era5_mdp[da_name] += - da_var_lag
            
            # Add attributes to each DataArray within Dataset.
            
            logging.info("Adding: attributes for each DataArray within output "+
                         f"Dataset from {func_cur}.")
            ds_era5_mdp[da_name].attrs = copy.deepcopy(attrs_da[da_name])
            # Append entries to front of attributes if dvar is given.
            ds_era5_mdp[da_name].attrs["abbreviation"] = (
                "d" + ds_era5_mdp[da_name].attrs["abbreviation"])
            ds_era5_mdp[da_name].attrs["full_name"] = (
                "Hourly Change in " + ds_era5_mdp[da_name].attrs["full_name"])
            
            # Rename DataArray names to include "d" in front.
            
            ds_era5_mdp = ds_era5_mdp.rename({da_name: "d" + da_name})
    
        # Add attributes to Dataset.
    
        add_ds_attrs(ds_era5_mdp, time_exec, func_cur, func_1up, 
                     args_cur, args_cur_values)
    
    # Create output file in data_processed folder.
    
    create_output_file(ds_era5_mdp, path_output_mdp_clim, func_1up)
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def calc_era5_mdp_clim_stats_given_var_or_dvar(region, period_start, period_end, 
                                               months_subset, var_or_dvar, 
                                               glass_source_pref=None):
    """
    Calculate the mean diurnal profile (MDP) climatology statistics for a particular
    variable or change in variable (compared to previous hour) using ERA5 data.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        var_or_dvar (str): Variable or value of change in variable to perform
            calculation over. Must be one of: ['u10', 'v10', 'ws10', 'wv10', 'u100', 
            'v100', 'ws100', 'wv100', 'mslp', 't2', 'slhf', 'sshf', 'nse', 'vidmf', 
            'viec', 'vipile', 'vike', 'tcclw', 'tcwv', 'nac', 'blh', 'fa', 'cbh', 'tcc', 
            'cape', 'ci', 'du10', 'dv10', 'dws10', 'dwv10', 'du100', 'dv100', 'dws100', 
            'dwv100', 'dmslp', 'dt2', 'dslhf', 'dsshf', 'dnse', 'dvidmf', 'dviec', 
            'dvipile', 'dvike', 'dtcclw', 'dtcwv', 'dnac', 'dblh', 'dfa', 'dcbh', 
            'dtcc', 'dcape', 'dci'].
        glass_source_pref (None): This argument is not used for this analysis. It is used 
            for applying the calc_diff function over an arbitrary calc_func.
                        
    Returns:
        ../data_processed/era5_mdp_clim_stats_given_var_or_dvar/{calc_funcs_ver}_calc_
        {region}_{period_start}_{period_end}_{months_subset_str}_era5-mdp_
        {var_or_dvar}_stats.nc:
            Output netcdf4 file in data_processed folder containing the hour of maximum,
            hour of minimum, maximum, minimum, mean and range of the MDP for
            var_or_dvar. {calc_funcs_ver} is the version of the calc_funcs script
            being used. {months_subset_str} is a string representing the list of 
            selected months to use as a subset.
    
    For each grid cell, calculate the max, min, mean, range, hour_max and hour_min of 
    the MDP for the selected var_or_dvar. The MDP values are computed over the period 
    from period_start to period_end (inclusive), and only using a subset of data within 
    this period (if months_subset not "all" is specified). If var_or_dvar is one of:
    ["wv10", "wv100", "dwv10", "dwv100"], then hour_max and hour_min refer to the hours
    for when the magnitude of these vectors are at a maximum or minimum respectively.
    Max and min then refer to the vector quantities at hour_max and hour_min respectively.
    Range refers to largest magnitude of all possible subtractions (corresponding to all 
    possible hourly combinations) between the 24 vector MDP values. The calculation uses
    monthly averaged reanalysis by hour of day netcdf4 data from the data_raw folder as 
    input, then outputs the result as a netcdf4 file into the data_processed folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    # Assert that input arguments are valid, and create path for months_subset.
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(region=region, period_start=period_start, period_end=period_end,
               months_subset=months_subset, var_or_dvar=var_or_dvar)
    
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Executing: {func_cur} to obtain {months_subset_str} " +
                     f"climatology of {var_or_dvar} MDP stats between " +
                     f"{period_start} and {period_end}.")
    else:
        logging.info(f"Executing: {func_cur} to obtain {months_subset_str} " +
                     f"climatology of {var_or_dvar} MDP stats between " +
                     f"{period_start} and {period_end} for use in {func_1up}.")
    
    # Define the output path, and create intermediate output file for MDP of var_or_dvar 
    # if it doesn't yet exist.
    
    path_output_mdp_clim_stats = get_path_for_calc_func(
        calc_func_name=func_cur, region=region, period_start=period_start, 
        period_end=period_end, months_subset=months_subset, var_or_dvar=var_or_dvar)
    terminate_if_file_exists(path_output_mdp_clim_stats, func_cur, func_1up)
    
    path_era5_mdp = get_path_for_calc_func(
        calc_func_name="calc_era5_mdp_clim_given_var_or_dvar", region=region, 
        period_start=period_start, period_end=period_end, months_subset=months_subset, 
        var_or_dvar=var_or_dvar)
    if Path(path_era5_mdp).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_era5_mdp}."
        logging.info(msg_open)
        print(msg_open)
    else:
        calc_era5_mdp_clim_given_var_or_dvar(region, period_start, period_end,
                                             months_subset, var_or_dvar)
        
    # Convert period_start and period_end to datetime.datetime objects, 
    # and months_subset to list if a str was used as input.
    
    period_start, period_end, months_subset = convert_period_data_types(
        period_start=period_start, period_end=period_end, months_subset=months_subset)
    
    # Compute MDP stats for given var_or_dvar, treating vector values separately.
    logging.debug(f"Computing: {var_or_dvar} MDP stats for use in {func_cur}.")    
    if var_or_dvar in params_vector:
        if var_or_dvar in vars_era5_all:
            var_or_dvar_attrs = copy.deepcopy(attrs_da[var_or_dvar])
        if var_or_dvar in dvars_era5_all:
            # Append entries to front of attributes if dvar is given.
            var = var_or_dvar[1:]
            var_or_dvar_attrs = copy.deepcopy(attrs_da[var])
            var_or_dvar_attrs["abbreviation"] = "d" + var_or_dvar_attrs["abbreviation"]
            var_or_dvar_attrs["full_name"] = ("Hourly Change in " +
                                              var_or_dvar_attrs["full_name"])
        ds_era5_mdp = xr.open_dataset(path_era5_mdp, engine = "netcdf4")
        if priority == "speed":
            ds_era5_mdp = ds_era5_mdp.persist()
        da_u = ds_era5_mdp[var_or_dvar.replace("wv", "u")]
        da_v = ds_era5_mdp[var_or_dvar.replace("wv", "v")]
        da_mag = get_magnitude(da_u, da_v)
        da_hour_max = xr.DataArray(da_mag.idxmax("hour"), name = "hour_max")
        da_hour_max = (da_hour_max + regions[region]["tz"]) % 24
        da_hour_min = xr.DataArray(da_mag.idxmin("hour"), name = "hour_min")
        da_hour_min = (da_hour_min + regions[region]["tz"]) % 24
        da_max_u = xr.DataArray(da_u.sel(hour = da_hour_max, drop = True), name = "max_u")
        da_max_v = xr.DataArray(da_v.sel(hour = da_hour_max, drop = True), name = "max_v")
        da_min_u = xr.DataArray(da_u.sel(hour = da_hour_min, drop = True), name = "min_u")
        da_min_v = xr.DataArray(da_v.sel(hour = da_hour_min, drop = True), name = "min_v")
        da_mean_u = xr.DataArray(da_u.mean("hour"), name = "mean_u")
        da_mean_v = xr.DataArray(da_v.mean("hour"), name = "mean_v")
        da_range = get_da_range_for_vector_mdp_values(ds_era5_mdp=ds_era5_mdp, 
                                                      var_or_dvar=var_or_dvar)
        logging.debug(f"Merging: {var_or_dvar} MDP stat DataArray's into Dataset.")
        ds_era5_mdp_stats = xr.merge([da_hour_max, da_hour_min, da_max_u, da_max_v, 
                                      da_min_u, da_min_v, da_mean_u, da_mean_v, da_range])
    else:
        da_era5_mdp = xr.open_dataset(path_era5_mdp, engine = "netcdf4")[var_or_dvar]
        var_or_dvar_attrs = copy.deepcopy(da_era5_mdp.attrs)
        if priority == "speed":
            da_era5_mdp = da_era5_mdp.persist()
        da_hour_max = xr.DataArray(da_era5_mdp.idxmax("hour"), name = "hour_max")
        da_hour_max = (da_hour_max + regions[region]["tz"]) % 24
        da_hour_min = xr.DataArray(da_era5_mdp.idxmin("hour"), name = "hour_min")
        da_hour_min = (da_hour_min + regions[region]["tz"]) % 24
        da_max = xr.DataArray(da_era5_mdp.max("hour"), name = "max")
        da_min = xr.DataArray(da_era5_mdp.min("hour"), name = "min")
        da_mean = xr.DataArray(da_era5_mdp.mean("hour"), name = "mean")
        da_range = xr.DataArray(da_max - da_min, name = "range")
        logging.debug(f"Merging: {var_or_dvar} MDP stat DataArray's into Dataset.")
        ds_era5_mdp_stats = xr.merge([da_hour_max, da_hour_min, da_max, da_min, 
                                      da_mean, da_range])
    
    # Obtain string for timezone relative to UTC (used for hour_max and hour_min).
    
    tz = regions[region]["tz"]
    tz_str = str(tz) if tz < 0 else "+" + str(tz)
    
    # Add attributes to each DataArray within Dataset.
    
    logging.info("Adding: attributes for each DataArray within output "+
                 f"Dataset from {func_cur}.")
    for da_name in [*ds_era5_mdp_stats.keys()]:
        ds_era5_mdp_stats[da_name].attrs = copy.deepcopy(attrs_da[da_name])
        # Append name of each statistic to front of var_or_dvar attributes.
        ds_era5_mdp_stats[da_name].attrs["abbreviation"] = (
            ds_era5_mdp_stats[da_name].attrs["abbreviation"]
            .format(var_or_dvar_attrs["abbreviation"])
        )
        ds_era5_mdp_stats[da_name].attrs["full_name"] = (
            ds_era5_mdp_stats[da_name].attrs["full_name"]
            .format(var_or_dvar_attrs["full_name"])
        )
        if (da_name == "hour_max") | (da_name == "hour_min"):
            ds_era5_mdp_stats[da_name].attrs["units"] = (
                ds_era5_mdp_stats[da_name].attrs["units"]
                .format(tz_str)
            )
        else:
            ds_era5_mdp_stats[da_name].attrs["units"] = (
                ds_era5_mdp_stats[da_name].attrs["units"]
                .format(var_or_dvar_attrs["units"])
            )
            
    # Add attributes to Dataset.
    
    add_ds_attrs(ds_era5_mdp_stats, time_exec, func_cur, func_1up, 
                 args_cur, args_cur_values)
    
    # Create output file in data_processed folder.
    
    create_output_file(ds_era5_mdp_stats, path_output_mdp_clim_stats, func_1up)
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def calc_era5_wsd_clim(region, period_start, period_end, months_subset, 
                       glass_source_pref=None, var_or_dvar=None):

    """
    Calculate climatology of wind speed distribution (WSD) properties using 
    ERA5 data for heights of 10 m and 100 m above surface.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        glass_source_pref (None): This argument is not used for this analysis. It is used 
            for applying the calc_diff function over an arbitrary calc_func.
        var_or_dvar (None): This argument is not used for this analysis. It is used 
            for applying the calc_diff function over an arbitrary calc_func.
                        
    Returns:
        ../data_processed/era5_wsd_clim/{calc_funcs_ver}_calc_{region}_{period_start}_
        {period_end}_{months_subset_str}_era5-wsd.nc:
            Output netcdf4 file in data_processed folder containing Weibull shape and
            scale parameters for wind speed at 10 m and 100 m above surface, expected
            rate of exceedance for a particular wind speed at 100 m above surface,
            and gross capacity factor for a typical turbine at 100 m above surface.
            {calc_funcs_ver} is the version of the calc_funcs script being used.
            {months_subset_str} is a string representing the list of selected months 
            to use as a subset.
    
    For each grid cell, calculate the Weibull scale and shape parameter for wind speed
    at 10 m above surface (C10 and K10), the Weibull scale and shape parameter for
    wind speed at 100 m above surface (C100 and K100), the expected rate of exceedance
    for a particular wind speed at 100 m above surface (EROE100) and the gross
    capacity factor for a typical wind turbine at 100 m above surface (TGCF100). The
    wind speed distributions (WSDs) are computed over the period between period_start
    and period_end (inclusive), and only using a subset of data within this period
    (if a months_subset not "all" is specified). The calculation uses hourly ERA5 
    netcdf4 data from the data_raw folder as input, then outputs the result as a
    netcdf4 file into the data_processed folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    # Assert that input arguments are valid, and create path for months_subset.
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(region=region, period_start=period_start, period_end=period_end,
               months_subset=months_subset)
    
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Executing: {func_cur} to obtain wind speed distribution " +
                     f"parameters over {months_subset_str} months between " + 
                     f"{period_start} and {period_end}.")
    else:
        logging.info(f"Executing: {func_cur} to obtain wind speed distribution " +
                     f"parameters over {months_subset_str} months between " +
                     f"{period_start} and {period_end} for use in {func_1up}.")
    
    # Define the output path, convert period_start and period_end to
    # datetime.datetime objects, and months_subset to list if a str was used as input.
    
    path_output_wsd_clim = get_path_for_calc_func(
        calc_func_name=func_cur, region=region, period_start=period_start, 
        period_end=period_end, months_subset=months_subset, var_or_dvar=var_or_dvar)
    terminate_if_file_exists(path_output_wsd_clim, func_cur, func_1up)
    
    period_start, period_end, months_subset = convert_period_data_types(
        period_start=period_start, period_end=period_end, months_subset=months_subset)
    if period_start + relativedelta(years=5) > period_end:
        msg_years = ("WARNING: It is recommended to use at least 5 years of data " +
                     "for the wind speed distribution analysis.")
        logging.warning(msg_years)
        print(msg_years)
    
    # The two functions below are used with xarray's open_mfdataset for parallel
    # computing using dask. Together they select out the relevant files to read
    # and persist in memory only the data which is necessary for the computation.
    
    def filter_era5_hour_files(file_name):
        # This function is used as a mask in conjunction with the default python
        # filter function later, in order to select out the raw data files with
        # years within the input period. The following preprocess function also
        # selects out the relevant years (as well as months) but by applying a
        # filter on the list of file names first we can avoid preprocessing a 
        # lot of files and hence save on memory.
        year = int(file_name[-7:-3])
        if period_start.year <= year <= period_end.year:
            return True
        else:
            return False
    
    def preprocess_era5_hour(ds):
        # This function is used for the preprocess argument in open_mfdataset.
        # It selects out only the wind speed data, and only for
        # the subset months. This is done for persist scalability.
        file_name = ds.encoding["source"]
        logging.debug(f"Preprocessing: file for use in {func_cur}: {file_name}.")
        ds = (regrid_era5(ds=ds)[["u10", "v10", "u100", "v100"]]
              .sel(time = ds.time.dt.month.isin(months_subset))
             )
        return ds
    
    # The following code opens the relevant hourly ERA5 files then uses the
    # u and v components of wind velocity to compute the magnitude (wind speed).
    
    files_era5_hour = glob(f"../data_raw/{region}_era5-slv-sfc_hour/*.nc")
    files_era5_hour.sort()
    if len(files_era5_hour) != number_of_era5_hour_files:
        msg_files = (
            f"WARNING: Expected {number_of_era5_hour_files} files in " +
            f"../data_raw/{region}_era5-slv-sfc_hour/ but got {len(files_era5_hour)}. " +
            "This could be because the data_download.ipynb notebook was not run " + 
            "properly. Or it could be that the user has selected a different " +
            "number of years to retrieve data for in the data_download.ipynb " +
            "notebook as compared with the original analysis. Or it may be that " +
            "the user has changed some files in this folder."
        )
        logging.warning(msg_files)
        print(msg_files)
        
    logging.debug(f"Filtering: ERA5 atm files from data_raw folder for use in {func_cur}.")
    files_era5_hour_filtered = list(filter(filter_era5_hour_files, files_era5_hour))
    files_era5_hour_filtered.sort()
    
    logging.debug(f"Opening: ERA5 atm files from data_raw folder for use in {func_cur}.")
    ds_era5_hour = (xr.open_mfdataset(files_era5_hour_filtered, engine = "netcdf4",
                                     preprocess=preprocess_era5_hour, parallel = True)
                    # We add an extra month to period_end here because period_end was
                    # specified as a month, and conversion into a datetime object
                    # defaults to the first (rather than last) day of that month. The
                    # -1 hr is to avoid selecting first hour of the following month.
                    .sel(time = slice(period_start, period_end +
                                relativedelta(months=1, hours = -1)))
                    # Rechunking after open_mfdataset here is actually bad practice
                    # since it requires extra computation, but the chunks argument
                    # for open_mfdataset doesn't seem to work here for some reason.
                    .chunk(chunks = {"time": chunksize})
                   )
    if priority == "speed":
        ds_era5_hour = ds_era5_hour.persist()
        
    da_ws10 = get_magnitude(ds_era5_hour["u10"], ds_era5_hour["v10"])
    da_ws100 = get_magnitude(ds_era5_hour["u100"], ds_era5_hour["v100"])
    da_ws10.name, da_ws100.name = "ws10", "ws100"
    
    # Compute the mean and standard deviation of wind speed, from these the 
    # Weibull scale and shape parameters, then the expected rate of exceedance
    # and typical gross capacity factor, then combine into a single dataset.
    
    da_ws10_mean = xr.DataArray(da_ws10.mean("time"), name = "ws10_mean")
    da_ws10_std = xr.DataArray(da_ws10.std("time"), name = "ws10_std")
    da_ws10_c, da_ws10_k = get_weibull_params(da_ws10_mean, da_ws10_std)
    da_ws10_c.name, da_ws10_k.name = "c10", "k10"
    da_ws100_mean = xr.DataArray(da_ws100.mean("time"), name = "ws100_mean")
    da_ws100_std = xr.DataArray(da_ws100.std("time"), name = "ws100_std")
    da_ws100_c, da_ws100_k = get_weibull_params(da_ws100_mean, da_ws100_std)
    da_ws100_c.name, da_ws100_k.name = "c100", "k100"
    da_ws100_eroe = get_weibull_eroe(da_ws100_c, da_ws100_k, speed_eroe)
    da_ws100_gcf = get_gcf(da_ws100, speeds_common, powers_avg, power_nameplate)
    da_ws100_eroe.name, da_ws100_gcf.name = "eroe100", "tgcf100"
    
    logging.debug("Merging: wind speed distribution parameter DataArray's into Dataset.")
    ds_era5_wsd = xr.merge([da_ws10_mean, da_ws10_std, da_ws10_c, da_ws10_k,
                            da_ws100_mean, da_ws100_std, da_ws100_c,
                            da_ws100_k, da_ws100_eroe, da_ws100_gcf])
    
    # Add attributes to each DataArray within Dataset.
    
    logging.info("Adding: attributes for each DataArray within output "+
                 f"Dataset from {func_cur}.")
    for da_name in [*ds_era5_wsd.keys()]:
        ds_era5_wsd[da_name].attrs = copy.deepcopy(attrs_da[da_name])
    
    # Add attributes to Dataset.
    
    add_ds_attrs(ds_era5_wsd, time_exec, func_cur, func_1up, args_cur, args_cur_values)
    
    # Create output file in data_processed folder.
    
    create_output_file(ds_era5_wsd, path_output_wsd_clim, func_1up)
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


## Main extra functions

def calc_diff(calc_func, region, period1_start, period1_end,
              period2_start, period2_end, months_subset, 
              glass_source_pref=None, var_or_dvar=None):
    
    """
    Calculates the difference in results for two separate periods which have
    each been outputted by the same calculation function.
    
    Arguments:
        calc_func (function): Calculation function to use in analysis. Must be one of: 
            [calc_glass_mean_clim,
            calc_era5_mdp_clim_given_var_or_dvar,
            calc_era5_mdp_clim_stats_given_var_or_dvar,
            calc_era5_wsd_clim].
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period1_start (str): Start of first period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period1_end (str): End of first period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period2_start (str): Start of second period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period2_end (str): End of second period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        var_or_dvar (str): Variable or value of change in variable to perform
            calculation over. Must be one of: ['u10', 'v10', 'ws10', 'wv10', 'u100', 
            'v100', 'ws100', 'wv100', 'mslp', 't2', 'slhf', 'sshf', 'nse', 'vidmf', 
            'viec', 'vipile', 'vike', 'tcclw', 'tcwv', 'nac', 'blh', 'fa', 'cbh', 'tcc', 
            'cape', 'ci', 'du10', 'dv10', 'dws10', 'dwv10', 'du100', 'dv100', 'dws100', 
            'dwv100', 'dmslp', 'dt2', 'dslhf', 'dsshf', 'dnse', 'dvidmf', 'dviec', 
            'dvipile', 'dvike', 'dtcclw', 'dtcwv', 'dnac', 'dblh', 'dfa', 'dcbh', 
            'dtcc', 'dcape', 'dci'].
    
    Returns:
        ../data_processed/glass_mean_clim/{calc_funcs_ver}_diff_{region}_{period1_start}_
            {period1_end}_{period2_start}_{period2_end}_{months_subset_str}_glass-mean_
            {glass_source}.nc OR
        ../data_processed/era5_mdp_clim_given_var_or_dvar/{calc_funcs_ver}_diff_{region}_
            {period1_start}_{period1_end}_{period2_start}_{period2_end}_
            {months_subset_str}_era5-mdp_{var_or_dvar}.nc OR
        ../data_processed/era5_mdp_clim_stats_given_var_or_dvar/{calc_funcs_ver}_diff_
            {region}_{period1_start}_{period1_end}_{period2_start}_{period2_end}_
            {months_subset_str}_era5-mdp_{var_or_dvar}_stats.nc OR
        ../data_processed/era5_wsd_clim/{calc_funcs_ver}_diff_{region}_{period1_start}_
            {period1_end}_{period2_start}_{period2_end}_{months_subset_str}_era5-wsd.nc:
                Output netcdf4 file in data_processed folder containing the difference
                in results, with name depending on calc_func being used.
                {calc_funcs_ver} is the version of the calc_funcs script being
                used. {months_subset_str} is a string representing the list of selected 
                months to use as a subset. {glass_source} is automatically selected 
                between ["avhrr", "modis"] based on the selected period. 
    
    First runs calc_func for each of the given periods if this has not already
    been done. Then calculates the difference in results as period2 - period1.
    For hour_max and hour_min stats, the result is expressed as a value between
    -12 (hour_max or hour_min for period2 is 12 hours behind of that for period1) and
    +12 (hour_max or hour_min for period2 is 12 hours ahead of that for period1).
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    # Assert that input arguments are valid.
    
    calc_func_name = calc_func.__name__
    check_args_for_none(calc_func_name, args_cur, args_cur_values)
    check_args(calc_func=calc_func, region=region, period1_start=period1_start,
               period1_end=period1_end, period2_start=period2_start, 
               period2_end=period2_end, months_subset=months_subset, 
               glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Executing: {func_cur} to obtain difference in outputs from " +
                     f"{calc_func_name} over {period1_start} to {period1_end} and " +
                     f"{period2_start} to {period2_end}.")
    else:
        logging.info(f"Executing: {func_cur} to obtain difference in outputs from " +
                     f"{calc_func_name} over {period1_start} to {period1_end} and " +
                     f"{period2_start} to {period2_end} for use in {func_1up}.")
    
    # Obtain paths for calc_diff output as well as intermediate calc_func outputs 
    # from each period.
    
    path_output_diff = get_path_for_calc_diff(
        calc_func_name=calc_func_name, region=region, period1_start=period1_start,
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end,
        months_subset=months_subset, glass_source_pref=glass_source_pref,
        var_or_dvar=var_or_dvar)
    terminate_if_file_exists(path_output_diff, func_cur, func_1up)
    
    # Create intermediate output files from each period if they don't already
    # exist, then read in these files as xarray datasets and compute difference.
    path_period1 = get_path_for_calc_func(
        calc_func_name=calc_func_name, region=region, period_start=period1_start, 
        period_end=period1_end, months_subset=months_subset, 
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    if Path(path_period1).exists():
        msg_open1 = f"Opening: existing file for use in {func_cur}: {path_period1}."
        logging.info(msg_open1)
        print(msg_open1)
    else:
        calc_func(region=region, period_start=period1_start, period_end=period1_end, 
                  months_subset=months_subset, glass_source_pref=glass_source_pref, 
                  var_or_dvar=var_or_dvar)
    ds_period1 = xr.open_dataset(path_period1, engine = "netcdf4")
    
    path_period2 = get_path_for_calc_func(
        calc_func_name=calc_func_name, region=region, period_start=period2_start, 
        period_end=period2_end, months_subset=months_subset,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    if Path(path_period2).exists():
        msg_open2 = f"Opening: existing file for use in {func_cur}: {path_period2}."
        logging.info(msg_open2)
        print(msg_open2)
    else:
        calc_func(region=region, period_start=period2_start, period_end=period2_end,
                  months_subset=months_subset, glass_source_pref=glass_source_pref,
                  var_or_dvar=var_or_dvar)
    ds_period2 = xr.open_dataset(path_period2, engine = "netcdf4")
    
    if priority == "speed":
        ds_period1 = ds_period1.persist()
        ds_period2 = ds_period2.persist()
    
    logging.debug(f"Computing: difference in results for {calc_func} over the " + 
                  "given periods.")
    ds_diff = ds_period2 - ds_period1
    
    # Treat hour_max and hour_min stats separately since the difference in these
    # should be between -12 hours (12 hours behind) and +12 hours (12 hours ahead).
    
    if calc_func_name == "calc_era5_mdp_clim_stats_given_var_or_dvar":
        logging.debug("Converting: difference in hour_max and hour_min values to be " + 
                      "between -12 and +12 hours.")
        ds_diff["hour_max"] = (ds_diff["hour_max"] + 12) % 24 - 12
        ds_diff["hour_min"] = (ds_diff["hour_min"] + 12) % 24 - 12
    
    # Add attributes to each DataArray within Dataset.
    
    logging.info("Adding: attributes for each DataArray within output "+
                 f"Dataset from {func_cur}.")
    for da_name in [*ds_diff.keys()]:
        # The subject over which the difference is being computed for appears right
        # before the opening bracket ( for statistics abbreviations, and at the end
        # of the abbreviation for variables. The following line accounts for both
        # in selecting the subject, then appends a "diff" superscript to the subject.
        ds_diff[da_name].attrs = copy.deepcopy(ds_period1[da_name].attrs)
        subject = ds_diff[da_name].attrs["abbreviation"].split("(")[0]
        ds_diff[da_name].attrs["abbreviation"] = (ds_diff[da_name].attrs["abbreviation"]
                                                  .replace(subject, subject + "$^{diff}$"))
        ds_diff[da_name].attrs["full_name"] = ("Difference in " + 
                                               ds_diff[da_name].attrs["full_name"])
        if (da_name == "hour_max") | (da_name == "hour_min"):
            ds_diff[da_name].attrs["units"] = "$h$"
        if calc_func_name == "calc_glass_mean_clim":
            ds_diff[da_name].attrs["source"] = path_output_diff[-8:-3]
    
    # Add attributes to Dataset.
    
    add_ds_attrs(ds_diff, time_exec, func_cur, func_1up, args_cur, args_cur_values)
    
    # Create output file in data_processed folder.
    
    create_output_file(ds_diff, path_output_diff, func_1up)
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def calc_era5_orog():
    
    """
    Calculate the elevation and slope of sub-gridscale orography for the global
    land surface using ERA5 data.
                        
    Returns:
        ../data_processed/era5_orog_slope/{calc_funcs_ver}_calc_
        global_static_orog.nc:
            Output netcdf4 file in data_processed folder containing the elevation and 
            slope of sub-gridscale orography for the global land surface. 
            {calc_funcs_ver} is the version of the calc_funcs script being used.
    
    For each grid cell, calculate the global land surface elevation (m above mean sea 
    level) using ERA5 geopotential data, obtain the global slope of sub-gridscale 
    orography directly from ERA5 data, apply a mask to both datasets, then change the 
    datasets' attributes. ERA5 land-sea mask data containing the fraction of each grid 
    cell which is land is then used to mask values above sea cover (by selecting only 
    the values where this fraction was greater than 0.5). The calculations use static 
    reanalysis netcdf4 data from the data_raw folder as input, then outputs the result 
    (static data) as a netcdf4 file into the data_processed folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Executing: {func_cur} to obtain global orographic variables.")
    else:
        logging.info(f"Executing: {func_cur} to obtain global orographic variables " +
                     f"for use in {func_1up}.")
    
    # Create paths, open raw datasets, calculate the land surface elevation,
    # then mask values and change dataset attributes
    
    path_output_orog = get_path_for_era5_orog()
    terminate_if_file_exists(path_output_orog, func_cur, func_1up)
    
    path_static = "../data_raw/global_era5-slv_static/global_era5-slv_static_all.nc"
    if Path(path_static).exists():
        msg_exist = f"Opening: existing file for use in {func_cur}: {path_static}."
        logging.info(msg_exist)
        print(msg_exist)
        ds_static = regrid_era5(xr.open_dataset(path_static, engine = "netcdf4"))
    else:
        msg_miss = (f"TERMINATED: file could not be found: {path_static}. This could " +
                    "be because the data_download.ipynb notebook was not run properly. " +
                    "Alternatively, the user may have changed some files in this folder.")
        logging.error(msg_miss)
        print(msg_miss)
        remove_handlers_if_directly_executed(func_1up)
        return None
    
    da_geop = ds_static["z"]
    da_ssgo_raw = ds_static["slor"]
    da_lsm = ds_static["lsm"]
    
    if priority == "speed":
        da_geop = da_geop.persist()
        da_ssgo_raw = da_ssgo_raw.persist()
        da_lsm = da_lsm.persist()
    
    logging.debug("Computing: land surface elevation.")
    da_era5_orog_eleva = (mpcalc.geopotential_to_height(da_geop)
                          .where(da_lsm > 0.5)
                          .where(da_geop >= 0)
                          .squeeze("time", drop=True)
                          .metpy.dequantify()
                          .rename("lse")
                         )
    
    logging.debug("Computing: slope of sub-gridscale orography.")
    da_era5_orog_slope = (da_ssgo_raw
                          .where(da_lsm > 0.5)
                          .where(da_ssgo_raw >= 0)
                          .squeeze("time", drop=True)
                          .rename("ssgo")
                         )
    
    ds_era5_orog = xr.merge([da_era5_orog_eleva, da_era5_orog_slope])
    
    # Add attributes to each DataArray within Dataset.
    
    logging.info("Adding: attributes for each DataArray within output "+
                 f"Dataset from {func_cur}.")
    for da_name in [*ds_era5_orog.keys()]:
        ds_era5_orog[da_name].attrs = copy.deepcopy(attrs_da[da_name])
    
    # Add attributes to Dataset.
    
    add_ds_attrs(ds_era5_orog, time_exec, func_cur, func_1up, args_cur, args_cur_values)
    
    # Create output file in data_processed folder.
    
    create_output_file(ds_era5_orog, path_output_orog, func_1up)
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def calc_glass_rolling_avg_of_annual_diff(region, year_start, year_end, months_subset,
                                          window_size, glass_source_pref):
    
    """
    Calculate the rolling average of the annual difference in mean leaf area index 
    (MLAI) and mean fraction of absorbed photosynthetically active radiation (MFAPAR)
    using GLASS data.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of ["ca", "sa", "wa"].
        year_start (int): Earliest year to compute the rolling average for.
        year_end (int): Latest year to compute the rolling average for.
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        window_size (int): Rolling window size (in years) to compute average for.
            Must be an odd number and greater than or equal to 3.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
                        
    Returns:
        ../data_processed/glass_rolling_avg_of_annual_diff/{calc_funcs_ver}_calc_
        {region}_{year_start}_{year_end}_{months_subset_str}_{window_size}-year_
        glass-rolling-diff_{glass_source_pref}.nc:
            Output netcdf4 file in data_processed folder containing the rolling average
            of the annual difference in MLAI and MFAPAR. {calc_funcs_ver} is the version 
            of the calc_funcs script being used. {months_subset_str} is a string 
            representing the list of selected months to use as a subset.
    
    For each grid cell, calculate the rolling average of the annual difference in 
    MLAI and MFAPAR, and only using a subset of data within this period (if a 
    months_subset not "all" is specified). These rolling averages are computed for each 
    year between year_start and year_end (inclusive). For example, the 3-year rolling 
    average of MLAI for the year 1992 would be the average of MLAI(1993)-MLAI(1992) and
    MLAI(1992)-MLAI(1991), which uses 2 annual differences across 3 years of data.
    The 5-year rolling average for the year 2002 would be the average of 
    MLAI(2004)-MLAI(2003), MLAI(2003)-MLAI(2002), MLAI(2002)-MLAI(2001) and
    MLAI(2001)-MLAI(2000), which uses 4 annual differences across 5 years of data.
    Thus, window_size must be odd for the rolling average of a year to equally
    weight years before and after it, and window_size must be greater than or equal to 3 
    for the rolling average to be well defined.
    
    This functions first runs calc_diff with calc_func = calc_glass_mean_clim for each 
    year between year_start - (window_size-1)/2 and year_end + (window_size-1)/2 
    (inclusive) which has not already been run (this is to capture all necessary years of 
    data). Then the rolling average of these annual differences is calculated for each
    year between year_start and year_end (inclusive). The calculations use 8-day satellite 
    HDF data from the data_raw folder as input, then outputs the result as a netcdf4 file 
    into the data_processed folder.
    
    Where an annual difference is completely contained within the time ranges of both AVHRR
    and MODIS data, glass_source_pref is selected as the data source for use. Otherwise, 
    AVHRR data is used where the annual difference is completely contained only within the 
    time range of AVHRR data, and conversely for MODIS data. Annual differences which 
    simultaneously cover both an AVHRR-only period (i.e. before Mar-2000) and a MODIS-only 
    period (i.e. after Dec-2018) use both AVHRR and MODIS data ("mixed").
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    # Assert that input arguments are valid.
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(region=region, year_start=year_start, year_end=year_end, 
               months_subset=months_subset, window_size=window_size, 
               glass_source_pref=glass_source_pref)
    
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Executing: {func_cur} to obtain the {months_subset_str} rolling " +
                     f"average of the annual difference in MLAI AND MFAPAR for years " +
                     f"{year_start} to {year_end}.")
    else:
        logging.info(f"Executing: {func_cur} to obtain the {months_subset_str} rolling " +
                     f"average of the annual difference in MLAI AND MFAPAR for years " +
                     f"{year_start} to {year_end} for use in {func_1up}.")
    
    # Define the output path.
    
    path_output_glass_roll = get_path_for_calc_glass_rolling(
        region=region, year_start=year_start, year_end=year_end, 
        months_subset=months_subset, window_size=window_size, 
        glass_source_pref=glass_source_pref)
    terminate_if_file_exists(path_output_glass_roll, func_cur, func_1up)
    
    # Years in which the annual difference MFAPAR DataArray is filled with NaNs.
    years_diff_nan = (set(range(avhrr_earliest_year + 1, fapar_earliest_year + 1))
                          .union(set(range(fapar_latest_year + 1, modis_latest_year + 1)))
                     )
    
    # Create intermediate output files for each annual difference if they don't already
    # exist, then read in these files, merge, and compute rolling average.
    
    datasets = []
    
    # This is a dictionary where elements will be appended onto in the following for loop.
    # It gives the glass source ("avhrr", "modis", or "mixed") for each annual difference.
    glass_source_diff = {}
    
    for year in range(int(year_start - (window_size-1)/2 + 1), 
                      int(year_end + (window_size-1)/2 + 1)
                     ):
        path_diff = get_path_for_calc_diff(
            calc_func_name="calc_glass_mean_clim", region=region, 
            period1_start="Jan-"+str(year-1), period1_end="Dec-"+str(year-1), 
            period2_start="Jan-"+str(year), period2_end="Dec-"+str(year), 
            months_subset=months_subset, glass_source_pref=glass_source_pref)
        
        glass_source_diff[str(year)] = path_diff[-8:-3]
        
        if Path(path_diff).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_diff}."
            logging.info(msg_open)
            print(msg_open)
            ds_diff = (xr.open_dataset(path_diff, engine = "netcdf4")
                       .expand_dims({"year": [year]})
                      )
            datasets.append(ds_diff)
        else:
            calc_diff(calc_func=calc_glass_mean_clim, region=region, 
                      period1_start="Jan-"+str(year-1), period1_end="Dec-"+str(year-1), 
                      period2_start="Jan-"+str(year), period2_end="Dec-"+str(year), 
                      months_subset=months_subset, glass_source_pref=glass_source_pref)
            ds_diff = (xr.open_dataset(path_diff, engine = "netcdf4")
                       .expand_dims({"year": [year]})
                      )
            datasets.append(ds_diff)
        
        if year in years_diff_nan:
            msg_avail = ("WARNING: A DataArray with NaN's was returned for the annual " +
                         f"difference from {year-1}-{year} because GLASS FAPAR data " +
                         f"is not available before {fapar_earliest_year} or after " +
                         f"{fapar_latest_year}.")
            logging.warning(msg_avail)
            print(msg_avail)
            
        if path_diff[-8:-3] == "mixed":
            msg_mixed_diff = (f"WARNING: Annual difference in mean from {year-1} to " +
                              f"{year} uses both AVHRR and MODIS data (mixed). Years " +
                              f"equal to or before {year-1} use AVHRR data while years " + 
                              f"equal to or after {year} use MODIS data.")
            logging.warning(msg_mixed_diff)
            print(msg_mixed_diff)
            year_mixed = year
        else:
            year_mixed = None
    
    logging.debug("Merging: datasets for annual difference in glass mean.")
    
    if priority == "speed":
        ds_diff_merged = xr.merge(datasets)
        ds_diff_merged = ds_diff_merged.persist()
        
    elif priority == "memory":
        ds_diff_merged = datasets[0]
        for i in range(1, len(datasets)):
            ds_diff_merged = xr.merge([ds_diff_merged, datasets[i]])
    
    # The code below computes the rolling avg by using xarray's rolling window but a more
    # efficient way would have been to recognise that the average of consecutive annual
    # differences actually produces a telescoping sum in the numerator, so that the 
    # rolling avg is in effect glass_mean(year_end) - glass_mean(year_start) divided
    # by (window_size - 1).
    
    logging.debug("Computing: rolling avg of annual difference in mean glass datasets.")
    ds_roll_diff = (ds_diff_merged
                    # xarray's rolling avg by default is assigned to coords at end of
                    # window. This line centres it on each year coord for ds_roll_diff.
                    .assign_coords({"year": ds_diff_merged.year - int((window_size-1)/2)})
                    .rolling(year = window_size - 1, min_periods = math.ceil(
                        (1-per_diff_nan_max/100) * (window_size-1)))
                    .mean(skipna = True)
                    .sel(year = slice(year_start, None))
                   )
    
    # A previous line of code found the annual difference DataArray's (for year minus 
    # year-1) which are completely filled with NaNs (i.e. the years for which year or 
    # year-1 includes data from a year for which there is no FAPAR data). The following 
    # code then computes the intersection of this set with the set of annual difference 
    # DataArray's used in calculating the rolling average centred upon each year between 
    # year_start and year_end (inclusive). This gives the number and hence percentage of
    # NaN DataArray's used in the computation of each year's rolling average.
    
    # This is a dictionary where elements will be appended onto in the following for loop.
    # It gives the glass source ("avhrr", "modis", or "mixed") for each rolling avg.
    glass_source_roll = {}
    
    for year in ds_roll_diff.year:
        year = int(year)
        # Years of annual difference DataArray's used in computing rolling avg for year.
        years_diff_used = set(range(year - int((window_size-1)/2) + 1, 
                                    year + int((window_size-1)/2) + 1)
                             )
        # Number of NaN annual difference DataArray's used.
        num_diff_nan = len(years_diff_nan
                           .intersection(years_diff_used)
                          )
        # Percentage of annual difference DataArray's used which were NaN DataArray's.
        per_diff_nan = round(num_diff_nan / (window_size - 1) * 100)
        
        if per_diff_nan > per_diff_nan_max:
            msg_nan = (f"WARNING: A DataArray with NaN's was returned for {year}'s " +
                       f"rolling avg MFAPAR since {per_diff_nan}% of necessary annual " +
                       f"difference DataArray's were NaN's (max is {per_diff_nan_max})%.")
            logging.warning(msg_nan)
            print(msg_nan)
        elif per_diff_nan > 0:
            msg_nan = (f"WARNING: {per_diff_nan}% of the annual difference DataArray's " +
                       f"used in computing {year}'s rolling avg MFAPAR were NaN's and " +
                       "ignored.")
            logging.warning(msg_nan)
            print(msg_nan)
            
        if year_mixed in years_diff_used:
            msg_mixed_roll = (f"WARNING: Rolling avg for year {year} uses both AVHRR " +
                              "and MODIS data (mixed).")
            logging.warning(msg_mixed_roll)
            print(msg_mixed_roll)
            
        if all(glass_source_diff[str(year_diff)] == "avhrr" 
               for year_diff in years_diff_used):
            glass_source_roll[str(year)] = "avhrr"
        elif all(glass_source_diff[str(year_diff)] == "modis" 
               for year_diff in years_diff_used):
            glass_source_roll[str(year)] = "modis"
        else:
            glass_source_roll[str(year)] = "mixed"
    
    # Add attributes to each DataArray within Dataset.
    
    logging.info("Adding: attributes for each DataArray within output "+
                 f"Dataset from {func_cur}.")
    for da_name in [*ds_roll_diff.keys()]:
        ds_roll_diff[da_name].attrs = copy.deepcopy(attrs_da[da_name])
        ds_roll_diff[da_name].attrs["abbreviation"] = (
            "$annual \ difference _{{{} \ year}} ^{{roll \ avg}}$({})"
            .format(window_size, ds_roll_diff[da_name].attrs["abbreviation"])
        )
        ds_roll_diff[da_name].attrs["full_name"] = (
            f"{window_size}-Year Rolling Average of Annual Difference in " + 
            ds_roll_diff[da_name].attrs["full_name"]
        )
        ds_roll_diff[da_name].attrs.update(glass_source_roll)
    
    # Add attributes to Dataset.
    
    add_ds_attrs(ds_roll_diff, time_exec, func_cur, func_1up, args_cur, args_cur_values)
    
    # Create output file in data_processed folder.
    
    create_output_file(ds_roll_diff, path_output_glass_roll, func_1up)
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


## High level calculation functions to create all possible data files.

def create_all_possible_calc_data_files(region, period_start, period_end, months_subset):
    
    """
    For the given inputs, run all possible calculation functions. For calculation
    functions which require a var_or_dvar input, the function is run for every
    var_or_dvar possible (i.e. all items in the vars_and_dvars_era5_all list
    in the global settings).
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(region=region, period_start=period_start, period_end=period_end,
               months_subset=months_subset)
    
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Executing: {func_cur} to create all possible calc data files " +
                     f"in {region} over {months_subset_str} months between " + 
                     f"{period_start} and {period_end}.")
    else:
        logging.info(f"Executing: {func_cur} to create all possible calc data files " +
                     f"in {region} over {months_subset_str} months between " + 
                     f"{period_start} and {period_end} for use in {func_1up}.")
    
    for calc_func_name in calc_func_names:
        calc_func = globals()[calc_func_name]
        
        if calc_func_name == "calc_era5_wsd_clim":
            path = get_path_for_calc_func(
                calc_func_name=calc_func_name, region=region, period_start=period_start, 
                period_end=period_end, months_subset=months_subset
            )
            if Path(path).exists():
                msg_exists = f"Skipped: file already exists: {path}."
                logging.info(msg_exists)
                print(msg_exists)
            else:
                calc_func(
                    region=region, period_start=period_start, period_end=period_end,
                    months_subset=months_subset
                )
        elif calc_func_name == "calc_glass_mean_clim":
            for glass_source in glass_sources_all:
                path = get_path_for_calc_func(
                    calc_func_name=calc_func_name, region=region, 
                    period_start=period_start, period_end=period_end, 
                    months_subset=months_subset, glass_source_pref=glass_source
                )                
                if Path(path).exists():
                    msg_exists = f"Skipped: file already exists: {path}."
                    logging.info(msg_exists)
                    print(msg_exists)
                else:
                    calc_func(
                        region=region, period_start=period_start, period_end=period_end,
                        months_subset=months_subset, glass_source_pref=glass_source
                    )
        else:
            for var_or_dvar in vars_and_dvars_era5_all:
                path = get_path_for_calc_func(
                    calc_func_name=calc_func_name, region=region, 
                    period_start=period_start, period_end=period_end, 
                    months_subset=months_subset, var_or_dvar=var_or_dvar
                )                
                if Path(path).exists():
                    msg_exists = f"Skipped: file already exists: {path}."
                    logging.info(msg_exists)
                    print(msg_exists)
                else:
                    calc_func(
                        region=region, period_start=period_start, period_end=period_end,
                        months_subset=months_subset, var_or_dvar=var_or_dvar
                    )
    
    remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def create_all_possible_diff_data_files(region, period1_start, period1_end,
                                        period2_start, period2_end, months_subset):
    
    """
    For the given inputs, run calc_diff over all possible calculation functions. 
    For calculation functions which require a var_or_dvar input, the function 
    is run for every var_or_dvar possible (i.e. all items in the vars_and_dvars_era5_all 
    list in the global settings).
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period1_start (str): Start of first period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period1_end (str): End of first period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period2_start (str): Start of second period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period2_end (str): End of second period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        months_subset (str or list): Subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                    args_cur, args_cur_values)
    
    check_args_for_none(func_cur, args_cur, args_cur_values)
    check_args(region=region, period1_start=period1_start, period1_end=period1_end, 
               period2_start=period2_start, period2_end=period2_end,
               months_subset=months_subset)
    
    months_subset_str = get_months_subset_str(months_subset=months_subset)
    
    if (func_1up == "<cell line: 1>") | (func_1up == "<module>"):
        logging.info(f"Executing: {func_cur} to create all possible {months_subset_str} " +
                     f"difference data files in {region} over {period1_start} to " +
                     f"{period1_end} and {period2_start} to {period2_end}.")
    else:
        logging.info(f"Executing: {func_cur} to create all possible {months_subset_str} " +
                     f"difference data files in {region} over {period1_start} to " +
                     f"{period1_end} and {period2_start} to {period2_end} " +
                     f"for use in {func_1up}.")
    
    for calc_func_name in calc_func_names:
        calc_func = globals()[calc_func_name]
        
        if calc_func_name == "calc_era5_wsd_clim":
            path = get_path_for_calc_diff(
                calc_func_name=calc_func_name, region=region, 
                period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                months_subset=months_subset
            )
            if Path(path).exists():
                msg_exists = f"Skipped: file already exists: {path}."
                logging.info(msg_exists)
                print(msg_exists)
            else:
                calc_diff(
                    calc_func=calc_func, region=region, 
                    period1_start=period1_start, period1_end=period1_end, 
                    period2_start=period2_start, period2_end=period2_end, 
                    months_subset=months_subset
                )
        elif calc_func_name == "calc_glass_mean_clim":
            for glass_source in glass_sources_all:       
                path = get_path_for_calc_diff(
                    calc_func_name=calc_func_name, region=region, 
                    period1_start=period1_start, period1_end=period1_end, 
                    period2_start=period2_start, period2_end=period2_end, 
                    months_subset=months_subset, glass_source_pref=glass_source
                )
                if Path(path).exists():
                    msg_exists = f"Skipped: file already exists: {path}."
                    logging.info(msg_exists)
                    print(msg_exists)
                else:
                    calc_diff(
                        calc_func=calc_func, region=region, 
                        period1_start=period1_start, period1_end=period1_end, 
                        period2_start=period2_start, period2_end=period2_end, 
                        months_subset=months_subset, glass_source_pref=glass_source
                    )
        else:
            for var_or_dvar in vars_and_dvars_era5_all:       
                path = get_path_for_calc_diff(
                    calc_func_name=calc_func_name, region=region, 
                    period1_start=period1_start, period1_end=period1_end, 
                    period2_start=period2_start, period2_end=period2_end, 
                    months_subset=months_subset, var_or_dvar=var_or_dvar
                )
                if Path(path).exists():
                    msg_exists = f"Skipped: file already exists: {path}."
                    logging.info(msg_exists)
                    print(msg_exists)
                else:
                    calc_diff(
                        calc_func=calc_func, region=region, 
                        period1_start=period1_start, period1_end=period1_end, 
                        period2_start=period2_start, period2_end=period2_end, 
                        months_subset=months_subset, var_or_dvar=var_or_dvar
                    )
    
    remove_handlers_if_directly_executed(func_1up)
