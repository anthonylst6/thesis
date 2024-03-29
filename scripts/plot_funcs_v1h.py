#!/usr/bin/env python
# coding: utf-8

# # Setup

# ## Import libraries for plotting

# In[ ]:


import cmocean
import inspect
import logging
import copy
import importlib
import math
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob
from pathlib import Path
from datetime import datetime
from dateutil.relativedelta import relativedelta
from textwrap import wrap


# ## Import calc_funcs module for use in plotting

# In[ ]:


# Choose calc_funcs_ver to use in plotting scripts.
cfv = "latest"

assert (isinstance(cfv, str) & (len(cfv) == 5) & (cfv[:3] == "cfv") & 
        cfv[3].isnumeric() & cfv[4].isalpha() & cfv[4].islower()) | (cfv == "latest"), \
    ("cfv must be 'latest' or of form 'cfvXY' where X is a single digit number " +
     "and Y is a lowercase alphabet character. eg. cfv1n")

if cfv == "latest":
    calc_funcs_scripts = glob("calc_funcs_*.py")
    calc_funcs_scripts.sort()
    calc_funcs_module = calc_funcs_scripts[-1][:-3]
    
else:
    calc_funcs_module = "calc_funcs_" + cfv[2:]
    
cf = importlib.import_module(calc_funcs_module)

print(f"Using: {calc_funcs_module}")


# ## Settings and global variables for plotting

# In[ ]:


# This is to get the version number of the plot_funcs script being used so that it
# can be appended to the file name of any outputs. The reason this is done is because
# the plot functions below and the calculation functions in the calc_funcs script both
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
    plot_funcs_ver = "pf" + Path(__file__).stem[-3:]
except:
    plot_funcs_ver = "pfv00"

# Set level of logging out of (in decreasing order of detail): 
# [logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL]
plot_log_level = logging.INFO
assert plot_log_level in cf.log_levels, \
    f"[plot_log_level (global variable in settings) must be one of: {cf.log_levels}"
cf.calc_log_level = plot_log_level

# Specify valid dimension names for data array input into plot functions.
# This is used in asserts later to prevent accidentally inputting the wrong data array.
da_dims_valid = ("latitude", "longitude")

# Categorise parameters for use in automating colourbars and extents in plot functions.
da_names_cyclic = ["hour_max", "hour_min"]
da_names_pos_with_vmin_0 = ["lse", "ssgo"] + cf.params_glass_mean
da_names_pos = ["range"] + cf.params_wsd
vars_pos_with_vmin_0 = []
vars_pos = ["ws10", "ws100", "mslp", "t2", 
            "vipile", "vike", "tcclw", "tcwv", 
            "blh", "fa", "cbh", "tcc", "ci"]

# Name of top level functions to create all possible plots. High level functions will
# read their call stack if it is being called from one of these top level functions,
# it will know to skip a plot if the file for it already exists (saves on computation).
funcs_create_all_plot = ["create_all_possible_calc_plot_files", 
                         "create_all_possible_diff_plot_files", 
                         "create_all_possible_comp_plot_files"]

# Standard figure width for figures which contain an individual plot. Figures with
# multiple subplots will have a width double this.
figwidth_standard = 10

# Standard width of subfigure titles (beyond this, text will wrap itself to a new line).
title_width = figwidth_standard * 6

# Customise size and shape of quivers in quiver plots.
quiver_scale_multiplier = 10
quiver_headwidth = 4.5

# Standardise bar widths for bars in climate indices plots.
bar_width = 31

# Set range of values where the symmetrical log plot is linear rather than log.
eroe100_linthresh = 1e-20

# Quantile of values to mask for percentage difference plots (this is used to avoid
# extremely large percentage differences if it was coming off a very low base value.
mask_perc_quantile_default = 10

# General plot settings (enable use of TeX and set resolution of output figures).
plt.rcParams['text.usetex'] = True
plt.rcParams['savefig.dpi'] = 300

# SMALL_SIZE = 14
# MEDIUM_SIZE = 16
# BIGGER_SIZE = 18

# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


# ## Dates for meteorological events

# In[ ]:


# Dates below are in [event_start, event_end] form.

# Dates for Negative and Positive Indian Ocean Dipole (IOD) according
# to historical data by the Japanese Meteorological Agency (JMA). Note
# that events are defined differently across meteorological agencies.
# JMA data is selected here since it was the only source of monthly
# data which the author could easily access and was up to date.
# https://ds.data.jma.go.jp/tcc/tcc/products/elnino/iodevents.html

dates_neg_iod = [
    ["Jul-1952", "Sep-1952"],
    ["Aug-1954", "Oct-1954"],
    ["Jun-1956", "Aug-1956"],
    ["Jun-1958", "Oct-1958"],
    ["Sep-1975", "Nov-1975"],
    ["Jul-1981", "Sep-1981"],
    ["Jun-1984", "Nov-1984"],
    ["Jun-1985", "Aug-1985"],
    ["Jul-1989", "Oct-1989"],
    ["Aug-1996", "Nov-1996"],
    ["Aug-1998", "Nov-1998"],
    ["Aug-2010", "Nov-2010"],
    ["Jun-2013", "Sep-2013"],
    ["Jun-2016", "Nov-2016"],
    ["Aug-2020", "Oct-2020"],
    ["Jun-2021", "Nov-2021"],
    # Negative IOD still in progress at time of writing
    ["Jun-2022", "Oct-2022"],
]

dates_pos_iod = [
    ["Jun-1961", "Nov-1961"],
    ["Jul-1963", "Nov-1963"],
    ["Jun-1967", "Nov-1967"],
    ["Jun-1972", "Nov-1972"],
    ["Aug-1982", "Nov-1982"],
    ["Jun-1994", "Nov-1994"],
    ["Jul-1997", "Nov-1997"],
    ["Aug-2006", "Nov-2006"],
    ["Jun-2007", "Sep-2007"],
    ["Jun-2008", "Aug-2008"],
    ["Aug-2011", "Oct-2011"],
    ["Jul-2012", "Oct-2012"],
    ["Aug-2015", "Nov-2015"],
    ["Jun-2019", "Nov-2019"],
]

# Dates for La Nina and El Nino according to historical data ny the
# Japanese Meteorological Agency (JMA). Note that La Nina and El Nino
# events are defined differently across meteorological agencies.
# JMA data is selected here for consistency with choice of IOD data.
# https://ds.data.jma.go.jp/tcc/tcc/products/elnino/ensoevents.html

dates_la_nina = [
    ["Apr-1954", "Feb-1956"],
    ["Apr-1964", "Jan-1965"],
    ["Sep-1967", "Apr-1968"],
    ["May-1970", "Dec-1971"],
    ["Jun-1973", "Mar-1974"],
    ["Apr-1975", "Mar-1976"],
    ["Jul-1984", "Sep-1985"],
    ["Apr-1988", "May-1989"],
    ["Jul-1995", "Feb-1996"],
    ["Aug-1998", "Apr-2000"],
    ["Oct-2005", "Mar-2006"],
    ["Apr-2007", "Apr-2008"],
    ["Jul-2010", "Mar-2011"],
    ["Sep-2017", "Mar-2018"],
    ["Jul-2020", "Apr-2021"],
    # La Nina still in progress at time of writing
    ["Sep-2021", "Oct-2022"],
]

dates_el_nino = [
    ["May-1951", "Feb-1952"],
    ["Apr-1953", "Nov-1953"],
    ["Apr-1957", "Apr-1958"],
    ["Jun-1963", "Jan-1964"],
    ["May-1965", "Feb-1966"],
    ["Sep-1968", "Feb-1970"],
    ["May-1972", "Mar-1973"],
    ["Jun-1976", "Mar-1977"],
    ["Apr-1982", "Aug-1983"],
    ["Sep-1986", "Jan-1988"],
    ["Apr-1991", "Jul-1992"],
    ["Apr-1997", "May-1998"],
    ["Jun-2002", "Feb-2003"],
    ["Jun-2009", "Mar-2010"],
    ["Jun-2014", "Apr-2016"],
    ["Sep-2018", "May-2019"],
]


# # Functions

# ## Supplementary functions for plotting

# In[ ]:


def get_plot_metadata(time_exec_1up, func_1up, args_1up, args_1up_values):
    
    """
    Obtain a string representing the function and input arguments used for creating a
    plot, as well as the date and time of creation. Designed to be called within the
    output section of plot functions to append this metadata to output PNG plots.
    Appended metadata can then be viewed using the exiftool command line tool.
    
    Arguments:
        time_exec_1up (datetime.datetime): Time of when func_1up was executed.
        func_1up (str): Name of function calling this function.
        args_1up (list): List of argument names for func_1up.
        args_1up_values (dict): Mapping of argument names for func_1up to their
            input values.
    
    Returns:
        {"History": "{func_1up}({args_1up_str})_{time_str}"} (dict):
            Dictionary mapping the "History" attribute of PNG plots to the metadata
            string. {args_1up_str} is the string representation of arguments input 
            into {func_1up}. {time_str} is a string giving the time of when 
            {func_1up} was executed, in "%Y-%m-%d-%H-%M-%S" format.
    """
    
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
    
    return {"History": f"{func_1up}({args_1up_str})_{time_str}"}


# In[ ]:


def apply_mask(da, frame_comp):
    
    """
    Return a DataArray where values have been masked according to specifications
    in a comparison plot function as to whether we should mask positive or negative
    summary statistics computed over each period.
    
    Arguments:
        da (xarray.DataArray): Input DataArray to be masked.
        frame_comp (frame): Python stack frame containing mask specifications for the 
        comparison plot function which is calling this function.
        
    Returns:
        da_masked (xarray.DataArray): DataArray where values have been masked.
    """
    
    da_masked = copy.deepcopy(da)
    
    # Obtain mask specifications by inspecting arguments of frame_comp.
    _, _, _, args_comp_values = inspect.getargvalues(frame_comp)
    calc_func_name = args_comp_values["calc_func"].__name__
    region = args_comp_values["region"]
    period1_start = args_comp_values["period1_start"]
    period1_end = args_comp_values["period1_end"]
    period2_start = args_comp_values["period2_start"]
    period2_end = args_comp_values["period2_end"]
    period1_months = args_comp_values["period1_months"]
    period2_months = args_comp_values["period2_months"]
    period1_hours = args_comp_values["period1_hours"]
    period2_hours = args_comp_values["period2_hours"]
    arg_extra = args_comp_values["arg_extra"]
    glass_source_pref = args_comp_values["glass_source_pref"]
    var_or_dvar = args_comp_values["var_or_dvar"]
    mask_period1 = args_comp_values["mask_period1"]
    mask_period2 = args_comp_values["mask_period2"]
    cfv_data = args_comp_values["cfv_data"]
    
    # Obtain paths to calc_func outputs (containing summary statistics) for each period.
    path_period1 = cf.get_path_for_calc_func(
        calc_func_name=calc_func_name, region=region, period_start=period1_start, 
        period_end=period1_end, period_months=period1_months, period_hours=period1_hours,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar
    )
    path_period2 = cf.get_path_for_calc_func(
        calc_func_name=calc_func_name, region=region, period_start=period2_start, 
        period_end=period2_end, period_months=period2_months, period_hours=period2_hours,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar
    )
    
    # Use data outputted from an older version of the calc_funcs script. This is useful
    # for results which required computationally intensive processing. And can also be
    # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
    if cfv_data:
        for path in [path_period1, path_period2]:
            path = path.replace(cf.calc_funcs_ver, cfv_data)
    
    # Open datasets using previously obtained calc_func output paths.
    ds_period1 = xr.open_dataset(path_period1, engine = "netcdf4")
    ds_period2 = xr.open_dataset(path_period2, engine = "netcdf4")
    
    # Select the relevant DataArray's contained within each opened dataset.
    if calc_func_name == "calc_era5_mdp_clim_given_var_or_dvar":
        if var_or_dvar in cf.params_vector:
            pass
        else:
            da_period1 = ds_period1[var_or_dvar].sel(hour=arg_extra)
            da_period2 = ds_period2[var_or_dvar].sel(hour=arg_extra)
    elif (var_or_dvar in cf.params_vector) & (arg_extra in ["max", "min", "mean"]):
        pass
    else:
        da_period1 = ds_period1[arg_extra]
        da_period2 = ds_period2[arg_extra]
    
    # Obtain name of main parameter we are plotting for.
    main_param = (da.attrs["abbreviation"]
                  .split("(")[-1]
                  .split(")")[0]
                  .split("^")[0]
                  .split("$")[0]
                  .lower())
    
    # Mask values only for parameters which are not cyclic and which can
    # take on both positive and negative values.
    if da.name in da_names_cyclic:
        pass
    elif da.name in da_names_pos_with_vmin_0:
        pass
    elif da.name in da_names_pos:
        pass
    elif main_param in vars_pos_with_vmin_0:
        pass
    elif main_param in vars_pos:
        pass
    else:
        if mask_period1 == "pos":
            da_masked = da_masked.where(da_period1 < 0)
        elif mask_period1 == "neg":
            da_masked = da_masked.where(da_period1 >= 0)
        if mask_period2 == "pos":
            da_masked = da_masked.where(da_period2 < 0)
        elif mask_period2 == "neg":
            da_masked = da_masked.where(da_period2 >= 0)
    
    return da_masked


# In[ ]:


def get_common_cbar_limits(
    calc_func, region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, arg_extra, 
    period1_hours=None, period2_hours=None, glass_source_pref=None, var_or_dvar=None, 
    extents=None, cfv_data=None
):
    
    """
    Obtain common colourbar limits for each period in a comparison plot.
    
    Arguments:
        calc_func (function): Calculation function to use in analysis. Must be one of: 
            [calc_glass_mean_clim,
            calc_era5_mdp_clim_given_var_or_dvar,
            calc_era5_mdp_clim_stats_given_var_or_dvar,
            calc_era5_mean_clim_given_var_or_dvar,
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        arg_extra (str or int): Extra plotting argument used to specify which GLASS 
            parameter to plot, which hour for the mean diurnal profile of an ERA5 
            parameter to plot, which statistic of the mean diurnal profile to plot, 
            or which parameter of the wind speed distribution to plot. Must be one of:
            ["mlai", "mfapar", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
            16, 17, 18, 19, 20, 21, 22, 23, "hour_max", "hour_min", "max", "max_u", 
            "max_v", "min", "min_u", "min_v", "mean", "mean_u", "mean_v", "range", 
            "ws10_mean", "ws10_std", "c10", "k10", "ws100_mean", "ws100_std", "c100", 
            "k100", "eroe100", "tgcf100"].
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
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
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        
    Returns:
        (vmin, vmax) (tuple): Common colourbar limits for each period in comparison.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    frame_1up = frame_cur.f_back
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    calc_func_name = calc_func.__name__
    
    # Assert that input arguments are valid, and obtain paths to calc_func outputs
    # (containing summary statistics) for each period.
    
    cf.check_args_for_none(calc_func_name, args_cur, args_cur_values)
    cf.check_args(calc_func=calc_func, region=region, 
                  period1_start=period1_start, period1_end=period1_end, 
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_months=period1_months, period2_months=period2_months, 
                  period1_hours=period1_hours, period2_hours=period2_hours, 
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar, 
                  arg_extra=arg_extra, extents=extents, cfv_data=cfv_data)
        
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
        
    path_period1 = cf.get_path_for_calc_func(
        calc_func_name=calc_func_name, region=region, period_start=period1_start, 
        period_end=period1_end, period_months=period1_months, period_hours=period1_hours,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar
    )
    path_period2 = cf.get_path_for_calc_func(
        calc_func_name=calc_func_name, region=region, period_start=period2_start, 
        period_end=period2_end, period_months=period2_months, period_hours=period2_hours,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar
    )
    
    # Use data outputted from an older version of the calc_funcs script. This is useful
    # for results which required computationally intensive processing. And can also be
    # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
    if cfv_data:
        for path in [path_period1, path_period2]:
            path = path.replace(cf.calc_funcs_ver, cfv_data)
            if Path(path).exists() == False:
                msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                           f"but could not find file: {path}")
                logging.error(msg_cfv)
                print(msg_cfv)
                cf.remove_handlers_if_directly_executed(func_1up)
                raise Exception(msg_cfv)
    
    # Open calc_func output files if they exist, otherwise create the output files.
    
    if Path(path_period1).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_period1}"
        logging.info(msg_open)
        print(msg_open)
    else:
        calc_func(region=region, period_start=period1_start, period_end=period1_end, 
                  period_months=period1_months, period_hours=period1_hours, 
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    
    ds_period1 = (xr.open_dataset(path_period1, engine = "netcdf4")
                  .sel(longitude=slice(extents[0], extents[1]), 
                       latitude=slice(extents[3], extents[2]))
                 )
    
    if Path(path_period2).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_period2}"
        logging.info(msg_open)
        print(msg_open)
    else:
        calc_func(region=region, period_start=period2_start, period_end=period2_end, 
                  period_months=period2_months, period_hours=period2_hours,
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    
    ds_period2 = (xr.open_dataset(path_period2, engine = "netcdf4")
                  .sel(longitude=slice(extents[0], extents[1]), 
                       latitude=slice(extents[3], extents[2]))
                 )
    
    # Select out relevant DataArray's from dataset, then assign colourbar limits based
    # on maximum and minimum values in each DataArray, depending on what category
    # of parameters the main parameter we are plotting for falls into.
    
    if calc_func_name == "calc_era5_mdp_clim_given_var_or_dvar":
        if var_or_dvar in cf.params_vector:
            da_u_period1 = ds_period1[var_or_dvar.replace("wv", "u")].sel(hour=arg_extra)
            da_v_period1 = ds_period1[var_or_dvar.replace("wv", "v")].sel(hour=arg_extra)
            da_u_period2 = ds_period2[var_or_dvar.replace("wv", "u")].sel(hour=arg_extra)
            da_v_period2 = ds_period2[var_or_dvar.replace("wv", "v")].sel(hour=arg_extra)
            da_mag_period1 = cf.get_magnitude(da_u_period1, da_v_period1)
            da_mag_period2 = cf.get_magnitude(da_u_period2, da_v_period2)
            vmin = float(min(da_mag_period1.min(), da_mag_period2.min()))
            vmax = float(max(da_mag_period1.max(), da_mag_period2.max()))
        else:
            da_period1 = ds_period1[var_or_dvar].sel(hour=arg_extra)
            da_period2 = ds_period2[var_or_dvar].sel(hour=arg_extra)
            if func_1up == "create_individual_comp_plot":
                da_period1 = apply_mask(da_period1, frame_1up)
                da_period2 = apply_mask(da_period2, frame_1up)
            main_param_period2 = (da_period2.attrs["abbreviation"]
                                  .split("(")[-1]
                                  .split(")")[0]
                                  .split("^")[0]
                                  .split("$")[0]
                                  .lower())
            if da_period2.name in da_names_cyclic:
                vmin = None
                vmax = None
            elif da_period2.name in da_names_pos_with_vmin_0:
                vmin = 0
                vmax = float(max(da_period1.max(), da_period2.max()))
            elif da_period2.name in da_names_pos:
                vmin = float(min(da_period1.min(), da_period2.min()))
                vmax = float(max(da_period1.max(), da_period2.max()))
            elif main_param_period2 in vars_pos_with_vmin_0:
                vmin = 0
                vmax = float(max(da_period1.max(), da_period2.max()))
            elif main_param_period2 in vars_pos:
                vmin = float(min(da_period1.min(), da_period2.min()))
                vmax = float(max(da_period1.max(), da_period2.max()))
            else:
                min_of_mins = float(min(da_period1.min(), da_period2.min()))
                max_of_maxs = float(max(da_period1.max(), da_period2.max()))
                vmin = min(-abs(min_of_mins), -abs(max_of_maxs))
                vmax = -vmin
    elif (var_or_dvar in cf.params_vector) & (arg_extra in ["max", "min", "mean"]):
        da_u_period1 = ds_period1[arg_extra + "_u"]
        da_v_period1 = ds_period1[arg_extra + "_v"]
        da_u_period2 = ds_period2[arg_extra + "_u"]
        da_v_period2 = ds_period2[arg_extra + "_v"]
        da_mag_period1 = cf.get_magnitude(da_u_period1, da_v_period1)
        da_mag_period2 = cf.get_magnitude(da_u_period2, da_v_period2)
        vmin = float(min(da_mag_period1.min(), da_mag_period2.min()))
        vmax = float(max(da_mag_period1.max(), da_mag_period2.max()))
    else:
        da_period1 = ds_period1[arg_extra]
        da_period2 = ds_period2[arg_extra]
        if func_1up == "create_individual_comp_plot":
            da_period1 = apply_mask(da_period1, frame_1up)
            da_period2 = apply_mask(da_period2, frame_1up)
        main_param_period2 = (da_period2.attrs["abbreviation"]
                              .split("(")[-1]
                              .split(")")[0]
                              .split("^")[0]
                              .split("$")[0]
                              .lower())
        if da_period2.name in da_names_cyclic:
            vmin = None
            vmax = None
        elif da_period2.name in da_names_pos_with_vmin_0:
            vmin = 0
            vmax = float(max(da_period1.max(), da_period2.max()))
        elif da_period2.name in da_names_pos:
            vmin = float(min(da_period1.min(), da_period2.min()))
            vmax = float(max(da_period1.max(), da_period2.max()))
        elif main_param_period2 in vars_pos_with_vmin_0:
            vmin = 0
            vmax = float(max(da_period1.max(), da_period2.max()))
        elif main_param_period2 in vars_pos:
            vmin = float(min(da_period1.min(), da_period2.min()))
            vmax = float(max(da_period1.max(), da_period2.max()))
        else:
            min_of_mins = float(min(da_period1.min(), da_period2.min()))
            max_of_maxs = float(max(da_period1.max(), da_period2.max()))
            vmin = min(-abs(min_of_mins), -abs(max_of_maxs))
            vmax = -vmin
    
    # If the minimum or maximum returns NaNs, then just set the limits to None
    # instead and let the default plotting algorithms for matplotlib handle it.
    if vmin != None:
        if math.isnan(vmin):
            vmin = None
    if vmax != None:
        if math.isnan(vmax):
            vmax = None
    
    cf.remove_handlers_if_directly_executed(func_1up)
    
    return vmin, vmax


# In[ ]:


def round_down_first_sig(value_to_be_rounded):
    
    """
    Round a value downwards to the nearest value with one significant figure.
    
    Arguments:
        value_to_be_rounded (float): Value to be rounded.
        
    Returns:
        value_rounded (float): Value which has been rounded.
    """
    
    # The goal here is to select out the decimal place and value for the first
    # significant figure, round the input value by that decimal place, then
    # divide by the decimal value. For example, for 0.0647, the decimal place
    # and value for the first sig fig is 2 and 6 respectively. Rounding 0.647
    # to 2nd decimal place then gives 0.6, and dividing this by 5 gives 0.1.
    
    # If the value for the second sig fig is greater than or equal to 5, then
    # we divide by the decimal value + 1 instead. Eg. for 0.0657, rounding gives 
    # 0.07, so we would need to divide by 7 rather 6 to obtain 0.01.
    
    dec_place_of_first_sig = -int(math.floor(math.log10(abs(value_to_be_rounded))))
    
    # We use the format function to express the input value as a string then
    # extract the value given the decimal place. However, this strategy fails
    # if all the sig figs happen to be 9 and the last digit is greater than or
    # equal to 5 (because the format function would round it up). Therefore,
    # we extend out the string to 40 sig figs, as the chances of 39 consecutive
    # 9's are negligible.
    
    str_with_40_sig = format(value_to_be_rounded, f".{dec_place_of_first_sig+39}f")
    dec_value_of_first_sig = int(str_with_40_sig[-40])
    dec_value_of_second_sig = int(str_with_40_sig[-39])

    if dec_value_of_second_sig >= 5:
        value_rounded = (round(value_to_be_rounded, dec_place_of_first_sig) / 
                         (dec_value_of_first_sig + 1))
    else:
        value_rounded = (round(value_to_be_rounded, dec_place_of_first_sig) / 
                         dec_value_of_first_sig)
        
    return value_rounded


# ## Low-level plotting functions

# In[ ]:


def create_pcolormesh(da, extents=None, vmin=None, vmax=None, ax=None):
    
    """
    Create a plot where grid cell values are coloured according to their value.
    
    Arguments:
        da (xarray.DataArray): DataArray to create plot for.
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        vmin (float or int): Minimum of colourbar extents for a calc plot.
        vmax (float or int): Maximum of colourbar extents for a calc plot.
        ax (cartopy.GeoAxesSubplot): Figure axis to create plot on.
        
    Returns:
        A pcolormesh plot.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid.
    
    assert ((str(type(da)) == "<class 'xarray.core.dataarray.DataArray'>") & 
            (da.dims == da_dims_valid)), \
        f"da must be an xarray.DataArray with da.dims == {da_dims_valid}"
    cf.check_args(extents=extents, vmin=vmin, vmax=vmax, ax=ax)
    
    # Obtain frame 2 steps up in call stack for to handle masking of values
    # in case this function is called for use in a comparison plot function.
    
    frame_2up = frame_cur.f_back.f_back
    func_2up = inspect.getframeinfo(frame_2up)[2]
    
    # If no extents specified, use extents of coordinates in DataArray.
    
    if extents == None:
        extents = []
        extents.append(da.longitude.min())
        extents.append(da.longitude.max())
        extents.append(da.latitude.min())
        extents.append(da.latitude.max())
    
    # Create copy of input ax to refer to when deciding whether to display plot.
    
    ax_input = ax
    
    # If no ax is specified (i.e. this function is not being used to create a subplot
    # within a figure for a different plot function), then create a figure.
    
    if ax == None:
        figrows = 1
        figcols = 1
        figwidth = figwidth_standard
        figheight = (figwidth * figrows/figcols * 
                     (extents[3]-extents[2])/(extents[1]-extents[0])
                    )
        fig, ax = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                               subplot_kw = {"projection": ccrs.PlateCarree()}
                              )
    
    # Obtain name of main parameter we will be plotting.
    
    main_param = (da.attrs["abbreviation"]
                  .split("(")[-1]
                  .split(")")[0]
                  .split("^")[0]
                  .split("$")[0]
                  .lower())
    
    # Specify colourmap, colourbar extents and levels in colourbar depending
    # on category of parameter we are plotting.
    
    levels = None
    
    if da.attrs["full_name"].split(" ")[1] == "Rolling":
        cmap = cmocean.cm.balance_r
        if (vmin == None) & (vmax == None):
            da_subset = da.sel(longitude=slice(extents[0], extents[1]), 
                               latitude=slice(extents[3], extents[2]))
            min_of_da = float(da_subset.min())
            max_of_da = float(da_subset.max())
            vmin = min(-abs(min_of_da), -abs(max_of_da))
            vmax = -vmin
    elif da.attrs["full_name"].split(" ")[0] == "Difference":
        if da.name in da_names_cyclic:
            cmap = "twilight_shifted_r"
            levels = np.arange(-12, 12+1)
        else:
            cmap = cmocean.cm.balance_r
            if func_2up == "create_individual_comp_plot":
                da = apply_mask(da, frame_2up)
            if (vmin == None) & (vmax == None):
                da_subset = da.sel(longitude=slice(extents[0], extents[1]), 
                                   latitude=slice(extents[3], extents[2]))
                min_of_da = float(da_subset.min())
                max_of_da = float(da_subset.max())
                vmin = min(-abs(min_of_da), -abs(max_of_da))
                vmax = -vmin
    else:
        if da.name in da_names_cyclic:
            cmap = cmocean.cm.phase
            levels = np.arange(0, 24+1)
        elif da.name in da_names_pos_with_vmin_0:
            da_subset = da.sel(longitude=slice(extents[0], extents[1]), 
                               latitude=slice(extents[3], extents[2]))
            vmin = vmin if vmin else 0
            vmax = vmax if vmax else da_subset.max()
            if da.name  == "mlai":
                cmap = cmocean.cm.algae
            elif da.name == "lse":
                cmap = cmocean.tools.crop(cmocean.cm.topo, vmin, vmax, 0)
            else:
                cmap = "viridis"
        elif da.name in da_names_pos:
            cmap = "viridis"
            da_subset = da.sel(longitude=slice(extents[0], extents[1]), 
                               latitude=slice(extents[3], extents[2]))
            vmin = vmin if vmin else da_subset.min()
            vmax = vmax if vmax else da_subset.max()
        elif main_param in vars_pos_with_vmin_0:
            cmap = "viridis"
            da_subset = da.sel(longitude=slice(extents[0], extents[1]), 
                               latitude=slice(extents[3], extents[2]))
            vmin = vmin if vmin else 0
            vmax = vmax if vmax else da_subset.max()
        elif main_param in vars_pos:
            cmap = "viridis"
            da_subset = da.sel(longitude=slice(extents[0], extents[1]), 
                               latitude=slice(extents[3], extents[2]))
            vmin = vmin if vmin else da_subset.min()
            vmax = vmax if vmax else da_subset.max()
        else:
            cmap = "PuOr"
            if func_2up == "create_individual_comp_plot":
                da = apply_mask(da, frame_2up)
            if (vmin == None) & (vmax == None):
                da_subset = da.sel(longitude=slice(extents[0], extents[1]), 
                                   latitude=slice(extents[3], extents[2]))
                min_of_da = float(da_subset.min())
                max_of_da = float(da_subset.max())
                vmin = min(-abs(min_of_da), -abs(max_of_da))
                vmax = -vmin
    
    if vmin != None:
        if math.isnan(vmin):
            vmin = None
    if vmax != None:
        if math.isnan(vmax):
            vmax = None
    
    # Create plot
    
    ax.set_extent(extents=extents, crs=ccrs.PlateCarree())
    
    if da.attrs["full_name"].split(" ")[1] == "Rolling":
        cbar_label = "{}".format(da.attrs["abbreviation"])
    else:
        cbar_label = "{} [{}]".format(da.attrs["abbreviation"], da.attrs["units"])
        
    # The units = dimensionless condition is used so that this line is invoked only if
    # perc = False in the mid and high-level plotting function arguments.
    if (da.name == "eroe100") & (da.attrs["units"] == "dimensionless"):
        if da.attrs["full_name"].split(" ")[0] == "Difference":
            # The conditional reassignment of vmin and vmax below is to avoid cbar
            # extents falling within the linear range and producing ugly graphs.
            if (vmin != None) & (vmax != None):
                if vmax < eroe100_linthresh:
                    vmin = round_down_first_sig(vmin)
                    vmax = round_down_first_sig(vmax)
            da.plot.pcolormesh(ax = ax, cmap = cmap, transform = ccrs.PlateCarree(),
                               norm = colors.SymLogNorm(linthresh=eroe100_linthresh, 
                                                        vmin=vmin, vmax=vmax),
                               extend = "both", cbar_kwargs = {"label": cbar_label}
                              )
        else:
            da.plot.pcolormesh(ax = ax, cmap = cmap, transform = ccrs.PlateCarree(),
                               norm = colors.LogNorm(vmin=eroe100_linthresh, vmax=vmax),
                               cbar_kwargs = {"label": cbar_label}
                              )
    else:
        da.plot.pcolormesh(ax = ax, cmap = cmap, transform = ccrs.PlateCarree(),
                           vmin = vmin, vmax = vmax, levels = levels, 
                           cbar_kwargs = {"label": cbar_label}
                          )
    
    # Add the State Boundary Fence of Western Australia to plot.
    
    path_sbfwa = cf.get_path_for_sbfwa_def()
    if Path(path_sbfwa).exists() == False:
        cf.proc_sbfwa_def()
    gdf_sbfwa = gpd.read_file(path_sbfwa)
    ax.add_geometries(gdf_sbfwa.geometry, crs=ccrs.PlateCarree(), 
                      facecolor='none', edgecolor='k')
    
    ax.set_title(da.attrs["full_name"])
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    grid = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    grid.top_labels = False
    grid.right_labels = False
    
    # Display plot if there was no input ax.
    
    if ax_input == None:
        fig.tight_layout()
        plt.show()
        fig.clear()
        plt.close(fig)
        
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def create_quiver(da_u, da_v, extents=None, vmin=None, vmax=None, ax=None):
    
    """
    Create a plot with vector arrows representing vector values at each grid cell, 
    and backgroud shading indicating magnitude of vector values.
    
    Arguments:
        da_u (xarray.DataArray): DataArray containing zonal components of vector.
        da_v (xarray.DataArray): DataArray containing meridional components of vector.
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        vmin (float or int): Minimum of colourbar extents for a calc plot.
        vmax (float or int): Maximum of colourbar extents for a calc plot.
        ax (cartopy.GeoAxesSubplot): Figure axis to create plot on.
        
    Returns:
        A quiver plot (vectors) overlayed on a pcolormesh plot (magnitudes).
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid.
    
    assert ((str(type(da_u)) == "<class 'xarray.core.dataarray.DataArray'>") & 
            (da_u.dims == da_dims_valid)), \
        f"da_u must be an xarray.DataArray with da_u.dims == {da_dims_valid}"
    assert ((str(type(da_v)) == "<class 'xarray.core.dataarray.DataArray'>") & 
            (da_v.dims == da_dims_valid)), \
        f"da_v must be an xarray.DataArray with da_v.dims == {da_dims_valid}"
    cf.check_args(extents=extents, vmin=vmin, vmax=vmax, ax=ax)
    
    attrs_u = copy.deepcopy(da_u.attrs)
    attrs_u["abbreviation"] = attrs_u["abbreviation"].replace("_u", "")
    attrs_u["full_name"] = (attrs_u["full_name"]
                            .replace("Zonal Component of ", ""))
    
    attrs_v = copy.deepcopy(da_v.attrs)
    attrs_v["full_name"] = (attrs_v["full_name"]
                            .replace("Meridional Component of ", ""))
    
    assert attrs_u["full_name"] == attrs_v["full_name"], \
        ("da_u and da_v must be the zonal and meridional components " +
         "of the same variable")
    
    vector_test = (attrs_u["abbreviation"]
                   .split("(")[-1]
                   .split(")")[0]
                   .split("^")[0]
                   .split("$")[0]
                   .replace("U", "WV")
                   .lower())
    
    assert vector_test in cf.params_vector, \
        "da_u and da_v must be the components of a vector parameter"
    
    if extents:
        assert (isinstance(extents, list) & (len(extents) == 4) & 
                (extents[1] > extents[0]) & (extents[3] > extents[2])), \
            "extents must a 4 element list in [W, E, S, N] format or None"
    else:
        extents = []
        extents.append(da_u.longitude.min())
        extents.append(da_u.longitude.max())
        extents.append(da_u.latitude.min())
        extents.append(da_u.latitude.max())
    
    # Create copy of input ax to refer to when deciding whether to display plot.
    
    ax_input = ax
    
    # If no ax is specified (i.e. this function is not being used to create a subplot
    # within a figure for a different plot function), then create a figure.
    
    if ax == None:
        figrows = 1
        figcols = 1
        figwidth = figwidth_standard
        figheight = (figwidth * figrows/figcols * 
                     (extents[3]-extents[2])/(extents[1]-extents[0])
                    )
        fig, ax = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                               subplot_kw = {"projection": ccrs.PlateCarree()}
                              )
    
    # For study regions with very wide extents, coarsen the grid into a lower
    # resolution by averaging over values. This is to avoid too many quiver
    # arrows all overlapping each other and producing an unreadable plot.
    
    coarsen_window_size = math.ceil((extents[1]-extents[0]) / figwidth_standard)
    da_u = (da_u
            .coarsen(longitude = coarsen_window_size, boundary = "trim")
            .mean()
            .coarsen(latitude = coarsen_window_size, boundary = "trim")
            .mean()
           )
    da_v = (da_v
            .coarsen(longitude = coarsen_window_size, boundary = "trim")
            .mean()
            .coarsen(latitude = coarsen_window_size, boundary = "trim")
            .mean()
           )    
    da_mag = xr.DataArray(cf.get_magnitude(da_u, da_v), name = "mag")
    ds = xr.merge([da_u, da_v])   
    
    # Specify colourmap, colourbar extents and levels in colourbar depending
    # on category of parameter we are plotting.
    
    if attrs_u["full_name"].split(" ")[0] == "Difference":
        cmap = cmocean.cm.speed
    else:
        cmap = cmocean.cm.tempo
    da_mag_subset = da_mag.sel(longitude=slice(extents[0], extents[1]), 
                               latitude=slice(extents[3], extents[2]))
    vmin = vmin if vmin else da_mag_subset.min()
    vmax = vmax if vmax else da_mag_subset.max()
    cbar_label = "{} [{}]".format(attrs_u["abbreviation"], attrs_u["units"])
    scale = float(vmax) * quiver_scale_multiplier
    
    if vmin != None:
        if math.isnan(vmin):
            vmin = None
    if vmax != None:
        if math.isnan(vmax):
            vmax = None
    
    # Create plot
    
    ax.set_extent(extents=extents, crs=ccrs.PlateCarree())
    da_mag.plot.pcolormesh(ax = ax, cmap = cmap, transform = ccrs.PlateCarree(),
                           vmin = vmin, vmax = vmax, cbar_kwargs = {"label": cbar_label}
                          )
    ds.plot.quiver(x = "longitude", y = "latitude", ax = ax, 
                   u = da_u.name, v = da_v.name, vmin = vmin, vmax = vmax, 
                   scale = scale, headwidth = quiver_headwidth, 
                   transform = ccrs.PlateCarree(), add_guide = False
                  )
    
    # Add the State Boundary Fence of Western Australia to plot.
    
    path_sbfwa = cf.get_path_for_sbfwa_def()
    if Path(path_sbfwa).exists() == False:
        cf.proc_sbfwa_def()
    gdf_sbfwa = gpd.read_file(path_sbfwa)
    ax.add_geometries(gdf_sbfwa.geometry, crs=ccrs.PlateCarree(), 
                      facecolor='none', edgecolor='k')
    
    ax.set_title(attrs_u["full_name"])
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    grid = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    grid.top_labels = False
    grid.right_labels = False
    
    # Display plot if there was no input ax.
    
    if ax_input == None:
        fig.tight_layout()
        plt.show()
        fig.clear()
        plt.close(fig)
        
    cf.remove_handlers_if_directly_executed(func_1up)


# ## Main mid-level plotting functions

# In[ ]:


def create_individual_calc_plot(
    calc_func, region, period_start, period_end, period_months, arg_extra, 
    period_hours=None, glass_source_pref=None, var_or_dvar=None, extents=None,
    vmin=None, vmax=None, ax=None, cfv_data=None, output=False
):
    
    """
    Create a single plot displaying the results from a calc_func output file.
    
    Arguments:
        calc_func (function): Calculation function to use in analysis. Must be one of: 
            [calc_glass_mean_clim,
            calc_era5_mdp_clim_given_var_or_dvar,
            calc_era5_mdp_clim_stats_given_var_or_dvar,
            calc_era5_mean_clim_given_var_or_dvar,
            calc_era5_wsd_clim].
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_months (str or list): Months subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        arg_extra (str or int): Extra plotting argument used to specify which GLASS 
            parameter to plot, which hour for the mean diurnal profile of an ERA5 
            parameter to plot, which statistic of the mean diurnal profile to plot, 
            or which parameter of the wind speed distribution to plot. Must be one of:
            ["mlai", "mfapar", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
            16, 17, 18, 19, 20, 21, 22, 23, "hour_max", "hour_min", "max", "max_u", 
            "max_v", "min", "min_u", "min_v", "mean", "mean_u", "mean_v", "range", 
            "ws10_mean", "ws10_std", "c10", "k10", "ws100_mean", "ws100_std", "c100", 
            "k100", "eroe100", "tgcf100"].
        period_hours (str or list): Hours subset of period to perform calculation over.
            Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
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
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        vmin (float or int): Minimum of colourbar extents for a calc plot.
        vmax (float or int): Maximum of colourbar extents for a calc plot.
        ax (cartopy.GeoAxesSubplot): Figure axis to create plot on.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
            
    Returns (if output = True):
        ../data_final/glass_mean_clim/{plot_funcs_ver}_{cfv_used}_calc_
            {extents_used}_{period_start}_{period_end}_{period_months_str}_
            glass-mean_{glass_source}_{arg_extra}.png OR
        ../data_final/era5_mdp_clim_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_
            calc_{extents_used}_{period_start}_{period_end}_{period_months_str}_
            era5-mdp_{var_or_dvar}_{arg_extra}.png OR
        ../data_final/era5_mdp_clim_stats_given_var_or_dvar/{plot_funcs_ver}_
            {cfv_used}_calc_{extents_used}_{period_start}_{period_end}_
            {period_months_str}_era5-mdp_{var_or_dvar}_stats_{arg_extra}.png OR
        ../data_final/era5_mean_clim_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_
            calc_{extents_used}_{period_start}_{period_end}_{period_months_str}_
            era5-mean_{period_hours_str}_{var_or_dvar}_{arg_extra}.png OR
        ../data_final/era5_wsd_clim/{plot_funcs_ver}_{cfv_used}_calc_
            {extents_used}_{period_start}_{period_end}_{period_months_str}_
            era5-wsd_{period_hours_str}_{arg_extra}.png:
                Output PNG file in data_final folder for the calculation plot. 
                {plot_funcs_ver} is the version of the plot_funcs script being used. 
                {cfv_used} is the version of the calc_funcs script which outputted 
                the data used in making this plot. {extents_used} is a string 
                indicating the region, or the WESN coordinates for the extents argument 
                if this was specified. {period1_months_str} and {period2_months_str} 
                are strings representing the list of selected months to use as a 
                subset in each period.
                
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    calc_func_name = calc_func.__name__
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(calc_func_name, args_cur, args_cur_values)
    cf.check_args(calc_func=calc_func, region=region, period_start=period_start, 
                  period_end=period_end, period_months=period_months, 
                  arg_extra=arg_extra, period_hours=period_hours, 
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar, 
                  extents=extents, vmin=vmin, vmax=vmax,
                  ax=ax, cfv_data=cfv_data, output=output)
    
    # Obtain string representation for month and hour subsets.
    
    period_months_str = cf.get_period_months_str(period_months=period_months).upper()
    if period_hours:
        period_hours_str = cf.get_period_hours_str(period_hours=period_hours)
    
    # Create copy of input extents to later decide what to name output file.
    
    extents_input = copy.deepcopy(extents)
    
    # If no input extents are entered, use entire region defined in calc_funcs.
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    # Create copy of input ax to refer to when deciding whether to display plot.
    
    ax_input = ax
    
    # If no ax is specified (i.e. this function is not being used to create a subplot
    # within a figure for a different plot function), then create a figure.
    
    if ax == None:
        figrows = 1
        figcols = 1
        figwidth = figwidth_standard
        figheight = (figwidth * figrows/figcols * 
                     (extents[3]-extents[2])/(extents[1]-extents[0])
                    )
        fig, ax = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                               subplot_kw = {"projection": ccrs.PlateCarree()}
                              )
    
    # Obtain path to calc_func output.
    
    path_calc = cf.get_path_for_calc_func(
        calc_func_name=calc_func_name, region=region, period_start=period_start, 
        period_end=period_end, period_months=period_months, period_hours=period_hours,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar
    )
    
    # Use data outputted from an older version of the calc_funcs script. This is useful
    # for results which required computationally intensive processing. And can also be
    # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
    if cfv_data:
        path_calc = path_calc.replace(cf.calc_funcs_ver, cfv_data)
        if Path(path_calc).exists() == False:
            msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                       f"but could not find file: {path_calc}")
            logging.error(msg_cfv)
            print(msg_cfv)
            cf.remove_handlers_if_directly_executed(func_1up)
            raise Exception(msg_cfv)
    
    # Open calc_func output file if it exists, otherwise create the file.
    
    if Path(path_calc).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_calc}"
        logging.info(msg_open)
        print(msg_open)
    else:
        calc_func(region=region, period_start=period_start, period_end=period_end, 
                  period_months=period_months, period_hours=period_hours, 
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    
    ds_calc = xr.open_dataset(path_calc, engine = "netcdf4")
    
    # Open relevant DataArray's and create plots
    
    if calc_func_name == "calc_era5_mdp_clim_given_var_or_dvar":
        if var_or_dvar in cf.params_vector:
            da_u = ds_calc[var_or_dvar.replace("wv", "u")].sel(hour=arg_extra)
            da_v = ds_calc[var_or_dvar.replace("wv", "v")].sel(hour=arg_extra)
            create_quiver(da_u=da_u, da_v=da_v, extents=extents, 
                          vmin=vmin, vmax=vmax, ax=ax)
        else:
            da_calc = ds_calc[var_or_dvar].sel(hour=arg_extra)
            create_pcolormesh(da=da_calc, extents=extents, vmin=vmin, vmax=vmax, ax=ax)
    elif (var_or_dvar in cf.params_vector) & (arg_extra in ["max", "min", "mean"]):
        da_u = ds_calc[arg_extra + "_u"]
        da_v = ds_calc[arg_extra + "_v"]
        create_quiver(da_u=da_u, da_v=da_v, extents=extents, vmin=vmin, vmax=vmax, ax=ax)
    else:
        da_calc = ds_calc[arg_extra]
        create_pcolormesh(da=da_calc, extents=extents, vmin=vmin, vmax=vmax, ax=ax)
    
    # Rename plot title to include metadata.
    
    ax_title = ax.get_title()
    
    if calc_func_name == "calc_glass_mean_clim":
        source_str = ds_calc[arg_extra].attrs["source"].upper()
        ax.set_title("\n".join(wrap(f"{ax_title} [{period_start} to {period_end} " +
                                    f"(months={period_months_str}); source={source_str}]", 
                                    title_width)))
    elif calc_func_name == "calc_era5_mdp_clim_given_var_or_dvar":
        ax.set_title("\n".join(wrap(f"Hour={arg_extra} Value for {ax_title} " +
                                    f"[{period_start} to {period_end} " + 
                                    f"(months={period_months_str})]", title_width)))
    elif calc_func_name == "calc_era5_mdp_clim_stats_given_var_or_dvar":
        ax.set_title("\n".join(wrap(f"{ax_title} [{period_start} to {period_end} " + 
                                    f"(months={period_months_str})]", title_width)))
    else:
        ax.set_title("\n".join(wrap(f"{ax_title} [{period_start} to {period_end} " + 
                                    f"(months={period_months_str}, " +
                                    f"hours={period_hours_str})]", title_width)))
   
    # If no ax is specified (i.e. this function is not being used to create a subplot
    # within a figure for a different plot function), then edit figure layout.
    
    if ax_input == None:
        fig.tight_layout()
        
        if output == True:
            if cfv_data:
                cfv_used = copy.deepcopy(cfv_data)
            else:
                cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
            if extents_input:
                extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                                     extents[2], extents[3])
            else:
                extents_used = copy.deepcopy(region)
            
            path_output = (path_calc
                           .replace("data_processed", "data_final")
                           .replace(".nc", f"_{arg_extra}.png")
                           .replace(f"{cfv_used}_calc_{region}", 
                                    f"{plot_funcs_ver}_{cfv_used}_calc_{extents_used}")
                          )
            path_output_dir = "/".join(path_output.split("/")[:-1])
            Path(path_output_dir).mkdir(parents=True, exist_ok=True)
            path_output = path_output.replace(
                path_output_dir + "/", path_output_dir + "/" + f"{plot_funcs_ver}_")
            
            if Path(path_output).exists():
                msg_exist = ("WARNING: plot file already exists (and was " +
                             f"not overwritten): {path_output}")
                logging.warning(msg_exist)
                print(msg_exist)
            else:
                plt.savefig(path_output, metadata=get_plot_metadata(
                    time_exec, func_cur, args_cur, args_cur_values)
                           )
                msg_create = f"CREATED: plot file: {path_output}"
                logging.info(msg_create)
                print(msg_create)
                
        plt.show()
        fig.clear()
        plt.close(fig)
        
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def create_individual_diff_plot(
    calc_func, region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, arg_extra, 
    period1_hours=None, period2_hours=None, glass_source_pref=None, var_or_dvar=None,
    perc=False, mask_perc_quantile=mask_perc_quantile_default, extents=None, 
    vmin=None, vmax=None, ax=None, cfv_data=None, output=False
):
    
    """
    Create a single plot displaying the difference in results between two calc_func 
    output files (corresponding to difference in summary statistics over two periods).
    
    Arguments:
        calc_func (function): Calculation function to use in analysis. Must be one of: 
            [calc_glass_mean_clim,
            calc_era5_mdp_clim_given_var_or_dvar,
            calc_era5_mdp_clim_stats_given_var_or_dvar,
            calc_era5_mean_clim_given_var_or_dvar,
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        arg_extra (str or int): Extra plotting argument used to specify which GLASS 
            parameter to plot, which hour for the mean diurnal profile of an ERA5 
            parameter to plot, which statistic of the mean diurnal profile to plot, 
            or which parameter of the wind speed distribution to plot. Must be one of:
            ["mlai", "mfapar", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
            16, 17, 18, 19, 20, 21, 22, 23, "hour_max", "hour_min", "max", "max_u", 
            "max_v", "min", "min_u", "min_v", "mean", "mean_u", "mean_v", "range", 
            "ws10_mean", "ws10_std", "c10", "k10", "ws100_mean", "ws100_std", "c100", 
            "k100", "eroe100", "tgcf100"].
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
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
        perc (bool): Whether to plot the difference in values as a percentage of the
            value (magnitude if negative) in period1. This is used for the comp plots
            in the plot_funcs script. Must be one of: [True, False].
        mask_perc_quantile (int): If perc is True, specify the quantile of values
            (magnitude if negative) from period1 to mask for the difference plot.
            This is used because percentage differences may be particularly high
            for values which had a low magnitude as a base in period 1.
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        vmin (float or int): Minimum of colourbar extents for a calc plot.
        vmax (float or int): Maximum of colourbar extents for a calc plot.
        ax (cartopy.GeoAxesSubplot): Figure axis to create plot on.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
            
    Returns (if output = True):
        ../data_final/glass_mean_clim/{plot_funcs_ver}_{cfv_used}_diff_
            {extents_used}_{period1_start}_{period1_end}_{period2_start}_
            {period2_end}_{period1_months_str}_{period2_months_str}_glass-mean_
            {glass_source}_{arg_extra}_perc-{quantile_used}.png OR
        ../data_final/era5_mdp_clim_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_
            diff_{extents_used}_{period1_start}_{period1_end}_{period2_start}_
            {period2_end}_{period1_months_str}_{period2_months_str}_era5-mdp_
            {var_or_dvar}_{arg_extra}_perc-{quantile_used}.png OR
        ../data_final/era5_mdp_clim_stats_given_var_or_dvar/{plot_funcs_ver}_
            {cfv_used}_diff_{extents_used}_{period1_start}_{period1_end}_
            {period2_start}_{period2_end}_{period1_months_str}_{period2_months_str}_
            era5-mdp_{var_or_dvar}_stats_{arg_extra}_perc-{quantile_used}.png OR
        ../data_final/era5_mean_clim_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_
            diff_{extents_used}_{period1_start}_{period1_end}_{period2_start}_
            {period2_end}_{period1_months_str}_{period2_months_str}_era5-mean_
            {period1_hours_str}_{period2_hours_str}_{var_or_dvar}_{arg_extra}_
            perc-{quantile_used}.png OR
        ../data_final/era5_wsd_clim/{plot_funcs_ver}_{cfv_used}_diff_{extents_used}_
            {period1_start}_{period1_end}_{period2_start}_{period2_end}_
            {period1_months_str}_{period2_months_str}_era5-wsd_{period1_hours_str}_
            {period2_hours_str}_{arg_extra}_perc-{quantile_used}.png:
                Output PNG file in data_final folder for the difference plot. 
                {plot_funcs_ver} is the version of the plot_funcs script being used. 
                {cfv_used} is the version of the calc_funcs script which outputted 
                the data used in making this plot. {extents_used} is a string 
                indicating the region, or the WESN coordinates for the extents 
                argument if this was specified. {period1_months_str} and 
                {period2_months_str} are strings representing the list of 
                selected months to use as a subset in each period. 
                {quantile_used} is equal to mask_perc_quantile if perc = True, 
                otherwise it is set to None.
                
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    The mask applies only to scalar variables (i.e. quiver plots for vectors are not 
    masked by their magnitude).
                
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    calc_func_name = calc_func.__name__
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(calc_func_name, args_cur, args_cur_values)
    cf.check_args(calc_func=calc_func, region=region, period1_start=period1_start, 
                  period1_end=period1_end, period2_start=period2_start, 
                  period2_end=period2_end, period1_months=period1_months, 
                  period2_months=period2_months, period1_hours=period1_hours,
                  period2_hours=period2_hours, arg_extra=arg_extra, 
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar, 
                  perc=perc, mask_perc_quantile=mask_perc_quantile, extents=extents,
                  vmin=vmin, vmax=vmax, ax=ax, cfv_data=cfv_data, output=output)
    
    # Obtain string representation for month and hour subsets.
    
    period1_months_str = cf.get_period_months_str(period_months=period1_months).upper()
    period2_months_str = cf.get_period_months_str(period_months=period2_months).upper()
    if period1_hours:
        period1_hours_str = cf.get_period_hours_str(period_hours=period1_hours)
    if period2_hours:
        period2_hours_str = cf.get_period_hours_str(period_hours=period2_hours)
    
    # Create copy of input extents to later decide what to name output file.
    
    extents_input = copy.deepcopy(extents)
    
    # If no input extents are entered, use entire region defined in calc_funcs.
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    # Create copy of input ax to refer to when deciding whether to display plot.
    
    ax_input = ax
    
    # If no ax is specified (i.e. this function is not being used to create a subplot
    # within a figure for a different plot function), then create a figure.
    
    if ax == None:
        figrows = 1
        figcols = 1
        figwidth = figwidth_standard
        figheight = (figwidth * figrows/figcols * 
                     (extents[3]-extents[2])/(extents[1]-extents[0])
                    )
        fig, ax = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                               subplot_kw = {"projection": ccrs.PlateCarree()}
                              )
    
    # Obtain path to calc_func and calc_diff outputs.
    
    path_period1 = cf.get_path_for_calc_func(
        calc_func_name=calc_func_name, region=region, period_start=period1_start, 
        period_end=period1_end, period_months=period1_months, period_hours=period1_hours,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar
    )
    path_period2 = cf.get_path_for_calc_func(
        calc_func_name=calc_func_name, region=region, period_start=period2_start, 
        period_end=period2_end, period_months=period2_months, period_hours=period2_hours,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar
    )
    path_diff = cf.get_path_for_calc_diff(
        calc_func_name=calc_func_name, region=region, 
        period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end,
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=period1_hours, period2_hours=period2_hours,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar
    )
    
    # Use data outputted from an older version of the calc_funcs script. This is useful
    # for results which required computationally intensive processing. And can also be
    # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
    if cfv_data:
        for path in [path_period1, path_period2, path_diff]:
            path = path.replace(cf.calc_funcs_ver, cfv_data)
            if Path(path).exists() == False:
                msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                           f"but could not find file: {path}")
                logging.error(msg_cfv)
                print(msg_cfv)
                cf.remove_handlers_if_directly_executed(func_1up)
                raise Exception(msg_cfv)
    
    # Open calc_func and calc_diff output files if they exists, 
    # otherwise create the files.
    
    if Path(path_period1).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_period1}"
        logging.info(msg_open)
        print(msg_open)
    else:
        calc_func(region=region, period_start=period1_start, period_end=period1_end, 
                  period_months=period1_months, period_hours=period1_hours,
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    
    ds_period1 = xr.open_dataset(path_period1, engine = "netcdf4")
    
    if Path(path_period2).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_period2}"
        logging.info(msg_open)
        print(msg_open)
    else:
        calc_func(region=region, period_start=period2_start, period_end=period2_end, 
                  period_months=period2_months, period_hours=period2_hours, 
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    
    ds_period2 = xr.open_dataset(path_period2, engine = "netcdf4")
    
    if Path(path_diff).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_diff}"
        logging.info(msg_open)
        print(msg_open)
    else:
        cf.calc_diff(calc_func=calc_func, region=region, 
                     period1_start=period1_start, period1_end=period1_end, 
                     period2_start=period2_start, period2_end=period2_end, 
                     period1_months=period1_months, period2_months=period2_months,
                     period1_hours=period1_hours, period2_hours=period2_hours,
                     glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar)
    
    ds_diff = xr.open_dataset(path_diff, engine = "netcdf4")
    
    # Open relevant DataArray's and create plots
    
    if calc_func_name == "calc_era5_mdp_clim_given_var_or_dvar":
        if var_or_dvar in cf.params_vector:
            da_u = ds_diff[var_or_dvar.replace("wv", "u")].sel(hour=arg_extra)
            da_v = ds_diff[var_or_dvar.replace("wv", "v")].sel(hour=arg_extra)
            create_quiver(da_u=da_u, da_v=da_v, extents=extents, 
                          vmin=vmin, vmax=vmax, ax=ax)
        else:
            da_diff = ds_diff[var_or_dvar].sel(hour=arg_extra)
            if perc == True:
                da_period1 = ds_period1[var_or_dvar].sel(hour=arg_extra)
                da_period1_mag = cf.get_magnitude(da_period1, da_period1)
                da_period1_mag = da_period1_mag.where(
                    da_period1_mag > da_period1_mag.quantile(mask_perc_quantile / 100))
                da_diff_perc = da_diff / da_period1_mag * 100
                da_diff_perc.attrs = copy.deepcopy(da_diff.attrs)
                da_diff_perc.attrs["units"] = "\%"
                create_pcolormesh(da=da_diff_perc, extents=extents, 
                                  vmin=vmin, vmax=vmax, ax=ax)
            else:
                create_pcolormesh(da=da_diff, extents=extents, 
                                  vmin=vmin, vmax=vmax, ax=ax)
    elif (var_or_dvar in cf.params_vector) & (arg_extra in ["max", "min", "mean"]):
        da_u = ds_diff[arg_extra + "_u"]
        da_v = ds_diff[arg_extra + "_v"]
        create_quiver(da_u=da_u, da_v=da_v, extents=extents, vmin=vmin, vmax=vmax, ax=ax)
    else:
        da_diff = ds_diff[arg_extra]
        if perc == True:
            da_period1 = ds_period1[arg_extra]
            da_period1_mag = cf.get_magnitude(da_period1, da_period1)
            da_period1_mag = da_period1_mag.where(
                da_period1_mag > da_period1_mag.quantile(mask_perc_quantile / 100))
            da_diff_perc = da_diff / da_period1_mag * 100
            da_diff_perc.attrs = copy.deepcopy(da_diff.attrs)
            da_diff_perc.attrs["units"] = "\% of period 1"
            create_pcolormesh(da=da_diff_perc, extents=extents, 
                              vmin=vmin, vmax=vmax, ax=ax)
        else:
            create_pcolormesh(da=da_diff, extents=extents, 
                              vmin=vmin, vmax=vmax, ax=ax)
    
    # Rename plot title to include metadata.
    
    ax_title = ax.get_title()
    
    if calc_func_name == "calc_glass_mean_clim":
        source_str = ds_diff[arg_extra].attrs["source"].upper()
        ax.set_title("\n".join(wrap(f"{ax_title} [{period2_start} to {period2_end} " +
                                    f"(months={period2_months_str}) minus " +
                                    f"{period1_start} to {period1_end} (months=" +
                                    f"{period1_months_str}); source={source_str}]", 
                                    title_width)))
    elif calc_func_name == "calc_era5_mdp_clim_given_var_or_dvar":
        ax_title = ax_title.replace("Difference in ", 
                                    f"Difference in Hour={arg_extra} Value for ")
        ax.set_title("\n".join(wrap(f"{ax_title} [{period2_start} to {period2_end} " +
                                    f"(months={period2_months_str}) minus " +
                                    f"{period1_start} to {period1_end} (months=" +
                                    f"{period1_months_str})]", title_width)))
    elif calc_func_name == "calc_era5_mdp_clim_stats_given_var_or_dvar":
        ax.set_title("\n".join(wrap(f"{ax_title} [{period2_start} to {period2_end} " +
                                    f"(months={period2_months_str}) minus " +
                                    f"{period1_start} to {period1_end} (months=" +
                                    f"{period1_months_str})]", title_width)))
    else:
        ax.set_title("\n".join(wrap(f"{ax_title} [{period2_start} to {period2_end} " +
                                    f"(months={period2_months_str}, " + 
                                    f"hours={period2_hours_str}) minus "
                                    f"{period1_start} to {period1_end} " +
                                    f"(months={period1_months_str}, " +
                                    f"hours={period1_hours_str})]", title_width)))
    
    # If no ax is specified (i.e. this function is not being used to create a subplot
    # within a figure for a different plot function), then edit figure layout.
    
    if ax_input == None:
        fig.tight_layout()
        
        if output == True:
            if cfv_data:
                cfv_used = copy.deepcopy(cfv_data)
            else:
                cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
            if extents_input:
                extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                                     extents[2], extents[3])
            else:
                extents_used = copy.deepcopy(region)
            
            if perc == True:
                quantile_used = None
            else:
                quantile_used = copy.deepcopy(mask_perc_quantile)
            
            path_output = (path_diff
                           .replace("data_processed", "data_final")
                           .replace(".nc", f"_{arg_extra}_perc-{quantile_used}.png")
                           .replace(f"{cfv_used}_diff_{region}", 
                                    f"{plot_funcs_ver}_{cfv_used}_diff_{extents_used}")
                          )
            path_output_dir = "/".join(path_output.split("/")[:-1])
            Path(path_output_dir).mkdir(parents=True, exist_ok=True)
            
            if Path(path_output).exists():
                msg_exist = ("WARNING: plot file already exists (and was " +
                             f"not overwritten): {path_output}")
                logging.warning(msg_exist)
                print(msg_exist)
            else:
                plt.savefig(path_output, metadata=get_plot_metadata(
                    time_exec, func_cur, args_cur, args_cur_values)
                           )
                msg_create = f"CREATED: plot file: {path_output}"
                logging.info(msg_create)
                print(msg_create)
                
        plt.show()
        fig.clear()
        plt.close(fig)
        
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def create_individual_comp_plot(
    calc_func, region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, arg_extra, period1_hours=None, period2_hours=None, 
    glass_source_pref=None, var_or_dvar=None, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, 
    vmin_periods=None, vmax_periods=None, vmin_diff=None, vmax_diff=None, 
    ax_period1=None, ax_period2=None, ax_diff=None, cfv_data=None, output=False
):
    
    """
    Create a figure consisting of a row with 3 subplots displaying results from the 
    calc_func output files over each specified periods, as well as the difference 
    in results between these two periods.
    
    Arguments:
        calc_func (function): Calculation function to use in analysis. Must be one of: 
            [calc_glass_mean_clim,
            calc_era5_mdp_clim_given_var_or_dvar,
            calc_era5_mdp_clim_stats_given_var_or_dvar,
            calc_era5_mean_clim_given_var_or_dvar,
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        arg_extra (str or int): Extra plotting argument used to specify which GLASS 
            parameter to plot, which hour for the mean diurnal profile of an ERA5 
            parameter to plot, which statistic of the mean diurnal profile to plot, 
            or which parameter of the wind speed distribution to plot. Must be one of:
            ["mlai", "mfapar", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
            16, 17, 18, 19, 20, 21, 22, 23, "hour_max", "hour_min", "max", "max_u", 
            "max_v", "min", "min_u", "min_v", "mean", "mean_u", "mean_v", "range", 
            "ws10_mean", "ws10_std", "c10", "k10", "ws100_mean", "ws100_std", "c100", 
            "k100", "eroe100", "tgcf100"].
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
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
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        vmin (float or int): Minimum of colourbar extents for a calc plot.
        vmax (float or int): Maximum of colourbar extents for a calc plot.
        ax (cartopy.GeoAxesSubplot): Figure axis to create plot on.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
            
    Returns (if output = True):
        ../data_final/glass_mean_clim/{plot_funcs_ver}_{cfv_used}_diff_
            {extents_used}_{period1_start}_{period1_end}_{period2_start}_
            {period2_end}_{period1_months_str}_{period2_months_str}_glass-mean_
            {glass_source}_{arg_extra}_perc-{quantile_used}_mask1-{mask_period1}_
            mask2-{mask_period2}.png OR
        ../data_final/era5_mdp_clim_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_
            diff_{extents_used}_{period1_start}_{period1_end}_{period2_start}_
            {period2_end}_{period1_months_str}_{period2_months_str}_era5-mdp_
            {var_or_dvar}_{arg_extra}_perc-{quantile_used}_mask1-{mask_period1}_
            mask2-{mask_period2}.png OR
        ../data_final/era5_mdp_clim_stats_given_var_or_dvar/{plot_funcs_ver}_
            {cfv_used}_diff_{extents_used}_{period1_start}_{period1_end}_
            {period2_start}_{period2_end}_{period1_months_str}_{period2_months_str}_
            era5-mdp_{var_or_dvar}_stats_{arg_extra}_perc-{quantile_used}_
            mask1-{mask_period1}_mask2-{mask_period2}.png OR
        ../data_final/era5_mean_clim_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_
            diff_{extents_used}_{period1_start}_{period1_end}_{period2_start}_
            {period2_end}_{period1_months_str}_{period2_months_str}_era5-mean_
            {period1_hours_str}_{period2_hours_str}_{var_or_dvar}_{arg_extra}_
            perc-{quantile_used}_mask1-{mask_period1}_mask2-{mask_period2}.png OR
        ../data_final/era5_wsd_clim/{plot_funcs_ver}_{cfv_used}_diff_{extents_used}_
            {period1_start}_{period1_end}_{period2_start}_{period2_end}_
            {period1_months_str}_{period2_months_str}_era5-wsd_{period1_hours_str}_
            {period2_hours_str}_{arg_extra}_perc-{quantile_used}_mask1-{mask_period1}_
            mask2-{mask_period2}.png:
                Output PNG file in data_final folder for the difference plot. 
                {plot_funcs_ver} is the version of the plot_funcs script being used. 
                {cfv_used} is the version of the calc_funcs script which outputted 
                the data used in making this plot. {extents_used} is a string 
                indicating the region, or the WESN coordinates for the extents 
                argument if this was specified. {period1_months_str} and 
                {period2_months_str} are strings representing the list of 
                selected months to use as a subset in each period. 
                {quantile_used} is equal to mask_perc_quantile if perc = True, 
                otherwise it is set to None.
                
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    If mask_period1 = "pos" or "neg" is specified, all grid cells where var_or_dvar was 
    positive or negative respectively in period1 will be masked, but only if var_or_dvar 
    was a parameter which could take both positive and negative values in the first place 
    (for positive-only parameters, this mask argument is ignored). Similar comments
    apply for mask_period2 but for values in period2. All 3 masks apply only to scalar
    variables (i.e. quiver plots for vectors are not masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    calc_func_name = calc_func.__name__
    
    # Assert that input arguments are valid
    
    axes_input = [ax_period1, ax_period2, ax_diff]
    assert (all(ax_input == None for ax_input in axes_input) | 
            all(ax_input != None for ax_input in axes_input)), \
        "ax_period1, ax_period2, ax_diff must be either all None or all not None"
    cf.check_args_for_none(calc_func_name, args_cur, args_cur_values)
    cf.check_args(calc_func=calc_func, region=region, period1_start=period1_start, 
                  period1_end=period1_end, period2_start=period2_start, 
                  period2_end=period2_end, period1_months=period1_months, 
                  period2_months=period2_months, period1_hours=period1_hours, 
                  period2_hours=period2_hours, arg_extra=arg_extra, 
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar, 
                  perc=perc, mask_perc_quantile=mask_perc_quantile, 
                  mask_period1=mask_period1, mask_period2=mask_period2, extents=extents, 
                  vmin_periods=vmin_periods, vmax_periods=vmax_periods, 
                  vmin_diff=vmin_diff, vmax_diff=vmax_diff,
                  ax_period1=ax_period1, ax_period2=ax_period2, ax_diff=ax_diff, 
                  cfv_data=cfv_data, output=output)
    
    # Create copy of input extents to later decide what to name output file.
    
    extents_input = copy.deepcopy(extents)
    
    # If no input extents are entered, use entire region defined in calc_funcs.
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    # Make sure colourbar limits for subplots over each period are standardised.
    
    if (vmin_periods == None) & (vmax_periods == None):
        vmin_periods, vmax_periods = get_common_cbar_limits(
            calc_func=calc_func, region=region, period1_start=period1_start, 
            period1_end=period1_end, period2_start=period2_start, 
            period2_end=period2_end, period1_months=period1_months, 
            period2_months=period2_months, period1_hours=period1_hours, 
            period2_hours=period2_hours, arg_extra=arg_extra, 
            glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar, 
            extents=extents, cfv_data=cfv_data
        )
    
    # Create copy of input ax_diff to refer to when deciding whether to display plot.
    
    ax_diff_input = ax_diff
    
    # If no ax_diff is specified (i.e. this function is not being used to create a 
    # subplots within a figure for a different plot function), then create a figure.
    
    if ax_diff == None:
        figrows = 1
        figcols = 3
        figwidth = figwidth_standard * 2
        figheight = (figwidth * figrows/figcols * 
                     (extents[3]-extents[2])/(extents[1]-extents[0])
                    )
        fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                                 subplot_kw = {"projection": ccrs.PlateCarree()}
                                )
        ax_period1 = axes[0]
        ax_period2 = axes[1]
        ax_diff = axes[2]
    
    # Obtain path to calc_diff output (this will later be edited to obtain the
    # output path for this comparison plot).
    
    path_diff = cf.get_path_for_calc_diff(
        calc_func_name=calc_func_name, region=region, 
        period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end,
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=period1_hours, period2_hours=period2_hours,
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar
    )
    
    # Use data outputted from an older version of the calc_funcs script. This is useful
    # for results which required computationally intensive processing. And can also be
    # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
    if cfv_data:
        path_diff = path_diff.replace(cf.calc_funcs_ver, cfv_data)
        if Path(path_diff).exists() == False:
            msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                       f"but could not find file: {path_diff}")
            logging.error(msg_cfv)
            print(msg_cfv)
            cf.remove_handlers_if_directly_executed(func_1up)
            raise Exception(msg_cfv)
    
    # Create plots.
    
    create_individual_calc_plot(
        calc_func=calc_func, region=region, period_start=period1_start, 
        period_end=period1_end, period_months=period1_months, 
        period_hours=period1_hours, arg_extra=arg_extra, 
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar, extents=extents, 
        vmin=vmin_periods, vmax=vmax_periods, ax=ax_period1, cfv_data=cfv_data, 
        output=False
    )
    create_individual_calc_plot(
        calc_func=calc_func, region=region, period_start=period2_start, 
        period_end=period2_end, period_months=period2_months, 
        period_hours=period2_hours, arg_extra=arg_extra, 
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar, extents=extents, 
        vmin=vmin_periods, vmax=vmax_periods, ax=ax_period2, cfv_data=cfv_data,
        output=False
    )
    create_individual_diff_plot(
        calc_func=calc_func, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=period1_hours, period2_hours=period2_hours, arg_extra=arg_extra, 
        glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar, 
        perc=perc, mask_perc_quantile=mask_perc_quantile, extents=extents, 
        vmin=vmin_diff, vmax=vmax_diff, ax=ax_diff, cfv_data=cfv_data, output=False
    )
    
    # If no ax_diff is specified (i.e. this function is not being used to create 
    # subplots within a figure for a different plot function), then edit figure layout.
    
    if ax_diff_input == None:
        
        # for idx in range(0, 2+1):
        #     ax_title = axes[idx].get_title()
        #     axes[idx].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
        
        fig.tight_layout()
        
        if output == True:
            if cfv_data:
                cfv_used = copy.deepcopy(cfv_data)
            else:
                cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
            if extents_input:
                extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                                     extents[2], extents[3])
            else:
                extents_used = copy.deepcopy(region)
            
            if perc == True:
                quantile_used = None
            else:
                quantile_used = copy.deepcopy(mask_perc_quantile)
            
            path_output = (path_diff
                           .replace("data_processed", "data_final")
                           .replace(".nc", f"_{arg_extra}_perc-{quantile_used}_" +
                                    f"mask1-{mask_period1}_mask2-{mask_period2}.png")
                           .replace(f"{cfv_used}_diff_{region}", 
                                    f"{plot_funcs_ver}_{cfv_used}_comp_{extents_used}")
                          )
            path_output_dir = "/".join(path_output.split("/")[:-1])
            Path(path_output_dir).mkdir(parents=True, exist_ok=True)
            
            if Path(path_output).exists():
                msg_exist = ("WARNING: plot file already exists (and was " +
                             f"not overwritten): {path_output}")
                logging.warning(msg_exist)
                print(msg_exist)
            else:
                plt.savefig(path_output, metadata=get_plot_metadata(
                    time_exec, func_cur, args_cur, args_cur_values)
                           )
                msg_create = f"CREATED: plot file: {path_output}"
                logging.info(msg_create)
                print(msg_create)
                
        plt.show()
        fig.clear()
        plt.close(fig)
        
    cf.remove_handlers_if_directly_executed(func_1up)


# ## Extra mid-level plotting functions

# In[ ]:


def create_orog_static_plot(param_orog, region=None, extents=None, vmin=None, vmax=None, 
                            ax=None, cfv_data=None, output=False):
    
    """
    Create a single plot displaying the land surface elevation or slope of 
    sub-gridscale orography over land using ERA5 data.
    
    Arguments:
        param_orog (str): ERA5 orographic parameter to create plot for. Must be one 
            of: ["lse", "ssgo"].
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        vmin (float or int): Minimum of colourbar extents for a calc plot.
        vmax (float or int): Maximum of colourbar extents for a calc plot.
        ax (cartopy.GeoAxesSubplot): Figure axis to create plot on.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
                        
    Returns (if output = True):
        ../data_final/era5_orog/{plot_funcs_ver}_{cfv_used}_calc_{extents_used}_
        era5-orog_{param_orog}.png:
            Output PNG file in data_final folder displaying the land surface elevation 
            or slope of sub-gridscale orography over land. {cfv_used} is the version of 
            the calc_funcs script which outputted the data used in making this plot. 
            {extents_used} is a string indicating the WESN coordinates for the extents 
            argument if this was specified, otherwise it is set equal to region if
            region was specified but not extents. If neither region or extents was
            specified, then {extents_used} is set to "global" instead.
    
    ERA5 land-sea mask data containing the fraction of each grid cell which is land is 
    then used to mask values above sea cover (by selecting only the values where this 
    fraction was greater than 0.5). The computation uses data from the data_raw folder, 
    then outputs intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(param_orog=param_orog, region=region, extents=extents, vmin=vmin,
                  vmax=vmax, ax=ax, cfv_data=cfv_data, output=output)
    
    # Create copy of input extents to later decide what to name output file.
    
    extents_input = copy.deepcopy(extents)
        
    # If no input extents are entered, use input region if not none, otherwise
    # set extents and region to global values.
        
    if extents == None:
        if region:
            extents = copy.deepcopy(cf.regions[region]["extents"])
        else:
            extents = [-180, 180, -90, 90]
            region = "global"
    
    # Create copy of input ax to refer to when deciding whether to display plot.
    
    ax_input = ax
    
    # If no ax is specified (i.e. this function is not being used to create a 
    # subplots within a figure for a different plot function), then create a figure.
    
    if ax == None:
        figrows = 1
        figcols = 1
        figwidth = figwidth_standard
        figheight = (figwidth * figrows/figcols * 
                     (extents[3]-extents[2])/(extents[1]-extents[0])
                    )
        fig, ax = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                               subplot_kw = {"projection": ccrs.PlateCarree()}
                              )
    
    # Obtain path to calc_era5_orog output (this will later be edited to obtain the
    # output path for this plot).
    
    path_orog = cf.get_path_for_era5_orog()
    
    # Use data outputted from an older version of the calc_funcs script. This is useful
    # for results which required computationally intensive processing. And can also be
    # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
    if cfv_data:
        path_orog = path_orog.replace(cf.calc_funcs_ver, cfv_data)
        if Path(path_orog).exists() == False:
            msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                       f"but could not find file: {path_orog}")
            logging.error(msg_cfv)
            print(msg_cfv)
            cf.remove_handlers_if_directly_executed(func_1up)
            raise Exception(msg_cfv)
    
    # Open calc_era5_orog output file if it exists, otherwise create the file.
    
    if Path(path_orog).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_orog}"
        logging.info(msg_open)
        print(msg_open)
    else:
        cf.calc_era5_orog()
    
    da_orog = xr.open_dataset(path_orog, engine = "netcdf4")[param_orog]
    
    # Create plot
    
    create_pcolormesh(da_orog, extents=extents, vmin=vmin, vmax=vmax, ax=ax)
    
    # If no ax is specified (i.e. this function is not being used to create a subplot
    # within a figure for a different plot function), then edit figure layout.
    
    if ax_input == None:
        fig.tight_layout()
        
        if output == True:
            if cfv_data:
                cfv_used = copy.deepcopy(cfv_data)
            else:
                cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
            if extents_input:
                extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                                     extents[2], extents[3])
            else:
                extents_used = copy.deepcopy(region)
                
            path_output = (path_orog
                           .replace("data_processed", "data_final")
                           .replace(".nc", f"_{param_orog}.png")
                           .replace(f"{cfv_used}_calc_global", 
                                    f"{plot_funcs_ver}_{cfv_used}_calc_{extents_used}")
                          )
            path_output_dir = "/".join(path_output.split("/")[:-1])
            Path(path_output_dir).mkdir(parents=True, exist_ok=True)
            
            if Path(path_output).exists():
                msg_exist = ("WARNING: plot file already exists (and was " +
                             f"not overwritten): {path_output}")
                logging.warning(msg_exist)
                print(msg_exist)
            else:
                plt.savefig(path_output, metadata=get_plot_metadata(
                    time_exec, func_cur, args_cur, args_cur_values)
                           )
                msg_create = f"CREATED: plot file: {path_output}"
                logging.info(msg_create)
                print(msg_create)
                
        plt.show()
        fig.clear()
        plt.close(fig)
        
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def create_glass_rolling_plot(region, year_start, year_end, period_months, window_size, 
                              param_glass_mean, glass_source_pref, extents=None, 
                              vmin=None, vmax=None, cfv_data=None, output=False):
    
    """
    Create a series of plots displaying the rolling average of the annual difference in 
    mean leaf area index (MLAI) and mean fraction of absorbed photosynthetically active 
    radiation (MFAPAR) using GLASS data.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of ["ca", "sa", "wa"].
        year_start (int): Earliest year to compute the rolling average for.
        year_end (int): Latest year to compute the rolling average for.
        period_months (str or list): Months subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        window_size (int): Rolling window size (in years) to compute average for.
            Must be an odd number and greater than or equal to 3.
        param_glass_mean (str): GLASS mean parameter to create plot for. Must be one 
            of: ["mlai", "mfapar"].
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        vmin (float or int): Minimum of colourbar extents for a calc plot.
        vmax (float or int): Maximum of colourbar extents for a calc plot.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
                        
    Returns:
        ../data_final/glass_rolling_avg_of_annual_diff/{plot_funcs_ver}_{cfv_used}_calc_
        {extents_used}_{year_start}_{year_end}_{period_months_str}_{window_size}-year_
        glass-rolling-diff_pref-{glass_source_pref}_{param_glass_mean}.png:
            Output PNG file in data_final folder displaying the rolling average
            of the annual difference in MLAI or MFAPAR. {calc_funcs_ver} is the version 
            of the calc_funcs script being used. {period_months_str} is a string 
            representing the list of selected months to use as a subset.
    
    For each grid cell, plot the rolling average of the annual difference in 
    MLAI and MFAPAR, and only using a subset of data within this period (if a 
    period_months not "all" is specified). These rolling averages are computed for each 
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
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid and get string
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, year_start=year_start, year_end=year_end, 
                  period_months=period_months, window_size=window_size, 
                  param_glass_mean=param_glass_mean, 
                  glass_source_pref=glass_source_pref, extents=extents, 
                  vmin=vmin, vmax=vmax, cfv_data=cfv_data, output=output)
    
    # Obtain string representation for month subset.
    
    months_str = cf.get_period_months_str(period_months=period_months).upper()
    
    # Create copy of input extents to later decide what to name output file.
    
    extents_input = copy.deepcopy(extents)
    
    # If no input extents are entered, use entire region defined in calc_funcs.
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    # Obtain path to calc_glass_rolling_avg_of_annual_diff output (this will 
    # later be edited to obtain the output path for this plot).
    
    path_roll = cf.get_path_for_calc_glass_rolling(
        region=region, year_start=year_start, year_end=year_end, 
        period_months=period_months, window_size=window_size, 
        glass_source_pref=glass_source_pref)
    
    # Use data outputted from an older version of the calc_funcs script. This is useful
    # for results which required computationally intensive processing. And can also be
    # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
    if cfv_data:
        path_roll = path_roll.replace(cf.calc_funcs_ver, cfv_data)
        if Path(path_roll).exists() == False:
            msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                       f"but could not find file: {path_roll}")
            logging.error(msg_cfv)
            print(msg_cfv)
            cf.remove_handlers_if_directly_executed(func_1up)
            raise Exception(msg_cfv)
    
    # Open intermediate output data file if it exists, otherwise create the file.
    
    if Path(path_roll).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_roll}"
        logging.info(msg_open)
        print(msg_open)
    else:
        cf.calc_glass_rolling_avg_of_annual_diff(
            region=region, year_start=year_start, year_end=year_end, 
            period_months=period_months, window_size=window_size, 
            glass_source_pref=glass_source_pref)
    da_roll = xr.open_dataset(path_roll, engine = "netcdf4")[param_glass_mean]
    
    # Create plots
    
    if cf.priority == "speed":
        da_roll = da_roll.persist()
    
    if (vmin == None) & (vmax == None):
        min_of_da_roll = float(da_roll.min())
        max_of_da_roll = float(da_roll.max())
        vmin = min(-abs(min_of_da_roll), -abs(max_of_da_roll))
        vmax = -vmin    
    
    figcols = 4
    figrows = math.ceil(len(da_roll.year) / figcols)
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    for idx, year in enumerate(da_roll.year.data):
        row = math.floor(idx / figcols)
        col = idx % figcols
        da_roll_year = da_roll.sel(year=year)
        create_pcolormesh(da=da_roll_year, extents=extents, 
                          vmin=vmin, vmax=vmax, ax=axes[row, col])
        ax_title = (("{}-Year Rolling Average of Annual Difference in {} " +
                     "(centred upon {}; months={}; source={})")
                    .format(window_size, param_glass_mean.upper(), year, 
                            months_str, da_roll.attrs[str(year)].upper()))
        axes[row, col].set_title("\n".join(wrap(ax_title, title_width)))
    
    fig.tight_layout()
        
    if output == True:
        if cfv_data:
            cfv_used = copy.deepcopy(cfv_data)
        else:
            cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
        if extents_input:
            extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                                 extents[2], extents[3])
        else:
            extents_used = copy.deepcopy(region)
        
        path_output = (path_roll
                       .replace("data_processed", "data_final")
                       .replace(".nc", f"_{param_glass_mean}.png")
                       .replace(f"{cfv_used}_calc_{region}", 
                                f"{plot_funcs_ver}_{cfv_used}_calc_{extents_used}")
                      )
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
            
        if Path(path_output).exists():
            msg_exist = ("WARNING: plot file already exists (and was " +
                         f"not overwritten): {path_output}")
            logging.warning(msg_exist)
            print(msg_exist)
        else:
            plt.savefig(path_output, metadata=get_plot_metadata(
                time_exec, func_cur, args_cur, args_cur_values)
                       )
            msg_create = f"CREATED: plot file: {path_output}"
            logging.info(msg_create)
            print(msg_create)
                
    plt.show()
    fig.clear()
    plt.close(fig)
        
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def create_climate_indices_plot(
    year_start, year_end, window_size, period1_mid=None, period2_mid=None, 
    month1_mark=None, month2_mark=None, cfv_data=None, output=False
):
    
    """
    Create a series of 8 plots (one on each row) displaying the monthly values 
    and rolling average for 8 separate NOAA climate indices, as well as any 
    El Nino-Southern Oscillation (ENSO) or Indian Ocean Dipole (IOD) events.
    
    Arguments:
        year_start (int): Earliest year to compute the rolling average for.
        year_end (int): Latest year to compute the rolling average for.
        window_size (int): Rolling window size (in years) to compute average for.
            Must be an odd number and greater than or equal to 3.
        period1_mid (str): Month in middle of first period to plot rolling average 
            of climate indices for. Must be of form "%b-%Y" eg. "Jul-1990".
        period2_mid (str): Month in middle of second period to plot rolling average 
            of climate indices for. Must be of form "%b-%Y" eg. "Jul-1990".
        month1_mark (str): First month to plot monthly value for in the climate 
            indices plot. Must be of form "%b-%Y" eg. "Jul-1990".
        month2_mark (str): Second month to plot monthly value for in the climate 
            indices plot. Must be of form "%b-%Y" eg. "Jul-1990".
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
                        
    Returns:
        ../data_final/noaa_ind/{plot_funcs_ver}_{cfv_used}_proc_global_noaa-ind_
        {year_start}_{year_end}_{window_size}_mid1-{period1_mid_str}_
        mid2-{period2_mid_str}_mark1-{month1_mark_str}_mark2-{month2_mark_str}.png:
            Output PNG file in data_final folder containing the processed
            NOAA climate indices. {calc_funcs_ver} is the version of the 
            calc_funcs script being used. {period1_mid_str}, {period2_mid_str},
            {month1_mark_str} and {month2_mark_str} are the string representations
            of period1_mid, period2_mid, month1_mark and month2_mark respectively.
            So None inputs will be represented as a string.
    
    If period_mid is specified, then the rolling window centred at that point is
    displayed along with the rolling average value. If month_mark is specified,
    then the monthly value at that point is displayed. The 8 indices displayed are
    the Atlantic Multidecadal Oscillation Index (AMOI), Pacific Decadal Oscillation 
    Index (PDOI), Oceanic Nino Index (ONI), Dipole Mode Index (DMI), Antarctic 
    Oscillation Index (AAOI), Arctic Oscillation Index (AOI), North Atlantic 
    Oscillation Index (NAOI), Eastern Pacific Oscillation Index (EPOI). ENSO and
    IOD events are plotted according to Japanese Meteorological Agency (JMA) data
    (note that different meteorological agencies use different definitions).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(year_start=year_start, year_end=year_end, window_size=window_size, 
                  period1_mid=period1_mid, period2_mid=period2_mid, 
                  month1_mark=month1_mark, month2_mark=month2_mark,
                  cfv_data=cfv_data, output=output)
    
    # Obtain path to proc_noaa_ind output (this will later be edited to obtain
    # the output path for this plot).
    
    path_noaa = cf.get_path_for_noaa_ind()
    
    # Use data outputted from an older version of the calc_funcs script. This is useful
    # for results which required computationally intensive processing. And can also be
    # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
    if cfv_data:
        path_noaa = path_noaa.replace(cf.calc_funcs_ver, cfv_data)
        if Path(path_noaa).exists() == False:
            msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                       f"but could not find file: {path_noaa}")
            logging.error(msg_cfv)
            print(msg_cfv)
            cf.remove_handlers_if_directly_executed(func_1up)
            raise Exception(msg_cfv)
    
    # Open intermediate output data file if it exists, otherwise create the file.
    
    if Path(path_noaa).exists():
        msg_open = f"Opening: existing file for use in {func_cur}: {path_noaa}"
        logging.info(msg_open)
        print(msg_open)
    else:
        cf.proc_noaa_ind()
    ds_noaa = xr.open_dataset(path_noaa, engine = "netcdf4")
    
    if cf.priority == "speed":
        ds_noaa = ds_noaa.persist()
    
    # Obtain time_start and time_end of data used in computation of rolling
    # averages between year_start and year_end, then obtain the correct string 
    # representations for variables used in file output name. 
    
    time_start = (datetime.strptime(str(year_start), "%Y") + 
                  relativedelta(years=-(window_size-1)/2))
    time_end = (datetime.strptime(str(year_end), "%Y") + 
                relativedelta(years=(window_size-1)/2+1, months=-1))
    period1_mid_str = str(period1_mid)
    period2_mid_str = str(period2_mid)
    month1_mark_str = str(month1_mark)
    month2_mark_str = str(month2_mark)
    
    if period1_mid:
        period1_mid = datetime.strptime(period1_mid, "%b-%Y")
        period1_start = period1_mid + relativedelta(years=-(window_size-1)/2, 
                                                    months=-6)
        period1_end = period1_mid + relativedelta(years=(window_size-1)/2, 
                                                  months=5)
        period1_start_str = period1_start.strftime("%b-%Y")
        period1_end_str = period1_end.strftime("%b-%Y")
        
    if period2_mid:
        period2_mid = datetime.strptime(period2_mid, "%b-%Y")
        period2_start = period2_mid + relativedelta(years=-(window_size-1)/2, 
                                                    months=-6)
        period2_end = period2_mid + relativedelta(years=(window_size-1)/2, 
                                                  months=5)
        period2_start_str = period2_start.strftime("%b-%Y")
        period2_end_str = period2_end.strftime("%b-%Y")
        
    if month1_mark:
        month1_mark = datetime.strptime(month1_mark, "%b-%Y")
        
    if month2_mark:
        month2_mark = datetime.strptime(month2_mark, "%b-%Y")

    def filter_dates_event(dates):
        # Filter function to be used in process_dates_event function to select out
        # the list of ENSO and IOD events which are contained within the period
        # between year_start and year_end.
        if (datetime.strptime(dates[1], "%b-%Y") + relativedelta(months=1, hours=-1) < 
            time_start):
            return False
        elif datetime.strptime(dates[0], "%b-%Y") > time_end:
            return False
        else:
            return True
        
    def process_dates_event(dates):
        # Function to select out the list of ENSO and IOD events which are contained 
        # within the period between year_start and year_end, then truncate any parts of 
        # events which extend to before or after year_start and year_end respectively.
        dates_processed = list(filter(filter_dates_event, copy.deepcopy(dates)))
        if len(dates_processed) == 0:
            return dates_processed
        for idx in range(0, len(dates_processed)):
            dates_processed[idx][0] = datetime.strptime(dates_processed[idx][0], "%b-%Y")
            dates_processed[idx][1] = (datetime.strptime(dates_processed[idx][1], 
                                                         "%b-%Y"))
        if dates_processed[0][0] < time_start:
            dates_processed[0][0] = time_start
        if dates_processed[-1][1] > time_end:
            dates_processed[-1][1] = time_end
        return dates_processed

    # Select out the list of ENSO and IOD events which are contained within the
    # period between year_start and year_end, then truncate any parts of events 
    # which extend to before or after year_start and year_end respectively.
    
    dates_la_nina_processed = process_dates_event(dates_la_nina)
    dates_el_nino_processed = process_dates_event(dates_el_nino)
    dates_neg_iod_processed = process_dates_event(dates_neg_iod)
    dates_pos_iod_processed = process_dates_event(dates_pos_iod)
    
    # Compute rolling averages for indices.
    
    ds_noaa_roll = (ds_noaa
                    # The use of min_periods and skipna below is to get around problem
                    # with EPOI data having missing values for December.
                    .rolling(time = 12 * window_size, center = True, 
                             min_periods = 11 * window_size)
                    .mean(skipna = True)
                    .sel(time = slice(time_start, time_end))
                   )
    
    # Obtain dataframe representation for rolling averages of indices (this
    # makes it easier to extract date labels for plotting).
    
    df_noaa = (ds_noaa
               .sel(time = slice(time_start, time_end))
               .to_dataframe()
              )

    # Set minor xticks to have annual frequency, while major xticks is automated
    # based on how many years are being covered in the plot.
    
    xticks_minor = df_noaa.index[::12]
    xticks_major = xticks_minor[::math.ceil((year_end - year_start) / 20)]
    
    # Create plot
    
    figcols = 1
    figrows = len(ds_noaa.keys())
    figwidth = figwidth_standard * 2
    figheight = figwidth / 3 * figrows
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight))

    for row, index in enumerate(ds_noaa.keys()):
        ax = axes[row]
        index_attrs = ds_noaa[index].attrs

        # Create bar plot overlayed with line plot for monthly values and the rolling
        # average over these monthly values respectively.
        
        ax.bar(df_noaa.index, df_noaa[index], width=bar_width, color="gray", alpha=0.5, 
               label="Monthly Values")
        ds_noaa_roll[index].plot(ax=ax, color="k", 
                                 label=f"{window_size}-Year Rolling Average (centred)")

        ax.set_xticks(xticks_minor, minor = True)
        ax.set_xticks(xticks_major)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b-%Y'))
        ax.set_xlabel(None)
        ax.set_ylabel("{} [{}]"
                      .format(index_attrs["abbreviation"], index_attrs["units"]))
        # ax.set_title(chr(ord('`')+(row+1)) + ") {} ({} data)"
        #              .format(index_attrs["full_name"], index_attrs["source"]))
        ax.set_title("{} ({} data)"
                     .format(index_attrs["full_name"], index_attrs["source"]))
        
        # Add highlights for ENSO events.
        
        if index == "oni":
            for idx, dates_list in enumerate(dates_la_nina_processed):
                ax.axvspan(dates_list[0], dates_list[1], color="blue", alpha=0.05, 
                           label="_"*idx+"La Nina (JMA data)")
            for idx, dates_list in enumerate(dates_el_nino_processed):
                ax.axvspan(dates_list[0], dates_list[1], color="red", alpha=0.05, 
                           label="_"*idx+"El Nino (JMA data)")
        
        # Add highlights for IOD events.
        
        if index == "dmi":
            for idx, dates_list in enumerate(dates_neg_iod_processed):
                ax.axvspan(dates_list[0], dates_list[1], color="blue", alpha=0.05, 
                           label="_"*idx+"Negative IOD (JMA data)")
            for idx, dates_list in enumerate(dates_pos_iod_processed):
                ax.axvspan(dates_list[0], dates_list[1], color="red", alpha=0.05, 
                           label="_"*idx+"Positive IOD (JMA data)")
                
        # Highlight the rolling window and rolling average centred on period1_mid.
                
        if period1_mid:
            period1_avg = float(ds_noaa_roll[index].sel(time=period1_mid).data)
            ax.axvspan(period1_start, period1_end, color="green", alpha=0.15, 
                       label=f"Period 1: {period1_start_str} to " +
                       f"{period1_end_str} (inclusive)")
            ax.plot(period1_mid, period1_avg, marker="X", markersize=10, color="green",
                    label="{}-Year Average over Period 1 = {}"
                    .format(window_size, round(period1_avg, 3)))
        
        # Highlight the rolling window and rolling average centred on period2_mid.
        
        if period2_mid:
            period2_avg = float(ds_noaa_roll[index].sel(time=period2_mid).data)
            ax.axvspan(period2_start, period2_end, color="green", alpha=0.15, 
                       label=f"Period 2: {period2_start_str} to " +
                       f"{period2_end_str} (inclusive)")
            ax.plot(period2_mid, period2_avg, marker="X", markersize=10, color="green",
                    label="{}-Year Average over Period 2 = {}"
                    .format(window_size, round(period2_avg, 3)))
            
        # Add marker for monthly value at month1_mark.
            
        if month1_mark:
            month1_value = float(ds_noaa[index].sel(time=month1_mark).data)
            ax.plot(month1_mark, month1_value, marker="P", markersize=10, color="purple", 
                    label="Value for {} = {}"
                    .format(month1_mark_str, round(month1_value, 3)))
        
        # Add marker for monthly value at month2_mark.
        
        if month2_mark:
            month2_value = float(ds_noaa[index].sel(time=month2_mark).data)
            ax.plot(month2_mark, month2_value, marker="P", markersize=10, color="purple", 
                    label="Value for {} = {}"
                    .format(month2_mark_str, round(month2_value, 3)))

        ax.legend(loc="upper right")
            
    fig.tight_layout()
        
    if output == True:
        if cfv_data:
            cfv_used = copy.deepcopy(cfv_data)
        else:
            cfv_used = copy.deepcopy(cf.calc_funcs_ver)
        
        path_output = (path_noaa
                       .replace("data_processed", "data_final")
                       .replace(".nc", f"_{year_start}_{year_end}_{window_size}_" +
                                f"mid1-{period1_mid_str}_mid2-{period2_mid_str}_" +
                                f"mark1-{month1_mark_str}_mark2-{month2_mark_str}.png")
                       .replace(f"{cfv_used}_proc_", 
                                f"{plot_funcs_ver}_{cfv_used}_proc_")
                      )
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
            
        if Path(path_output).exists():
            msg_exist = ("WARNING: plot file already exists (and was " +
                         f"not overwritten): {path_output}")
            logging.warning(msg_exist)
            print(msg_exist)
        else:
            plt.savefig(path_output, metadata=get_plot_metadata(
                time_exec, func_cur, args_cur, args_cur_values)
                       )
            msg_create = f"CREATED: plot file: {path_output}"
            logging.info(msg_create)
            print(msg_create)
                
    plt.show()
    fig.clear()
    plt.close(fig)
        
    cf.remove_handlers_if_directly_executed(func_1up)


# ## High-level plotting functions

# ### Calc plots

# In[ ]:


def plot_mdp_clim_stats_given_var_or_dvar(
    region, period_start, period_end, period_months, glass_source_pref, var_or_dvar, 
    extents=None, cfv_data=None, output=False
):
    
    """
    Create plot for the mean diurnal profile (MDP) stats over a period.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_months (str or list): Months subset of period to perform calculation over.
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
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/mdp_clim_stats_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_calc_
        {extents_used}_{period_start}_{period_end}_{period_months_str}_stats_
        {var_or_dvar}.png:
            Output PNG file in data_final folder for the plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period_months_str} and
            is a string representing the list of selected months to use as a subset.
    
    For each period, plot the max, min, mean, range, hour of max (24-hour Local Time) 
    and hour of min (24-hour Local Time), for the mean diurnal profile (MDP) of a given 
    variable for each grid cell. The MDPs are computed over the period between 
    period1_start and period1_end, and period2_start and period2_end (inclusive), 
    and only using a subset of data within these periods (if period1_months and 
    period2_months respectively not "all" is specified). Also included for reference 
    are subplots for mean leaf area index (MLAI) and mean fraction of absorbed 
    photosynthetically active radiation (MFAPAR).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period_start=period_start, period_end=period_end, 
                  period_months=period_months, glass_source_pref=glass_source_pref,
                  var_or_dvar=var_or_dvar, extents=extents, cfv_data=cfv_data, 
                  output=output)
    
    # Obtain string representations for month subsets.
    
    period_months_str = cf.get_period_months_str(period_months=period_months)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)
            
    path_output = (f"../data_final/mdp_clim_stats_given_var_or_dvar/" +
                   f"{plot_funcs_ver}_{cfv_used}_calc_{extents_used}_" +
                   f"{period_start}_{period_end}_{period_months_str}_" +
                   f"stats_{var_or_dvar}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months=period_months, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, extents=extents, 
        ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months=period_months, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, extents=extents, 
        ax=axes[1][1], cfv_data=cfv_data
    )
    
    # Remove all stats until we have the remaining 6 we want to create subplots for.
    
    stats_to_plot = copy.deepcopy(cf.params_stat)
    
    for stat in ["max_u", "max_v", "min_u", "min_v", "mean_u", "mean_v"]:
        try:
            stats_to_plot.remove(stat)
        except:
            pass
    
    # Plot these stats in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, stat in enumerate(stats_to_plot):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_calc_plot(
            calc_func=cf.calc_era5_mdp_clim_stats_given_var_or_dvar, region=region, 
            period_start=period_start, period_end=period_end, 
            period_months=period_months, arg_extra=stat, glass_source_pref=None, 
            var_or_dvar=var_or_dvar, extents=extents, ax=axes[row][col], 
            cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_means_given_layer_and_type(
    region, period_start, period_end, period_months, period_hours, 
    glass_source_pref, var_or_dvar_layer, var_or_dvar_type, 
    extents=None, cfv_data=None, output=False
):
    
    """
    Create plot for the means over a period for an assortment of different parameters 
    specified by parameter layer and type.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_months (str or list): Months subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period_hours (str or list): Hours subset of period to perform calculation over.
            Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        var_or_dvar_layer (str): Spatial layer from which to draw ERA5 parameters for 
            analysis. This is used for the plot_funcs script. Must be one of: 
            ["sfc", "atm", "cld"].
        var_or_dvar_type (str): Whether to analyse the variables themselves or the 
            change in their mean diurnal profile values as compared with their values
            in the previous hour. This is used for the plot_funcs script.
            Must be one of: ["vars", "dvars"].
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/means_given_layer_and_type/{plot_funcs_ver}_{cfv_used}_calc_
        {extents_used}_{period_start}_{period_end}_means_{period_months_str}_
        {period_hours_str}_{var_or_dvar_layer}_{var_or_dvar_type}.png:
            Output PNG file in data_final folder for the plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period_months_str} is a
            string representing the list of selected months to use as a subset. 
            {period_hours_str} is a string representing the list of selected 
            hours to use as a subset.
    
    For each period, plot the mean values for an assortment of parameters which are 
    categorised within the surface ("sfc"), atmosphere ("atm") or cloud ("cld") layers, 
    and are either variable values themselves or an hourly change in the variable values. 
    The means are computed over the period between period1_start and period1_end, and 
    period2_start and period2_end (inclusive), using a subset of months (if period1_months 
    and period2_months respectively not "all" is  specified) and hours 
    (if period1_hours and period2_hours respectively not "all" is specified) within
    each period. Also included for reference are subplots for mean leaf area index (MLAI) 
    and mean fraction of absorbed photosynthetically active radiation (MFAPAR).
    
    If var_or_dvar_layer == “sfc” and var_or_dvar_type == “vars”, then the wind speed at 
    100 m above surface (WS100), wind velocity at 100 m above surface (WV100), mean sea 
    level pressure (MSLP), temperature at 2 m above surface (T2), surface latent heat flux 
    (SLHF) and surface sensible heat flux (SSHF) will be plotted. If var_or_dvar_layer == 
    “sfc” and var_or_dvar_type == “dvars”, then the change in each of these variables as 
    compared with the value in the previous hour (dWS100, dWV100, dMSLP, dT2, dSLHF, 
    dSSHF) will be plotted instead.
    
    If var_or_dvar_layer == “atm” and var_or_dvar_type == “vars”, then the vertical integral 
    of energy conversion (VIEC), vertical integral of potential, internal and latent energy 
    (VIPILE), vertical integral of kinetic energy (VIKE), total column cloud liquid water 
    (TCCLW), total column water vapour (TCWV) and net atmospheric condensation (NAC) will be 
    plotted. If var_or_dvar_layer == “atm” and var_or_dvar_type == “dvars”, then the change 
    in each of these variables as compared with the value in the previous hour (dVIEC, 
    dVIPILE, dVIKE, dTCCLW, dTCWV, dNAC) will be plotted instead.
    
    If var_or_dvar_layer == “cld” and var_or_dvar_type == “vars”, then the boundary layer
    height (BLH), forecast albedo (FA), cloud base height (CBH), total cloud cover (TCC),
    convective available potential energy (CAPE) and convective inhibition (CI) will be
    plotted. If var_or_dvar_layer == “cld” and var_or_dvar_type == “dvars”, then the change 
    in each of these variables as compared with the value in the previous hour (dBLH, dFA, 
    dCBH, dTCC, dCAPE, dCI) will be plotted instead.
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period_start=period_start, period_end=period_end, 
                  period_months=period_months, period_hours=period_hours,
                  glass_source_pref=glass_source_pref,
                  var_or_dvar_layer=var_or_dvar_layer, 
                  var_or_dvar_type=var_or_dvar_type, extents=extents,
                  cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month and hour subsets.
    
    period_months_str = cf.get_period_months_str(period_months=period_months)
    period_hours_str = cf.get_period_hours_str(period_hours=period_hours)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    path_output = (f"../data_final/means_given_layer_and_type/" +
                   f"{plot_funcs_ver}_{cfv_used}_calc_{extents_used}_" +
                   f"{period_start}_{period_end}_means_{period_months_str}_" +
                   f"{period_hours_str}_{var_or_dvar_layer}_" +
                   f"{var_or_dvar_type}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months=period_months, period_hours=None,
        arg_extra="mlai", glass_source_pref=glass_source_pref, var_or_dvar=None, 
        extents=extents, ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months=period_months, period_hours=None,
        arg_extra="mfapar", glass_source_pref=glass_source_pref, var_or_dvar=None, 
        extents=extents, ax=axes[1][1], cfv_data=cfv_data
    )
    
    # Remove all params until we have the remaining 6 we want to create subplots for.
    
    params_to_plot = copy.deepcopy(
        cf.vars_and_dvars_era5[var_or_dvar_type][var_or_dvar_layer])
    
    for param in ["u10", "v10", "ws10", "wv10", "u100", "v100", "nse", "vidmf", "du10", 
                  "dv10", "dws10", "dwv10", "du100", "dv100", "dnse", "dvidmf"]:
        try:
            params_to_plot.remove(param)
        except:
            pass
    
    # Plot these params in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, param in enumerate(params_to_plot):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_calc_plot(
            calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
            period_start=period_start, period_end=period_end, 
            period_months=period_months, period_hours=period_hours,
            arg_extra="mean", glass_source_pref=None, var_or_dvar=param, 
            extents=extents, ax=axes[row][col], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_hourly_means_given_var_or_dvar(
    region, period_start, period_end, period_months, glass_source_pref, 
    var_or_dvar, hours_to_plot, extents=None, cfv_data=None, output=False
):
    
    """
    Create plot for the hourly mean values of var_or_dvar over a period.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_months (str or list): Months subset of period to perform calculation over.
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
        hours_to_plot (str): Which hours (in local time) to display results for in 
            plot. Must be one of: ["0-5", "6-11", "12-17", "18-23"].
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/hourly_means_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_calc_
        {extents_used}_{period_start}_{period_end}_means-hourly_{period_months_str}_
        {var_or_dvar}_{hours_to_plot}.png:
            Output PNG file in data_final folder for the plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period_months_str} is a
            string representing the list of selected months to use as a subset.
    
    For each period, plot the hourly mean values for a given variable or its change since 
    the previous hour for Local Time: 0000 to 0500, 0600 to 1100, 1200 to 1700, 
    or 1800-2300 (inclusive) at each grid cell. The hourly mean values are computed over 
    the period between period1_start and period1_end, and period2_start and period2_end 
    (inclusive), and only using a subset of months within these periods (if 
    period1_months and period2_months respectively not "all" is specified). 
    Also included for reference are subplots for mean leaf area index (MLAI) 
    and mean fraction of absorbed photosynthetically active radiation (MFAPAR).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period_start=period_start, period_end=period_end, 
                  period_months=period_months, glass_source_pref=glass_source_pref,
                  var_or_dvar=var_or_dvar, hours_to_plot=hours_to_plot, 
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month subsets.
    
    period_months_str = cf.get_period_months_str(period_months=period_months)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    path_output = (f"../data_final/hourly_means_given_var_or_dvar/" +
                   f"{plot_funcs_ver}_{cfv_used}_calc_{extents_used}_" +
                   f"{period_start}_{period_end}_means-hourly_" +
                   f"{period_months_str}_{var_or_dvar}_{hours_to_plot}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months=period_months, period_hours=None,
        arg_extra="mlai", glass_source_pref=glass_source_pref, var_or_dvar=None, 
        extents=extents, ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months=period_months, period_hours=None,
        arg_extra="mfapar", glass_source_pref=glass_source_pref, var_or_dvar=None, 
        extents=extents, ax=axes[1][1], cfv_data=cfv_data
    )
    
    # For each hour in hours_to_plot, open the relevant datasets then use them to
    # determine colourbar and colourbar extents. If the intermediate output data
    # files don't yet exist, then create them.
    
    hours = cf.hour_subsets[hours_to_plot]
    datasets = []
    for hour in cf.hour_subsets["all"]:
    
        path_calc = cf.get_path_for_calc_func(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period_start=period_start, period_end=period_end, period_months=period_months,
            period_hours=[hour], glass_source_pref=None, var_or_dvar=var_or_dvar
        )
    
        # Use data outputted from an older version of the calc_funcs script. This is useful
        # for results which required computationally intensive processing. And can also be
        # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
        if cfv_data:
            path_calc = path_calc.replace(cf.calc_funcs_ver, cfv_data)
            if Path(path_calc).exists() == False:
                msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                           f"but could not find file: {path_calc}")
                logging.error(msg_cfv)
                print(msg_cfv)
                cf.remove_handlers_if_directly_executed(func_1up)
                raise Exception(msg_cfv)
    
        if Path(path_calc).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_calc}"
            logging.info(msg_open)
        else:
            cf.calc_era5_mean_clim_given_var_or_dvar(
                region=region, period_start=period_start, period_end=period_end, 
                period_months=period_months, period_hours=[hour],
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
    
        ds_hour = (xr.open_dataset(path_calc, engine = "netcdf4")
                   .sel(longitude=slice(extents[0], extents[1]), 
                        latitude=slice(extents[3], extents[2]))
                   .expand_dims({"hour": [hour]})
                  )
        datasets.append(ds_hour)
    
    ds_hour = xr.merge(datasets)
    
    if var_or_dvar in cf.params_vector:
        da_hour_u = ds_hour["mean_u"]
        da_hour_v = ds_hour["mean_v"]
        da_hour_mag = cf.get_magnitude(da_hour_u, da_hour_v)
        vmin = float(da_hour_mag.min())
        vmax = float(da_hour_mag.max())
    else:
        da_hour = ds_hour["mean"]
        if var_or_dvar in vars_pos_with_vmin_0:
            vmin = 0
            vmax = float(da_hour.max())
        elif var_or_dvar in vars_pos:
            vmin = float(da_hour.min())
            vmax = float(da_hour.max())
        else:
            min_da = float(da_hour.min())
            max_da = float(da_hour.max())
            vmin = min(-abs(min_da), -abs(max_da))
            vmax = -vmin
    
    if vmin != None:
        if math.isnan(vmin):
            vmin = None
    if vmax != None:
        if math.isnan(vmax):
            vmax = None
    
    # Plot results for these hours in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, hour in enumerate(hours):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_calc_plot(
            calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
            period_start=period_start, period_end=period_end, 
            period_months=period_months, period_hours=[hour], arg_extra="mean", 
            glass_source_pref=None, var_or_dvar=var_or_dvar, extents=extents, 
            vmin=vmin, vmax=vmax, ax=axes[row][col], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_monthly_means_given_var_or_dvar(
    region, period_start, period_end, period_hours, glass_source_pref, 
    var_or_dvar, months_to_plot, extents=None, cfv_data=None, output=False
):
    
    """
    Create plot for the monthly mean values of var_or_dvar over a period.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_hours (str or list): Hours subset of period to perform calculation over.
            Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
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
        months_to_plot (str): Which months to display results for in plot. Must be
            one of: ["1-6", "7-12"].
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/monthly_means_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_calc_
        {extents_used}_{period_start}_{period_end}_means-monthly_{period_hours_str}_
        {var_or_dvar}_{months_to_plot}.png:
            Output PNG file in data_final folder for the plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period_hours_str} is a
            string representing the list of selected hours to use as a subset.
    
    For each period, plot the monthly mean values for a given variable or its change since 
    the previous hour for months: Jan-Jun (1-6) or Jul-Dec (7-12) inclusive at each grid 
    cell. The monthly mean values are computed over the period between period1_start and 
    period1_end, and period2_start and period2_end (inclusive), and only using a subset 
    of hours within these periods (if period1_hours and period2_hours respectively not 
    "all" is specified). Also included for reference are subplots for mean leaf area 
    index (MLAI) and mean fraction of absorbed photosynthetically active radiation (MFAPAR).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period_start=period_start, period_end=period_end, 
                  period_hours=period_hours, glass_source_pref=glass_source_pref,
                  var_or_dvar=var_or_dvar, months_to_plot=months_to_plot, 
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for hour subsets.
    
    period_hours_str = cf.get_period_hours_str(period_hours=period_hours)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    path_output = (f"../data_final/monthly_means_given_var_or_dvar/" +
                   f"{plot_funcs_ver}_{cfv_used}_calc_{extents_used}_" +
                   f"{period_start}_{period_end}_means-monthly_" +
                   f"{period_hours_str}_{var_or_dvar}_{months_to_plot}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months="all", period_hours=None,
        arg_extra="mlai", glass_source_pref=glass_source_pref, var_or_dvar=None, 
        extents=extents, ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months="all", period_hours=None,
        arg_extra="mfapar", glass_source_pref=glass_source_pref, var_or_dvar=None, 
        extents=extents, ax=axes[1][1], cfv_data=cfv_data
    )
    
    # For each month in months_to_plot, open the relevant datasets then use them to
    # determine colourbar and colourbar extents. If the intermediate output data
    # files don't yet exist, then create them.
    
    months = cf.month_subsets[months_to_plot]
    datasets = []
    for month in cf.month_subsets["all"]:
    
        path_calc = cf.get_path_for_calc_func(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period_start=period_start, period_end=period_end, period_months=[month],
            period_hours=period_hours, glass_source_pref=None, var_or_dvar=var_or_dvar
        )
    
        # Use data outputted from an older version of the calc_funcs script. This is useful
        # for results which required computationally intensive processing. And can also be
        # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
        if cfv_data:
            path_calc = path_calc.replace(cf.calc_funcs_ver, cfv_data)
            if Path(path_calc).exists() == False:
                msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                           f"but could not find file: {path_calc}")
                logging.error(msg_cfv)
                print(msg_cfv)
                cf.remove_handlers_if_directly_executed(func_1up)
                raise Exception(msg_cfv)
    
        if Path(path_calc).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_calc}"
            logging.info(msg_open)
        else:
            cf.calc_era5_mean_clim_given_var_or_dvar(
                region=region, period_start=period_start, period_end=period_end, 
                period_months=[month], period_hours=period_hours,
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
    
        ds_month = (xr.open_dataset(path_calc, engine = "netcdf4")
                    .sel(longitude=slice(extents[0], extents[1]), 
                         latitude=slice(extents[3], extents[2]))
                    .expand_dims({"month": [month]})
                   )
        datasets.append(ds_month)
    
    ds_month = xr.merge(datasets)
    
    if var_or_dvar in cf.params_vector:
        da_month_u = ds_month["mean_u"]
        da_month_v = ds_month["mean_v"]
        da_month_mag = cf.get_magnitude(da_month_u, da_month_v)
        vmin = float(da_month_mag.min())
        vmax = float(da_month_mag.max())
    else:
        da_month = ds_month["mean"]
        if var_or_dvar in vars_pos_with_vmin_0:
            vmin = 0
            vmax = float(da_month.max())
        elif var_or_dvar in vars_pos:
            vmin = float(da_month.min())
            vmax = float(da_month.max())
        else:
            min_da = float(da_month.min())
            max_da = float(da_month.max())
            vmin = min(-abs(min_da), -abs(max_da))
            vmax = -vmin
    
    if vmin != None:
        if math.isnan(vmin):
            vmin = None
    if vmax != None:
        if math.isnan(vmax):
            vmax = None
    
    # Plot results for these months in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, month in enumerate(months):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_calc_plot(
            calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
            period_start=period_start, period_end=period_end, 
            period_months=[month], period_hours=period_hours, arg_extra="mean", 
            glass_source_pref=None, var_or_dvar=var_or_dvar, extents=extents, 
            vmin=vmin, vmax=vmax, ax=axes[row][col], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_wsd_clim(
    region, period_start, period_end, period_months, period_hours, glass_source_pref,  
    extents=None, cfv_data=None, output=False
):
    
    """
    Create plot for the 100 m wind speed distribution parameters over a period.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of: ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_months (str or list): Months subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period_hours (str or list): Hours subset of period to perform calculation over.
            Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/wsd_clim/{plot_funcs_ver}_{cfv_used}_calc_{extents_used}_
        {period_start}_{period_end}_{period_months_str}_wsd_{period_hours_str}.png:
            Output PNG file in data_final folder for the plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period_months_str} is a
            string representing the list of selected months to use as a subset. 
            {period_hours_str} is a string representing the list of selected hours 
            to use as a subset.
    
    For each period, plot the mean and standard deviation of wind speed at 100 m above 
    surface, the Weibull scale and shape parameter for wind speed at 100 m above surface 
    (C100 and K100), the expected rate of exceedance for a particular wind speed at 100 m 
    above surface (EROE100) and the gross capacity factor for a typical wind turbine at 
    100 m above surface (TGCF100) for each grid cell. The wind speed distributions (WSDs) 
    are computed over the period between period1_start and period1_end, and period2_start 
    and period2_end (inclusive), using a subset of months (if period1_months and 
    period2_months respectively not "all" is  specified) and hours (if period1_hours and 
    period2_hours respectively not "all" is specified) within each period. Also included 
    for reference are subplots for mean leaf area index (MLAI) and mean fraction of 
    absorbed photosynthetically active radiation (MFAPAR).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period_start=period_start, period_end=period_end, 
                  period_months=period_months, period_hours=period_hours,
                  glass_source_pref=glass_source_pref, 
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month and hour subsets.
    
    period_months_str = cf.get_period_months_str(period_months=period_months)
    period_hours_str = cf.get_period_hours_str(period_hours=period_hours)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    path_output = (f"../data_final/wsd_clim/{plot_funcs_ver}_{cfv_used}_" +
                   f"calc_{extents_used}_{period_start}_{period_end}_" +
                   f"{period_months_str}_wsd_{period_hours_str}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months=period_months, period_hours=period_hours,
        arg_extra="mlai", glass_source_pref=glass_source_pref, var_or_dvar=None, 
        extents=extents, ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_calc_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period_start=period_start, 
        period_end=period_end, period_months=period_months, period_hours=period_hours,
        arg_extra="mfapar", glass_source_pref=glass_source_pref, var_or_dvar=None, 
        extents=extents, ax=axes[1][1], cfv_data=cfv_data
    )
    
    # Remove all params until we have the remaining 6 we want to create subplots for.
    
    params_to_plot = copy.deepcopy(cf.params_wsd)
    
    for param in ["ws10_mean", "ws10_std", "c10", "k10"]:
        try:
            params_to_plot.remove(param)
        except:
            pass
    
    # Plot these params in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, param in enumerate(params_to_plot):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_calc_plot(
            calc_func=cf.calc_era5_wsd_clim, region=region, period_start=period_start, 
            period_end=period_end, period_months=period_months, period_hours=period_hours,
            arg_extra=param, glass_source_pref=None, var_or_dvar=None, 
            extents=extents, ax=axes[row][col], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# ### Diff plots

# In[ ]:


def plot_diff_mdp_clim_stats_given_var_or_dvar(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, glass_source_pref, var_or_dvar, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    extents=None, cfv_data=None, output=False
):
    
    """
    Create difference plot for the mean diurnal profile (MDP) stats between
    two different periods.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
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
        perc (bool): Whether to plot the difference in values as a percentage of the
            value (magnitude if negative) in period1. This is used for the comp plots
            in the plot_funcs script. Must be one of: [True, False].
        mask_perc_quantile (int): If perc is True, specify the quantile of values
            (magnitude if negative) from period1 to mask for the difference plot.
            This is used because percentage differences may be particularly high
            for values which had a low magnitude as a base in period 1.
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/mdp_clim_stats_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_diff_
        {extents_used}_{period1_start}_{period1_end}_{period2_start}_{period2_end}_
        {period1_months_str}_{period2_months_str}_stats_{var_or_dvar}_
        perc-{quantile_used}.png:
            Output PNG file in data_final folder for the difference plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_months_str} and
            {period2_months_str} are strings representing the list of selected months 
            to use as a subset in each period. {quantile_used} is equal to 
            mask_perc_quantile if perc = True, otherwise it is set to None.
    
    For each period, calculate the max, min, mean, range, hour of max (24-hour Local Time) 
    and hour of min (24-hour Local Time), for the mean diurnal profile (MDP) of a given 
    variable for each grid cell, then plot the difference in results between the two 
    periods. The MDPs are computed over the period between period1_start and period1_end, 
    and period2_start and period2_end (inclusive), and only using a subset of data within 
    these periods (if period1_months and period2_months respectively not "all" is 
    specified). Also included for reference are subplots for mean leaf area index (MLAI) 
    and mean fraction of absorbed photosynthetically active radiation (MFAPAR).
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    The mask applies only to scalar variables (i.e. quiver plots for vectors are not 
    masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end, 
                  period2_start=period2_start, period2_end=period2_end,
                  period1_months=period1_months, period2_months=period2_months, 
                  glass_source_pref=glass_source_pref, var_or_dvar=var_or_dvar, 
                  perc=perc, mask_perc_quantile=mask_perc_quantile, extents=extents, 
                  cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month subsets.
    
    period1_months_str = cf.get_period_months_str(period_months=period1_months)
    period2_months_str = cf.get_period_months_str(period_months=period2_months)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/mdp_clim_stats_given_var_or_dvar/" +
                   f"{plot_funcs_ver}_{cfv_used}_diff_{extents_used}_" +
                   f"{period1_start}_{period1_end}_{period2_start}_" +
                   f"{period2_end}_{period1_months_str}_{period2_months_str}_" +
                   f"stats_{var_or_dvar}_perc-{quantile_used}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, extents=extents, 
        ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, extents=extents, 
        ax=axes[1][1], cfv_data=cfv_data
    )
    
    # Remove all stats until we have the remaining 6 we want to create subplots for.
    
    stats_to_plot = copy.deepcopy(cf.params_stat)
    
    for stat in ["max_u", "max_v", "min_u", "min_v", "mean_u", "mean_v"]:
        try:
            stats_to_plot.remove(stat)
        except:
            pass
    
    # Plot these stats in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, stat in enumerate(stats_to_plot):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_diff_plot(
            calc_func=cf.calc_era5_mdp_clim_stats_given_var_or_dvar, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months, 
            arg_extra=stat, glass_source_pref=None, var_or_dvar=var_or_dvar, 
            perc=perc, mask_perc_quantile=mask_perc_quantile,
            extents=extents, ax=axes[row][col], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_diff_means_given_layer_and_type(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, period1_hours, period2_hours, 
    glass_source_pref, var_or_dvar_layer, var_or_dvar_type, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    extents=None, cfv_data=None, output=False
):
    
    """
    Create difference plot for the means between two different periods for an
    assortment of different parameters specified by parameter layer and type.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
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
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/means_given_layer_and_type/{plot_funcs_ver}_{cfv_used}_diff_
        {extents_used}_{period1_start}_{period1_end}_{period2_start}_{period2_end}_
        means_{period1_months_str}_{period2_months_str}_{period1_hours_str}_
        {period2_hours_str}_{var_or_dvar_layer}_{var_or_dvar_type}_
        perc-{quantile_used}.png:
            Output PNG file in data_final folder for the difference plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_months_str} and
            {period2_months_str} are strings representing the list of selected months 
            to use as a subset in each period. {period1_hours_str} and {period2_hours_str} 
            are strings representing the list of selected hours to use as a subset in 
            each period. {quantile_used} is equal to mask_perc_quantile if perc = True, 
            otherwise it is set to None.
    
    Plot the difference in the mean values over period1 to period2 for an assortment 
    of parameters which are categorised within the surface ("sfc"), atmosphere ("atm") 
    or cloud ("cld") layers, and are either variable values themselves or an hourly change 
    in the variable values. The means are computed over the period between period1_start 
    and period1_end, and period2_start and period2_end (inclusive), using a subset of months 
    (if period1_months and period2_months respectively not "all" is specified) and hours 
    (if period1_hours and period2_hours respectively not "all" is specified) within
    each period. Also included for reference are subplots for mean leaf area index (MLAI) 
    and mean fraction of absorbed photosynthetically active radiation (MFAPAR).
    
    If var_or_dvar_layer == “sfc” and var_or_dvar_type == “vars”, then the wind speed at 
    100 m above surface (WS100), wind velocity at 100 m above surface (WV100), mean sea 
    level pressure (MSLP), temperature at 2 m above surface (T2), surface latent heat flux 
    (SLHF) and surface sensible heat flux (SSHF) will be plotted. If var_or_dvar_layer == 
    “sfc” and var_or_dvar_type == “dvars”, then the change in each of these variables as 
    compared with the value in the previous hour (dWS100, dWV100, dMSLP, dT2, dSLHF, 
    dSSHF) will be plotted instead.
    
    If var_or_dvar_layer == “atm” and var_or_dvar_type == “vars”, then the vertical integral 
    of energy conversion (VIEC), vertical integral of potential, internal and latent energy 
    (VIPILE), vertical integral of kinetic energy (VIKE), total column cloud liquid water 
    (TCCLW), total column water vapour (TCWV) and net atmospheric condensation (NAC) will be 
    plotted. If var_or_dvar_layer == “atm” and var_or_dvar_type == “dvars”, then the change 
    in each of these variables as compared with the value in the previous hour (dVIEC, 
    dVIPILE, dVIKE, dTCCLW, dTCWV, dNAC) will be plotted instead.
    
    If var_or_dvar_layer == “cld” and var_or_dvar_type == “vars”, then the boundary layer
    height (BLH), forecast albedo (FA), cloud base height (CBH), total cloud cover (TCC),
    convective available potential energy (CAPE) and convective inhibition (CI) will be
    plotted. If var_or_dvar_layer == “cld” and var_or_dvar_type == “dvars”, then the change 
    in each of these variables as compared with the value in the previous hour (dBLH, dFA, 
    dCBH, dTCC, dCAPE, dCI) will be plotted instead.
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    The mask applies only to scalar variables (i.e. quiver plots for vectors are not 
    masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end, 
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_months=period1_months, period2_months=period2_months, 
                  period1_hours=period1_hours, period2_hours=period2_hours,
                  glass_source_pref=glass_source_pref,
                  var_or_dvar_layer=var_or_dvar_layer, var_or_dvar_type=var_or_dvar_type, 
                  perc=perc, mask_perc_quantile=mask_perc_quantile,
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month and hour subsets.
    
    period1_months_str = cf.get_period_months_str(period_months=period1_months)
    period2_months_str = cf.get_period_months_str(period_months=period2_months)
    period1_hours_str = cf.get_period_hours_str(period_hours=period1_hours)
    period2_hours_str = cf.get_period_hours_str(period_hours=period2_hours)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/means_given_layer_and_type/" +
                   f"{plot_funcs_ver}_{cfv_used}_diff_{extents_used}_" +
                   f"{period1_start}_{period1_end}_{period2_start}_{period2_end}_" +
                   f"means_{period1_months_str}_{period2_months_str}_" +
                   f"{period1_hours_str}_{period2_hours_str}_" +
                   f"{var_or_dvar_layer}_{var_or_dvar_type}_perc-{quantile_used}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, 
        period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=None, period2_hours=None, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, 
        perc=perc, mask_perc_quantile=mask_perc_quantile,
        extents=extents, ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, 
        period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=None, period2_hours=None, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, 
        perc=perc, mask_perc_quantile=mask_perc_quantile,
        extents=extents, ax=axes[1][1], cfv_data=cfv_data
    )
    
    # Remove all params until we have the remaining 6 we want to create subplots for.
    
    params_to_plot = copy.deepcopy(
        cf.vars_and_dvars_era5[var_or_dvar_type][var_or_dvar_layer])
    
    for param in ["u10", "v10", "ws10", "wv10", "u100", "v100", "nse", "vidmf", "du10", 
                  "dv10", "dws10", "dwv10", "du100", "dv100", "dnse", "dvidmf"]:
        try:
            params_to_plot.remove(param)
        except:
            pass
    
    # Plot these params in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, param in enumerate(params_to_plot):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_diff_plot(
            calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months, 
            period1_hours=period1_hours, period2_hours=period2_hours, 
            arg_extra="mean", glass_source_pref=None, var_or_dvar=param, 
            perc=perc, mask_perc_quantile=mask_perc_quantile,
            extents=extents, ax=axes[row][col], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_diff_hourly_means_given_var_or_dvar(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, glass_source_pref, var_or_dvar, hours_to_plot, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    extents=None, cfv_data=None, output=False
):
    
    """
    Create difference plot for the hourly mean values of var_or_dvar between
    two different periods.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
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
        hours_to_plot (str): Which hours (in local time) to display results for in 
            plot. Must be one of: ["0-5", "6-11", "12-17", "18-23"].
        perc (bool): Whether to plot the difference in values as a percentage of the
            value (magnitude if negative) in period1. This is used for the comp plots
            in the plot_funcs script. Must be one of: [True, False].
        mask_perc_quantile (int): If perc is True, specify the quantile of values
            (magnitude if negative) from period1 to mask for the difference plot.
            This is used because percentage differences may be particularly high
            for values which had a low magnitude as a base in period 1.
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/hourly_means_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_diff_
        {extents_used}_{period1_start}_{period1_end}_{period2_start}_{period2_end}_
        means-hourly_{period1_months_str}_{period2_months_str}_{var_or_dvar}_
        {hours_to_plot}_perc-{quantile_used}.png:
            Output PNG file in data_final folder for the difference plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_months_str} and
            {period2_months_str} are strings representing the list of selected months 
            to use as a subset in each period. {quantile_used} is equal to 
            mask_perc_quantile if perc = True, otherwise it is set to None.
    
    For each period, calculate the hourly mean values for a given variable or its change 
    since the previous hour for Local Time: 0000 to 0500, 0600 to 1100, 1200 to 1700, 
    or 1800-2300 (inclusive) at each grid cell, then plot the difference in results between 
    the two periods. The hourly mean values are computed over the period between 
    period1_start and period1_end, and period2_start and period2_end (inclusive), and 
    only using a subset of months within these periods (if period1_months and period2_months 
    respectively not "all" is specified). Also included for reference are subplots for 
    mean leaf area index (MLAI) and mean fraction of absorbed photosynthetically active 
    radiation (MFAPAR).
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    The mask applies only to scalar variables (i.e. quiver plots for vectors are not 
    masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end, 
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_months=period1_months, period2_months=period2_months, 
                  glass_source_pref=glass_source_pref, 
                  var_or_dvar=var_or_dvar, hours_to_plot=hours_to_plot,
                  perc=perc, mask_perc_quantile=mask_perc_quantile,
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month subsets.
    
    period1_months_str = cf.get_period_months_str(period_months=period1_months)
    period2_months_str = cf.get_period_months_str(period_months=period2_months)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/hourly_means_given_var_or_dvar/" +
                   f"{plot_funcs_ver}_{cfv_used}_diff_{extents_used}_" +
                   f"{period1_start}_{period1_end}_{period2_start}_{period2_end}_" +
                   f"means-hourly_{period1_months_str}_{period2_months_str}_" +
                   f"{var_or_dvar}_{hours_to_plot}_perc-{quantile_used}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, 
        period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=None, period2_hours=None, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, 
        perc=perc, mask_perc_quantile=mask_perc_quantile,
        extents=extents, ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, 
        period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=None, period2_hours=None, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, 
        perc=perc, mask_perc_quantile=mask_perc_quantile,
        extents=extents, ax=axes[1][1], cfv_data=cfv_data
    )
    
    # For each hour in hours_to_plot, open the relevant datasets then use them to
    # determine colourbar and colourbar extents. If the intermediate output data
    # files don't yet exist, then create them.
    
    hours = cf.hour_subsets[hours_to_plot]
    datasets = []
    for hour in cf.hour_subsets["all"]:
    
        path_diff = cf.get_path_for_calc_diff(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period1_start=period1_start, period1_end=period1_end,
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months,
            period1_hours=[hour], period2_hours=[hour], 
            glass_source_pref=None, var_or_dvar=var_or_dvar
        )
    
        # Use data outputted from an older version of the calc_funcs script. This is useful
        # for results which required computationally intensive processing. And can also be
        # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
        if cfv_data:
            path_diff = path_diff.replace(cf.calc_funcs_ver, cfv_data)
            if Path(path_diff).exists() == False:
                msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                           f"but could not find file: {path_diff}")
                logging.error(msg_cfv)
                print(msg_cfv)
                cf.remove_handlers_if_directly_executed(func_1up)
                raise Exception(msg_cfv)
    
        if Path(path_diff).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_diff}"
            logging.info(msg_open)
        else:
            cf.calc_diff(
                calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
                period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_months=period1_months, period2_months=period2_months, 
                period1_hours=[hour], period2_hours=[hour], 
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
    
        ds_hour = (xr.open_dataset(path_diff, engine = "netcdf4")
                   .sel(longitude=slice(extents[0], extents[1]), 
                        latitude=slice(extents[3], extents[2]))
                   .expand_dims({"hour": [hour]})
                  )
        datasets.append(ds_hour)
    
    ds_hour = xr.merge(datasets)
    
    if var_or_dvar in cf.params_vector:
        da_hour_u = ds_hour["mean_u"]
        da_hour_v = ds_hour["mean_v"]
        da_hour_mag = cf.get_magnitude(da_hour_u, da_hour_v)
        vmin = float(da_hour_mag.min())
        vmax = float(da_hour_mag.max())
    else:
        da_hour = ds_hour["mean"]
        if var_or_dvar in vars_pos_with_vmin_0:
            vmin = 0
            vmax = float(da_hour.max())
        elif var_or_dvar in vars_pos:
            vmin = float(da_hour.min())
            vmax = float(da_hour.max())
        else:
            min_da = float(da_hour.min())
            max_da = float(da_hour.max())
            vmin = min(-abs(min_da), -abs(max_da))
            vmax = -vmin
    
    if vmin != None:
        if math.isnan(vmin):
            vmin = None
    if vmax != None:
        if math.isnan(vmax):
            vmax = None
    
    # Plot results for these hours in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, hour in enumerate(hours):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_diff_plot(
            calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months, 
            period1_hours=[hour], period2_hours=[hour], arg_extra="mean", 
            glass_source_pref=None, var_or_dvar=var_or_dvar, 
            perc=perc, mask_perc_quantile=mask_perc_quantile, extents=extents, 
            vmin=vmin, vmax=vmax, ax=axes[row][col], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_diff_monthly_means_given_var_or_dvar(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_hours, period2_hours, glass_source_pref, var_or_dvar, months_to_plot, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    extents=None, cfv_data=None, output=False
):
    
    """
    Create difference plot for the monthly mean values of var_or_dvar between
    two different periods.
    
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
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
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
        months_to_plot (str): Which months to display results for in plot. Must be
            one of: ["1-6", "7-12"].
        perc (bool): Whether to plot the difference in values as a percentage of the
            value (magnitude if negative) in period1. This is used for the comp plots
            in the plot_funcs script. Must be one of: [True, False].
        mask_perc_quantile (int): If perc is True, specify the quantile of values
            (magnitude if negative) from period1 to mask for the difference plot.
            This is used because percentage differences may be particularly high
            for values which had a low magnitude as a base in period 1.
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/monthly_means_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_diff_
        {extents_used}_{period1_start}_{period1_end}_{period2_start}_{period2_end}_
        means-monthly_{period1_hours_str}_{period2_hours_str}_{var_or_dvar}_
        {months_to_plot}_perc-{quantile_used}.png:
            Output PNG file in data_final folder for the difference plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_hours_str} and
            {period2_hours_str} are strings representing the list of selected hours 
            to use as a subset in each period. {quantile_used} is equal to 
            mask_perc_quantile if perc = True, otherwise it is set to None.
    
    For each period, calculate the monthly mean values for a given variable or its change 
    since the previous hour for months: Jan-Jun (1-6) or Jul-Dec (7-12) inclusive at each 
    grid cell, then plot the difference in results between the two periods. The monthly 
    mean values are computed over the period between period1_start and period1_end, and 
    period2_start and period2_end (inclusive), and only using a subset of hours within 
    these periods (if period1_hours and period2_hours respectively not "all" is specified). 
    Also included for reference are subplots for mean leaf area index (MLAI) and mean 
    fraction of absorbed photosynthetically active radiation (MFAPAR).
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    The mask applies only to scalar variables (i.e. quiver plots for vectors are not 
    masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end, 
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_hours=period1_hours, period2_hours=period2_hours, 
                  glass_source_pref=glass_source_pref, 
                  var_or_dvar=var_or_dvar, months_to_plot=months_to_plot,
                  perc=perc, mask_perc_quantile=mask_perc_quantile,
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for hour subsets.
    
    period1_hours_str = cf.get_period_hours_str(period_hours=period1_hours)
    period2_hours_str = cf.get_period_hours_str(period_hours=period2_hours)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/monthly_means_given_var_or_dvar/" +
                   f"{plot_funcs_ver}_{cfv_used}_diff_{extents_used}_" +
                   f"{period1_start}_{period1_end}_{period2_start}_{period2_end}_" +
                   f"means-monthly_{period1_hours_str}_{period2_hours_str}_" +
                   f"{var_or_dvar}_{months_to_plot}_perc-{quantile_used}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, 
        period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end, 
        period1_months="all", period2_months="all", 
        period1_hours=None, period2_hours=None, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, 
        perc=perc, mask_perc_quantile=mask_perc_quantile,
        extents=extents, ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, 
        period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end, 
        period1_months="all", period2_months="all", 
        period1_hours=None, period2_hours=None, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, 
        perc=perc, mask_perc_quantile=mask_perc_quantile,
        extents=extents, ax=axes[1][1], cfv_data=cfv_data
    )
    
    # For each month in months_to_plot, open the relevant datasets then use them to
    # determine colourbar and colourbar extents. If the intermediate output data
    # files don't yet exist, then create them.
    
    months = cf.month_subsets[months_to_plot]
    datasets = []
    for month in cf.month_subsets["all"]:
    
        path_diff = cf.get_path_for_calc_diff(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period1_start=period1_start, period1_end=period1_end,
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=[month], period2_months=[month],
            period1_hours=period1_hours, period2_hours=period2_hours, 
            glass_source_pref=None, var_or_dvar=var_or_dvar
        )
    
        # Use data outputted from an older version of the calc_funcs script. This is useful
        # for results which required computationally intensive processing. And can also be
        # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
        if cfv_data:
            path_diff = path_diff.replace(cf.calc_funcs_ver, cfv_data)
            if Path(path_diff).exists() == False:
                msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                           f"but could not find file: {path_diff}")
                logging.error(msg_cfv)
                print(msg_cfv)
                cf.remove_handlers_if_directly_executed(func_1up)
                raise Exception(msg_cfv)
    
        if Path(path_diff).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_diff}"
            logging.info(msg_open)
        else:
            cf.calc_diff(
                calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
                period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_months=[month], period2_months=[month],
                period1_hours=period1_hours, period2_hours=period2_hours, 
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
    
        ds_month = (xr.open_dataset(path_diff, engine = "netcdf4")
                    .sel(longitude=slice(extents[0], extents[1]), 
                         latitude=slice(extents[3], extents[2]))
                    .expand_dims({"month": [month]})
                   )
        datasets.append(ds_month)
    
    ds_month = xr.merge(datasets)
    
    if var_or_dvar in cf.params_vector:
        da_month_u = ds_month["mean_u"]
        da_month_v = ds_month["mean_v"]
        da_month_mag = cf.get_magnitude(da_month_u, da_month_v)
        vmin = float(da_month_mag.min())
        vmax = float(da_month_mag.max())
    else:
        da_month = ds_month["mean"]
        if var_or_dvar in vars_pos_with_vmin_0:
            vmin = 0
            vmax = float(da_month.max())
        elif var_or_dvar in vars_pos:
            vmin = float(da_month.min())
            vmax = float(da_month.max())
        else:
            min_da = float(da_month.min())
            max_da = float(da_month.max())
            vmin = min(-abs(min_da), -abs(max_da))
            vmax = -vmin
    
    if vmin != None:
        if math.isnan(vmin):
            vmin = None
    if vmax != None:
        if math.isnan(vmax):
            vmax = None
    
    # Plot results for these months in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, month in enumerate(months):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_diff_plot(
            calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=[month], period2_months=[month],
            period1_hours=period1_hours, period2_hours=period2_hours, arg_extra="mean", 
            glass_source_pref=None, var_or_dvar=var_or_dvar, 
            perc=perc, mask_perc_quantile=mask_perc_quantile, extents=extents, 
            vmin=vmin, vmax=vmax, ax=axes[row][col], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_diff_wsd_clim(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, period1_hours, period2_hours, glass_source_pref, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    extents=None, cfv_data=None, output=False
):
    
    """
    Create difference plot for the 100 m wind speed distribution parameters.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        perc (bool): Whether to plot the difference in values as a percentage of the
            value (magnitude if negative) in period1. This is used for the comp plots
            in the plot_funcs script. Must be one of: [True, False].
        mask_perc_quantile (int): If perc is True, specify the quantile of values
            (magnitude if negative) from period1 to mask for the difference plot.
            This is used because percentage differences may be particularly high
            for values which had a low magnitude as a base in period 1.
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/wsd_clim/{plot_funcs_ver}_{cfv_used}_diff_{extents_used}_
        {period1_start}_{period1_end}_{period2_start}_{period2_end}_
        {period1_months_str}_{period2_months_str}_wsd_{period1_hours_str}_
        {period2_hours_str}_perc-{quantile_used}.png:
            Output PNG file in data_final folder for the difference plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_months_str} and
            {period2_months_str} are strings representing the list of selected months 
            to use as a subset in each period. {period1_hours_str} and {period2_hours_str} 
            are strings representing the list of selected hours to use as a subset in 
            each period. {quantile_used} is equal to mask_perc_quantile if perc = True, 
            otherwise it is set to None.
    
    For each period, calculate the mean and standard deviation of wind speed at 100 m above 
    surface, the Weibull scale and shape parameter for wind speed at 100 m above surface 
    (C100 and K100), the expected rate of exceedance for a particular wind speed at 100 m 
    above surface (EROE100) and the gross capacity factor for a typical wind turbine at 
    100 m above surface (TGCF100) for each grid cell, then plot the difference in results 
    between the two periods. The wind speed distributions (WSDs) are computed over the 
    period between period1_start and period1_end, and period2_start and period2_end 
    (inclusive), using a subset of months (if period1_months and period2_months 
    respectively not "all" is  specified) and hours (if period1_hours and period2_hours 
    respectively not "all" is specified) within each period. Also included for reference 
    are subplots for mean leaf area index (MLAI) and mean fraction of absorbed 
    photosynthetically active radiation (MFAPAR).
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    The mask applies only to scalar variables (i.e. quiver plots for vectors are not 
    masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end,
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_months=period1_months, period2_months=period2_months, 
                  period1_hours=period1_hours, period2_hours=period2_hours,
                  glass_source_pref=glass_source_pref, 
                  perc=perc, mask_perc_quantile=mask_perc_quantile, 
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month and hour subsets.
    
    period1_months_str = cf.get_period_months_str(period_months=period1_months)
    period2_months_str = cf.get_period_months_str(period_months=period2_months)
    period1_hours_str = cf.get_period_hours_str(period_hours=period1_hours)
    period2_hours_str = cf.get_period_hours_str(period_hours=period2_hours)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/wsd_clim/{plot_funcs_ver}_{cfv_used}_" +
                   f"diff_{extents_used}_{period1_start}_{period1_end}_" +
                   f"{period2_start}_{period2_end}_{period1_months_str}_" +
                   f"{period2_months_str}_wsd_{period1_hours_str}_" +
                   f"{period2_hours_str}_perc-{quantile_used}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 5
    figcols = 2
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for orography and glass mean variables in first 2 rows.
    
    create_orog_static_plot(param_orog="lse", region=region, extents=extents, 
                            ax=axes[0][0], cfv_data=cfv_data)
    create_orog_static_plot(param_orog="ssgo", region=region, extents=extents, 
                            ax=axes[0][1], cfv_data=cfv_data)
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=period1_hours, period2_hours=period2_hours, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, extents=extents, 
        ax=axes[1][0], cfv_data=cfv_data
    )
    create_individual_diff_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=period1_hours, period2_hours=period2_hours, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, extents=extents, 
        ax=axes[1][1], cfv_data=cfv_data
    )
    
    # Remove all params until we have the remaining 6 we want to create subplots for.
    
    params_to_plot = copy.deepcopy(cf.params_wsd)
    
    for param in ["ws10_mean", "ws10_std", "c10", "k10"]:
        try:
            params_to_plot.remove(param)
        except:
            pass
    
    # Plot these params in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, param in enumerate(params_to_plot):
        row = math.floor(idx / figcols) + rows_to_skip
        col = idx % figcols
        create_individual_diff_plot(
            calc_func=cf.calc_era5_wsd_clim, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months, 
            period1_hours=period1_hours, period2_hours=period2_hours,
            arg_extra=param, glass_source_pref=None, var_or_dvar=None, 
            perc=perc, mask_perc_quantile=mask_perc_quantile, extents=extents, 
            ax=axes[row][col], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# ### Comp plots

# In[ ]:


def plot_comp_mdp_clim_stats_given_var_or_dvar(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, glass_source_pref, var_or_dvar, perc=False, 
    mask_perc_quantile=mask_perc_quantile_default, mask_period1=None, 
    mask_period2=None, extents=None, cfv_data=None, output=False
):
    
    """
    Create comparison plot for the mean diurnal profile (MDP) stats between
    two different periods.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
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
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/mdp_clim_stats_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_comp_
        {extents_used}_{period1_start}_{period1_end}_{period2_start}_{period2_end}_
        {period1_months_str}_{period2_months_str}_stats_{var_or_dvar}_
        perc-{quantile_used}_mask1-{mask_period1}_mask2-{mask_period2}.png:
            Output PNG file in data_final folder for the comparison plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_months_str} and
            {period2_months_str} are strings representing the list of selected months 
            to use as a subset in each period. {quantile_used} is equal to 
            mask_perc_quantile if perc = True, otherwise it is set to None.
    
    For each period, plot the max, min, mean, range, hour of max (24-hour Local Time) 
    and hour of min (24-hour Local Time), for the mean diurnal profile (MDP) of a given 
    variable for each grid cell, as well as the difference in results between the two 
    periods. The MDPs are computed over the period between period1_start and period1_end, 
    and period2_start and period2_end (inclusive), and only using a subset of data within 
    these periods (if period1_months and period2_months respectively not "all" is 
    specified). Also included for reference are subplots for mean leaf area index (MLAI) 
    and mean fraction of absorbed photosynthetically active radiation (MFAPAR).
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    If mask_period1 = "pos" or "neg" is specified, all grid cells where var_or_dvar was 
    positive or negative respectively in period1 will be masked, but only if var_or_dvar 
    was a parameter which could take both positive and negative values in the first place 
    (for positive-only parameters, this mask argument is ignored). Similar comments
    apply for mask_period2 but for values in period2. All 3 masks apply only to scalar
    variables (i.e. quiver plots for vectors are not masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end, 
                  period2_start=period2_start, period2_end=period2_end,
                  period1_months=period1_months, period2_months=period2_months, 
                  glass_source_pref=glass_source_pref,
                  var_or_dvar=var_or_dvar, perc=perc, 
                  mask_perc_quantile=mask_perc_quantile, 
                  mask_period1=mask_period1, mask_period2=mask_period2, 
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month subsets.
    
    period1_months_str = cf.get_period_months_str(period_months=period1_months)
    period2_months_str = cf.get_period_months_str(period_months=period2_months)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/mdp_clim_stats_given_var_or_dvar/" +
                   f"{plot_funcs_ver}_{cfv_used}_comp_{extents_used}_" +
                   f"{period1_start}_{period1_end}_{period2_start}_" +
                   f"{period2_end}_{period1_months_str}_{period2_months_str}_" +
                   f"stats_{var_or_dvar}_perc-{quantile_used}_" +
                   f"mask1-{mask_period1}_mask2-{mask_period2}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 8
    figcols = 3
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for glass mean variables in first 2 rows.
    
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[0][0], 
        ax_period2=axes[0][1], ax_diff=axes[0][2], cfv_data=cfv_data
    )
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[1][0], 
        ax_period2=axes[1][1], ax_diff=axes[1][2], cfv_data=cfv_data
    )
    
    # Remove all stats until we have the remaining 6 we want to create subplots for.
    
    stats_to_plot = copy.deepcopy(cf.params_stat)
    
    for stat in ["max_u", "max_v", "min_u", "min_v", "mean_u", "mean_v"]:
        try:
            stats_to_plot.remove(stat)
        except:
            pass
    
    # Plot these stats in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, stat in enumerate(stats_to_plot):
        row = idx + rows_to_skip
        create_individual_comp_plot(
            calc_func=cf.calc_era5_mdp_clim_stats_given_var_or_dvar, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months, 
            arg_extra=stat, glass_source_pref=None, var_or_dvar=var_or_dvar, 
            perc=perc, mask_perc_quantile=mask_perc_quantile,
            mask_period1=mask_period1, mask_period2=mask_period2, extents=extents, 
            ax_period1=axes[row][0], ax_period2=axes[row][1], ax_diff=axes[row][2], 
            cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_comp_means_given_layer_and_type(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, period1_hours, period2_hours, 
    glass_source_pref, var_or_dvar_layer, var_or_dvar_type, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=False
):
    
    """
    Create comparison plot for the means between two different periods for an
    assortment of different parameters specified by parameter layer and type.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
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
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/means_given_layer_and_type/{plot_funcs_ver}_{cfv_used}_comp_
        {extents_used}_{period1_start}_{period1_end}_{period2_start}_{period2_end}_
        means_{period1_months_str}_{period2_months_str}_{period1_hours_str}_
        {period2_hours_str}_{var_or_dvar_layer}_{var_or_dvar_type}_
        perc-{quantile_used}_mask1-{mask_period1}_mask2-{mask_period2}.png:
            Output PNG file in data_final folder for the comparison plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_months_str} and
            {period2_months_str} are strings representing the list of selected months 
            to use as a subset in each period. {period1_hours_str} and {period2_hours_str} 
            are strings representing the list of selected hours to use as a subset in 
            each period. {quantile_used} is equal to mask_perc_quantile if perc = True, 
            otherwise it is set to None.
    
    For each period, plot the mean values, as well as their differences, for an assortment 
    of parameters which are categorised within the surface ("sfc"), atmosphere ("atm") 
    or cloud ("cld") layers, and are either variable values themselves or an hourly change 
    in the variable values. The means are computed over the period between period1_start 
    and period1_end, and period2_start and period2_end (inclusive), using a subset of months 
    (if period1_months and period2_months respectively not "all" is specified) and hours 
    (if period1_hours and period2_hours respectively not "all" is specified) within
    each period. Also included for reference are subplots for mean leaf area index (MLAI) 
    and mean fraction of absorbed photosynthetically active radiation (MFAPAR).
    
    If var_or_dvar_layer == “sfc” and var_or_dvar_type == “vars”, then the wind speed at 
    100 m above surface (WS100), wind velocity at 100 m above surface (WV100), mean sea 
    level pressure (MSLP), temperature at 2 m above surface (T2), surface latent heat flux 
    (SLHF) and surface sensible heat flux (SSHF) will be plotted. If var_or_dvar_layer == 
    “sfc” and var_or_dvar_type == “dvars”, then the change in each of these variables as 
    compared with the value in the previous hour (dWS100, dWV100, dMSLP, dT2, dSLHF, 
    dSSHF) will be plotted instead.
    
    If var_or_dvar_layer == “atm” and var_or_dvar_type == “vars”, then the vertical integral 
    of energy conversion (VIEC), vertical integral of potential, internal and latent energy 
    (VIPILE), vertical integral of kinetic energy (VIKE), total column cloud liquid water 
    (TCCLW), total column water vapour (TCWV) and net atmospheric condensation (NAC) will be 
    plotted. If var_or_dvar_layer == “atm” and var_or_dvar_type == “dvars”, then the change 
    in each of these variables as compared with the value in the previous hour (dVIEC, 
    dVIPILE, dVIKE, dTCCLW, dTCWV, dNAC) will be plotted instead.
    
    If var_or_dvar_layer == “cld” and var_or_dvar_type == “vars”, then the boundary layer
    height (BLH), forecast albedo (FA), cloud base height (CBH), total cloud cover (TCC),
    convective available potential energy (CAPE) and convective inhibition (CI) will be
    plotted. If var_or_dvar_layer == “cld” and var_or_dvar_type == “dvars”, then the change 
    in each of these variables as compared with the value in the previous hour (dBLH, dFA, 
    dCBH, dTCC, dCAPE, dCI) will be plotted instead.
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    If mask_period1 = "pos" or "neg" is specified, all grid cells where var_or_dvar was 
    positive or negative respectively in period1 will be masked, but only if var_or_dvar 
    was a parameter which could take both positive and negative values in the first place 
    (for positive-only parameters, this mask argument is ignored). Similar comments
    apply for mask_period2 but for values in period2. All 3 masks apply only to scalar
    variables (i.e. quiver plots for vectors are not masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end, 
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_months=period1_months, period2_months=period2_months, 
                  period1_hours=period1_hours, period2_hours=period2_hours,
                  glass_source_pref=glass_source_pref,
                  var_or_dvar_layer=var_or_dvar_layer, var_or_dvar_type=var_or_dvar_type, 
                  perc=perc, mask_perc_quantile=mask_perc_quantile,
                  mask_period1=mask_period1, mask_period2=mask_period2, 
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month and hour subsets.
    
    period1_months_str = cf.get_period_months_str(period_months=period1_months)
    period2_months_str = cf.get_period_months_str(period_months=period2_months)
    period1_hours_str = cf.get_period_hours_str(period_hours=period1_hours)
    period2_hours_str = cf.get_period_hours_str(period_hours=period2_hours)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/means_given_layer_and_type/" +
                   f"{plot_funcs_ver}_{cfv_used}_comp_{extents_used}_" +
                   f"{period1_start}_{period1_end}_{period2_start}_{period2_end}_" +
                   f"means_{period1_months_str}_{period2_months_str}_" +
                   f"{period1_hours_str}_{period2_hours_str}_" +
                   f"{var_or_dvar_layer}_{var_or_dvar_type}_perc-{quantile_used}_" +
                   f"mask1-{mask_period1}_mask2-{mask_period2}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 8
    figcols = 3
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for glass mean variables in first 2 rows.
    
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[0][0], 
        ax_period2=axes[0][1], ax_diff=axes[0][2], cfv_data=cfv_data
    )
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[1][0], 
        ax_period2=axes[1][1], ax_diff=axes[1][2], cfv_data=cfv_data
    )
    
    # Remove all params until we have the remaining 6 we want to create subplots for.
    
    params_to_plot = copy.deepcopy(
        cf.vars_and_dvars_era5[var_or_dvar_type][var_or_dvar_layer])
    
    for param in ["u10", "v10", "ws10", "wv10", "u100", "v100", "nse", "vidmf", "du10", 
                  "dv10", "dws10", "dwv10", "du100", "dv100", "dnse", "dvidmf"]:
        try:
            params_to_plot.remove(param)
        except:
            pass
    
    # Plot these params in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, param in enumerate(params_to_plot):
        row = idx + rows_to_skip
        create_individual_comp_plot(
            calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months, 
            period1_hours=period1_hours, period2_hours=period2_hours, 
            arg_extra="mean", glass_source_pref=None, var_or_dvar=param, 
            perc=perc, mask_perc_quantile=mask_perc_quantile,
            mask_period1=mask_period1, mask_period2=mask_period2, extents=extents, 
            ax_period1=axes[row][0], ax_period2=axes[row][1], ax_diff=axes[row][2], 
            cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_comp_hourly_means_given_var_or_dvar(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, glass_source_pref, var_or_dvar, hours_to_plot, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=False
):
    
    """
    Create comparison plot for the hourly mean values of var_or_dvar between
    two different periods.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
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
        hours_to_plot (str): Which hours (in local time) to display results for in 
            plot. Must be one of: ["0-5", "6-11", "12-17", "18-23"].
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
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/hourly_means_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_comp_
        {extents_used}_{period1_start}_{period1_end}_{period2_start}_{period2_end}_
        means-hourly_{period1_months_str}_{period2_months_str}_{var_or_dvar}_
        {hours_to_plot}_perc-{quantile_used}_mask1-{mask_period1}_
        mask2-{mask_period2}.png:
            Output PNG file in data_final folder for the comparison plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_months_str} and
            {period2_months_str} are strings representing the list of selected months 
            to use as a subset in each period. {quantile_used} is equal to 
            mask_perc_quantile if perc = True, otherwise it is set to None.
    
    For each period, plot the hourly mean values for a given variable or its change since 
    the previous hour for Local Time: 0000 to 0500, 0600 to 1100, 1200 to 1700, 
    or 1800-2300 (inclusive) at each grid cell, as well as the difference in results between 
    the two periods. The hourly mean values are computed over the period between 
    period1_start and period1_end, and period2_start and period2_end (inclusive), and 
    only using a subset of months within these periods (if period1_months and period2_months 
    respectively not "all" is specified). Also included for reference are subplots for 
    mean leaf area index (MLAI) and mean fraction of absorbed photosynthetically active 
    radiation (MFAPAR).
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    If mask_period1 = "pos" or "neg" is specified, all grid cells where var_or_dvar was 
    positive or negative respectively in period1 will be masked, but only if var_or_dvar 
    was a parameter which could take both positive and negative values in the first place 
    (for positive-only parameters, this mask argument is ignored). Similar comments
    apply for mask_period2 but for values in period2. All 3 masks apply only to scalar
    variables (i.e. quiver plots for vectors are not masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end, 
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_months=period1_months, period2_months=period2_months, 
                  glass_source_pref=glass_source_pref, 
                  var_or_dvar=var_or_dvar, hours_to_plot=hours_to_plot,
                  perc=perc, mask_perc_quantile=mask_perc_quantile,
                  mask_period1=mask_period1, mask_period2=mask_period2, 
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month subsets.
    
    period1_months_str = cf.get_period_months_str(period_months=period1_months)
    period2_months_str = cf.get_period_months_str(period_months=period2_months)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/hourly_means_given_var_or_dvar/" +
                   f"{plot_funcs_ver}_{cfv_used}_comp_{extents_used}_" +
                   f"{period1_start}_{period1_end}_{period2_start}_{period2_end}_" +
                   f"means-hourly_{period1_months_str}_{period2_months_str}_" +
                   f"{var_or_dvar}_{hours_to_plot}_perc-{quantile_used}_" +
                   f"mask1-{mask_period1}_mask2-{mask_period2}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 8
    figcols = 3
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for glass mean variables in first 2 rows.
    
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[0][0], 
        ax_period2=axes[0][1], ax_diff=axes[0][2], cfv_data=cfv_data
    )
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[1][0], 
        ax_period2=axes[1][1], ax_diff=axes[1][2], cfv_data=cfv_data
    )
    
    # For each hour in hours_to_plot, open the relevant datasets then use them to
    # determine colourbar and colourbar extents. If the intermediate output data
    # files don't yet exist, then create them.
    
    hours = cf.hour_subsets[hours_to_plot]
    datasets_period1 = []
    datasets_period2 = []
    datasets_diff = []
    for hour in cf.hour_subsets["all"]:
    
        path_period1 = cf.get_path_for_calc_func(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period_start=period1_start, period_end=period1_end, 
            period_months=period1_months, period_hours=[hour], 
            glass_source_pref=None, var_or_dvar=var_or_dvar
        )
        path_period2 = cf.get_path_for_calc_func(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period_start=period2_start, period_end=period2_end, 
            period_months=period2_months, period_hours=[hour], 
            glass_source_pref=None, var_or_dvar=var_or_dvar
        )
        path_diff = cf.get_path_for_calc_diff(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period1_start=period1_start, period1_end=period1_end,
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months,
            period1_hours=[hour], period2_hours=[hour], 
            glass_source_pref=None, var_or_dvar=var_or_dvar
        )
    
        # Use data outputted from an older version of the calc_funcs script. This is useful
        # for results which required computationally intensive processing. And can also be
        # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
        if cfv_data:
            for path in [path_period1, path_period2, path_diff]:
                path = path.replace(cf.calc_funcs_ver, cfv_data)
                if Path(path).exists() == False:
                    msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                               f"but could not find file: {path}")
                    logging.error(msg_cfv)
                    print(msg_cfv)
                    cf.remove_handlers_if_directly_executed(func_1up)
                    raise Exception(msg_cfv)
    
        if Path(path_period1).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_period1}"
            logging.info(msg_open)
        else:
            cf.calc_era5_mean_clim_given_var_or_dvar(
                region=region, period_start=period1_start, period_end=period1_end, 
                period_months=period1_months, period_hours=[hour], 
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
            
        ds_period1 = (xr.open_dataset(path_period1, engine = "netcdf4")
                      .sel(longitude=slice(extents[0], extents[1]), 
                           latitude=slice(extents[3], extents[2]))
                      .expand_dims({"hour": [hour]})
                     )
        datasets_period1.append(ds_period1)
        
        if Path(path_period2).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_period2}"
            logging.info(msg_open)
        else:
            cf.calc_era5_mean_clim_given_var_or_dvar(
                region=region, period_start=period2_start, period_end=period2_end, 
                period_months=period2_months, period_hours=[hour], 
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
            
        ds_period2 = (xr.open_dataset(path_period2, engine = "netcdf4")
                      .sel(longitude=slice(extents[0], extents[1]), 
                           latitude=slice(extents[3], extents[2]))
                      .expand_dims({"hour": [hour]})
                     )
        datasets_period2.append(ds_period2)
        
        if Path(path_diff).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_diff}"
            logging.info(msg_open)
        else:
            cf.calc_diff(
                calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
                period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_months=period1_months, period2_months=period2_months, 
                period1_hours=[hour], period2_hours=[hour], 
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
            
        ds_diff = (xr.open_dataset(path_diff, engine = "netcdf4")
                   .sel(longitude=slice(extents[0], extents[1]), 
                        latitude=slice(extents[3], extents[2]))
                   .expand_dims({"hour": [hour]})
                  )
        datasets_diff.append(ds_diff)
    
    ds_period1 = xr.merge(datasets_period1)
    ds_period2 = xr.merge(datasets_period2)
    ds_diff = xr.merge(datasets_diff)
    
    if var_or_dvar in cf.params_vector:
        da_period1_u = ds_period1["mean_u"]
        da_period1_v = ds_period1["mean_v"]
        da_period1_mag = cf.get_magnitude(da_period1_u, da_period1_v)
        da_period2_u = ds_period2["mean_u"]
        da_period2_v = ds_period2["mean_v"]
        da_period2_mag = cf.get_magnitude(da_period2_u, da_period2_v)
        da_diff_u = ds_diff["mean_u"]
        da_diff_v = ds_diff["mean_v"]
        da_diff_mag = cf.get_magnitude(da_diff_u, da_diff_v)
        vmin_periods = float(min(da_period1_mag.min(), da_period2_mag.min()))
        vmax_periods = float(max(da_period1_mag.max(), da_period2_mag.max()))
        vmin_diff = float(da_diff_mag.min())
        vmax_diff = float(da_diff_mag.max())
    else:
        da_period1 = ds_period1["mean"]
        da_period2 = ds_period2["mean"]
        da_diff = ds_diff["mean"]
        min_da_diff = float(da_diff.min())
        max_da_diff = float(da_diff.max())
        vmin_diff = min(-abs(min_da_diff), -abs(max_da_diff))
        vmax_diff = -vmin_diff
        if var_or_dvar in vars_pos_with_vmin_0:
            vmin_periods = 0
            vmax_periods = float(max(da_period1.max(), da_period2.max()))
        elif var_or_dvar in vars_pos:
            vmin_periods = float(min(da_period1.min(), da_period2.min()))
            vmax_periods = float(max(da_period1.max(), da_period2.max()))
        else:
            min_da_periods = float(min(da_period1.min(), da_period2.min()))
            max_da_periods = float(max(da_period1.max(), da_period2.max()))
            vmin_periods = min(-abs(min_da_periods), -abs(max_da_periods))
            vmax_periods = -vmin_periods
    
    if vmin_diff != None:
        if math.isnan(vmin_diff):
            vmin_diff = None
    if vmax_diff != None:
        if math.isnan(vmax_diff):
            vmax_diff = None
    if vmin_periods != None:
        if math.isnan(vmin_periods):
            vmin_periods = None
    if vmax_periods != None:
        if math.isnan(vmax_periods):
            vmax_periods = None
    
    # Plot results for these hours in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, hour in enumerate(hours):
        row = idx + rows_to_skip
        create_individual_comp_plot(
            calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months, 
            period1_hours=[hour], period2_hours=[hour], arg_extra="mean", 
            glass_source_pref=None, var_or_dvar=var_or_dvar, 
            perc=perc, mask_perc_quantile=mask_perc_quantile, 
            mask_period1=mask_period1, mask_period2=mask_period2, extents=extents, 
            vmin_periods=vmin_periods, vmax_periods=vmax_periods, vmin_diff=vmin_diff, 
            vmax_diff=vmax_diff, ax_period1=axes[row][0], ax_period2=axes[row][1], 
            ax_diff=axes[row][2], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_comp_monthly_means_given_var_or_dvar(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_hours, period2_hours, glass_source_pref, var_or_dvar, months_to_plot, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=False
):
    
    """
    Create comparison plot for the monthly mean values of var_or_dvar between
    two different periods.
    
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
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
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
        months_to_plot (str): Which months to display results for in plot. Must be
            one of: ["1-6", "7-12"].
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
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/monthly_means_given_var_or_dvar/{plot_funcs_ver}_{cfv_used}_comp_
        {extents_used}_{period1_start}_{period1_end}_{period2_start}_{period2_end}_
        means-monthly_{period1_hours_str}_{period2_hours_str}_{var_or_dvar}_
        {months_to_plot}_perc-{quantile_used}_mask1-{mask_period1}_
        mask2-{mask_period2}.png:
            Output PNG file in data_final folder for the comparison plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_hours_str} and
            {period2_hours_str} are strings representing the list of selected hours 
            to use as a subset in each period. {quantile_used} is equal to 
            mask_perc_quantile if perc = True, otherwise it is set to None.
    
    For each period, plot the monthly mean values for a given variable or its change since 
    the previous hour for months: Jan-Jun (1-6) or Jul-Dec (7-12) inclusive at each grid 
    cell, as well as the difference in results between the two periods. The monthly mean 
    values are computed over the period between period1_start and period1_end, and 
    period2_start and period2_end (inclusive), and only using a subset of hours within 
    these periods (if period1_hours and period2_hours respectively not "all" is specified). 
    Also included for reference are subplots for mean leaf area index (MLAI) and mean 
    fraction of absorbed photosynthetically active radiation (MFAPAR).
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    If mask_period1 = "pos" or "neg" is specified, all grid cells where var_or_dvar was 
    positive or negative respectively in period1 will be masked, but only if var_or_dvar 
    was a parameter which could take both positive and negative values in the first place 
    (for positive-only parameters, this mask argument is ignored). Similar comments
    apply for mask_period2 but for values in period2. All 3 masks apply only to scalar
    variables (i.e. quiver plots for vectors are not masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end, 
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_hours=period1_hours, period2_hours=period2_hours, 
                  glass_source_pref=glass_source_pref, 
                  var_or_dvar=var_or_dvar, months_to_plot=months_to_plot,
                  perc=perc, mask_perc_quantile=mask_perc_quantile,
                  mask_period1=mask_period1, mask_period2=mask_period2, 
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for hour subsets.
    
    period1_hours_str = cf.get_period_hours_str(period_hours=period1_hours)
    period2_hours_str = cf.get_period_hours_str(period_hours=period2_hours)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/monthly_means_given_var_or_dvar/" +
                   f"{plot_funcs_ver}_{cfv_used}_comp_{extents_used}_" +
                   f"{period1_start}_{period1_end}_{period2_start}_{period2_end}_" +
                   f"means-monthly_{period1_hours_str}_{period2_hours_str}_" +
                   f"{var_or_dvar}_{months_to_plot}_perc-{quantile_used}_" +
                   f"mask1-{mask_period1}_mask2-{mask_period2}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 8
    figcols = 3
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for glass mean variables in first 2 rows.
    
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months="all", period2_months="all", arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[0][0], 
        ax_period2=axes[0][1], ax_diff=axes[0][2], cfv_data=cfv_data
    )
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months="all", period2_months="all", arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[1][0], 
        ax_period2=axes[1][1], ax_diff=axes[1][2], cfv_data=cfv_data
    )
    
    # For each month in months_to_plot, open the relevant datasets then use them to
    # determine colourbar and colourbar extents. If the intermediate output data
    # files don't yet exist, then create them.
    
    months = cf.month_subsets[months_to_plot]
    datasets_period1 = []
    datasets_period2 = []
    datasets_diff = []
    for month in cf.month_subsets["all"]:
    
        path_period1 = cf.get_path_for_calc_func(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period_start=period1_start, period_end=period1_end, 
            period_months=[month], period_hours=period1_hours, 
            glass_source_pref=None, var_or_dvar=var_or_dvar
        )
        path_period2 = cf.get_path_for_calc_func(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period_start=period2_start, period_end=period2_end, 
            period_months=[month], period_hours=period2_hours, 
            glass_source_pref=None, var_or_dvar=var_or_dvar
        )
        path_diff = cf.get_path_for_calc_diff(
            calc_func_name="calc_era5_mean_clim_given_var_or_dvar", region=region, 
            period1_start=period1_start, period1_end=period1_end,
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=[month], period2_months=[month],
            period1_hours=period1_hours, period2_hours=period2_hours, 
            glass_source_pref=None, var_or_dvar=var_or_dvar
        )
    
        # Use data outputted from an older version of the calc_funcs script. This is useful
        # for results which required computationally intensive processing. And can also be
        # set as "cfv00" to analyse test results output from a calc_funcs ipynb notebook.
    
        if cfv_data:
            for path in [path_period1, path_period2, path_diff]:
                path = path.replace(cf.calc_funcs_ver, cfv_data)
                if Path(path).exists() == False:
                    msg_cfv = (f"TERMINATED: cfv_data = {cfv_data} was specified " +
                               f"but could not find file: {path}")
                    logging.error(msg_cfv)
                    print(msg_cfv)
                    cf.remove_handlers_if_directly_executed(func_1up)
                    raise Exception(msg_cfv)
    
        if Path(path_period1).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_period1}"
            logging.info(msg_open)
        else:
            cf.calc_era5_mean_clim_given_var_or_dvar(
                region=region, period_start=period1_start, period_end=period1_end, 
                period_months=[month], period_hours=period1_hours, 
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
            
        ds_period1 = (xr.open_dataset(path_period1, engine = "netcdf4")
                      .sel(longitude=slice(extents[0], extents[1]), 
                           latitude=slice(extents[3], extents[2]))
                      .expand_dims({"month": [month]})
                     )
        datasets_period1.append(ds_period1)
        
        if Path(path_period2).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_period2}"
            logging.info(msg_open)
        else:
            cf.calc_era5_mean_clim_given_var_or_dvar(
                region=region, period_start=period2_start, period_end=period2_end, 
                period_months=[month], period_hours=period2_hours, 
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
            
        ds_period2 = (xr.open_dataset(path_period2, engine = "netcdf4")
                      .sel(longitude=slice(extents[0], extents[1]), 
                           latitude=slice(extents[3], extents[2]))
                      .expand_dims({"month": [month]})
                     )
        datasets_period2.append(ds_period2)
        
        if Path(path_diff).exists():
            msg_open = f"Opening: existing file for use in {func_cur}: {path_diff}"
            logging.info(msg_open)
        else:
            cf.calc_diff(
                calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
                period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_months=[month], period2_months=[month],
                period1_hours=period1_hours, period2_hours=period2_hours, 
                glass_source_pref=None, var_or_dvar=var_or_dvar
            )
            
        ds_diff = (xr.open_dataset(path_diff, engine = "netcdf4")
                   .sel(longitude=slice(extents[0], extents[1]), 
                        latitude=slice(extents[3], extents[2]))
                   .expand_dims({"month": [month]})
                  )
        datasets_diff.append(ds_diff)
    
    ds_period1 = xr.merge(datasets_period1)
    ds_period2 = xr.merge(datasets_period2)
    ds_diff = xr.merge(datasets_diff)
    
    if var_or_dvar in cf.params_vector:
        da_period1_u = ds_period1["mean_u"]
        da_period1_v = ds_period1["mean_v"]
        da_period1_mag = cf.get_magnitude(da_period1_u, da_period1_v)
        da_period2_u = ds_period2["mean_u"]
        da_period2_v = ds_period2["mean_v"]
        da_period2_mag = cf.get_magnitude(da_period2_u, da_period2_v)
        da_diff_u = ds_diff["mean_u"]
        da_diff_v = ds_diff["mean_v"]
        da_diff_mag = cf.get_magnitude(da_diff_u, da_diff_v)
        vmin_periods = float(min(da_period1_mag.min(), da_period2_mag.min()))
        vmax_periods = float(max(da_period1_mag.max(), da_period2_mag.max()))
        vmin_diff = float(da_diff_mag.min())
        vmax_diff = float(da_diff_mag.max())
    else:
        da_period1 = ds_period1["mean"]
        da_period2 = ds_period2["mean"]
        da_diff = ds_diff["mean"]
        min_da_diff = float(da_diff.min())
        max_da_diff = float(da_diff.max())
        vmin_diff = min(-abs(min_da_diff), -abs(max_da_diff))
        vmax_diff = -vmin_diff
        if var_or_dvar in vars_pos_with_vmin_0:
            vmin_periods = 0
            vmax_periods = float(max(da_period1.max(), da_period2.max()))
        elif var_or_dvar in vars_pos:
            vmin_periods = float(min(da_period1.min(), da_period2.min()))
            vmax_periods = float(max(da_period1.max(), da_period2.max()))
        else:
            min_da_periods = float(min(da_period1.min(), da_period2.min()))
            max_da_periods = float(max(da_period1.max(), da_period2.max()))
            vmin_periods = min(-abs(min_da_periods), -abs(max_da_periods))
            vmax_periods = -vmin_periods
    
    if vmin_diff != None:
        if math.isnan(vmin_diff):
            vmin_diff = None
    if vmax_diff != None:
        if math.isnan(vmax_diff):
            vmax_diff = None
    if vmin_periods != None:
        if math.isnan(vmin_periods):
            vmin_periods = None
    if vmax_periods != None:
        if math.isnan(vmax_periods):
            vmax_periods = None
    
    # Plot results for these months in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, month in enumerate(months):
        row = idx + rows_to_skip
        create_individual_comp_plot(
            calc_func=cf.calc_era5_mean_clim_given_var_or_dvar, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=[month], period2_months=[month],
            period1_hours=period1_hours, period2_hours=period2_hours, arg_extra="mean", 
            glass_source_pref=None, var_or_dvar=var_or_dvar, 
            perc=perc, mask_perc_quantile=mask_perc_quantile, 
            mask_period1=mask_period1, mask_period2=mask_period2, extents=extents, 
            vmin_periods=vmin_periods, vmax_periods=vmax_periods, vmin_diff=vmin_diff, 
            vmax_diff=vmax_diff, ax_period1=axes[row][0], ax_period2=axes[row][1], 
            ax_diff=axes[row][2], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def plot_comp_wsd_clim(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, period1_hours, period2_hours, glass_source_pref, 
    perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=False
):
    
    """
    Create comparison plot for the 100 m wind speed distribution parameters.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
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
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    
    Returns (if output = True):
        ../data_final/wsd_clim/{plot_funcs_ver}_{cfv_used}_comp_{extents_used}_
        {period1_start}_{period1_end}_{period2_start}_{period2_end}_
        {period1_months_str}_{period2_months_str}_wsd_{period1_hours_str}_
        {period2_hours_str}_perc-{quantile_used}_mask1-{mask_period1}_
        mask2-{mask_period2}.png:
            Output PNG file in data_final folder for the comparison plot. {plot_funcs_ver}
            is the version of the plot_funcs script being used. {cfv_used} is the version
            of the calc_funcs script which outputted the data used in making this plot.
            {extents_used} is a string indicating the region, or the WESN coordinates
            for the extents argument if this was specified. {period1_months_str} and
            {period2_months_str} are strings representing the list of selected months 
            to use as a subset in each period. {period1_hours_str} and {period2_hours_str} 
            are strings representing the list of selected hours to use as a subset in 
            each period. {quantile_used} is equal to mask_perc_quantile if perc = True, 
            otherwise it is set to None.
    
    For each period, plot the mean and standard deviation of wind speed at 100 m above 
    surface, the Weibull scale and shape parameter for wind speed at 100 m above surface 
    (C100 and K100), the expected rate of exceedance for a particular wind speed at 100 m 
    above surface (EROE100) and the gross capacity factor for a typical wind turbine at 
    100 m above surface (TGCF100) for each grid cell, as well as the difference in results 
    between the two periods. The wind speed distributions (WSDs) are computed over the 
    period between period1_start and period1_end, and period2_start and period2_end 
    (inclusive), using a subset of months (if period1_months and period2_months 
    respectively not "all" is specified) and hours (if period1_hours and period2_hours 
    respectively not "all" is specified) within each period. Also included for reference 
    are subplots for mean leaf area index (MLAI) and mean fraction of absorbed 
    photosynthetically active radiation (MFAPAR).
    
    If perc = True is specified, the differences in results between periods are plotted
    as percentage deviations from the magnitude of var_or_dvar in period1 at each 
    grid cell. Extremely high percentages may arise if coming off a low base magnitude
    in period1, so the mask_perc_quantile argument is used to specify at what magnitude
    quantile should a grid cell be masked. Eg. if mask_perc_quantile = 10, all grid cells
    where the magnitude of var_or_dvar was within the lowest 10% of values will be masked.
    If mask_period1 = "pos" or "neg" is specified, all grid cells where var_or_dvar was 
    positive or negative respectively in period1 will be masked, but only if var_or_dvar 
    was a parameter which could take both positive and negative values in the first place 
    (for positive-only parameters, this mask argument is ignored). Similar comments
    apply for mask_period2 but for values in period2. All 3 masks apply only to scalar
    variables (i.e. quiver plots for vectors are not masked by their magnitude).
    
    The computation uses data from the data_raw folder, then outputs 
    intermediate results as a netcdf4 file into the data_processed folder, 
    and a PNG plot into the data_final folder.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    # Assert that input arguments are valid
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end,
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_months=period1_months, period2_months=period2_months, 
                  period1_hours=period1_hours, period2_hours=period2_hours,
                  glass_source_pref=glass_source_pref, 
                  perc=perc, mask_perc_quantile=mask_perc_quantile,
                  mask_period1=mask_period1, mask_period2=mask_period2,
                  extents=extents, cfv_data=cfv_data, output=output)
    
    # Obtain string representations for month and hour subsets.
    
    period1_months_str = cf.get_period_months_str(period_months=period1_months)
    period2_months_str = cf.get_period_months_str(period_months=period2_months)
    period1_hours_str = cf.get_period_hours_str(period_hours=period1_hours)
    period2_hours_str = cf.get_period_hours_str(period_hours=period2_hours)
    
    # Create copy of inputs to use for output path.
    
    extents_input = copy.deepcopy(extents)
    
    if extents == None:
        extents = copy.deepcopy(cf.regions[region]["extents"])
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
            
    if extents_input:
        extents_used = "{}W{}E{}S{}N".format(extents[0], extents[1], 
                                             extents[2], extents[3])
    else:
        extents_used = copy.deepcopy(region)

    if perc == True:
        quantile_used = None
    else:
        quantile_used = copy.deepcopy(mask_perc_quantile)
        
    path_output = (f"../data_final/wsd_clim/{plot_funcs_ver}_{cfv_used}_" +
                   f"comp_{extents_used}_{period1_start}_{period1_end}_" +
                   f"{period2_start}_{period2_end}_{period1_months_str}_" +
                   f"{period2_months_str}_wsd_{period1_hours_str}_" +
                   f"{period2_hours_str}_perc-{quantile_used}_" +
                   f"mask1-{mask_period1}_mask2-{mask_period2}.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        # Skip the plotting if output plot already exists and we are
        # calling this function from a top level function.
        if func_1up in funcs_create_all_plot:
            return None
    
    # Create plot.
    
    figrows = 8
    figcols = 3
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols * 
                 (extents[3]-extents[2])/(extents[1]-extents[0])
                )
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    # Create subplots for glass mean variables in first 2 rows.
    
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=period1_hours, period2_hours=period2_hours, arg_extra="mlai", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[0][0], 
        ax_period2=axes[0][1], ax_diff=axes[0][2], cfv_data=cfv_data
    )
    create_individual_comp_plot(
        calc_func=cf.calc_glass_mean_clim, region=region, period1_start=period1_start, 
        period1_end=period1_end, period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months, 
        period1_hours=period1_hours, period2_hours=period2_hours, arg_extra="mfapar", 
        glass_source_pref=glass_source_pref, var_or_dvar=None, perc=perc, 
        mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1, 
        mask_period2=mask_period2, extents=extents, ax_period1=axes[1][0], 
        ax_period2=axes[1][1], ax_diff=axes[1][2], cfv_data=cfv_data
    )
    
    # Remove all params until we have the remaining 6 we want to create subplots for.
    
    params_to_plot = copy.deepcopy(cf.params_wsd)
    
    for param in ["ws10_mean", "ws10_std", "c10", "k10"]:
        try:
            params_to_plot.remove(param)
        except:
            pass
    
    # Plot these params in next 6 rows.
    
    rows_to_skip = 2
    
    for idx, param in enumerate(params_to_plot):
        row = idx + rows_to_skip
        create_individual_comp_plot(
            calc_func=cf.calc_era5_wsd_clim, region=region, 
            period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end, 
            period1_months=period1_months, period2_months=period2_months, 
            period1_hours=period1_hours, period2_hours=period2_hours,
            arg_extra=param, glass_source_pref=None, var_or_dvar=None, perc=perc,
            mask_perc_quantile=mask_perc_quantile, mask_period1=mask_period1,
            mask_period2=mask_period2, extents=extents, ax_period1=axes[row][0],
            ax_period2=axes[row][1], ax_diff=axes[row][2], cfv_data=cfv_data
        )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    if func_1up not in funcs_create_all_plot:
        plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# ## Top-level plotting functions to create all possible plot files

# In[ ]:


def create_comb_orog_static_plot(cfv_data=None, output=False):
    
    """
    Create a single plot with 6 subplots displaying the land surface elevation
    and slope of sub-gridscale orography for the Western Australia ("wa"),
    Central America ("ca") and South America ("sa") study regions defined
    in the calc_funcs script.
    
    Arguments:
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
        output (bool): Whether to output the plot as a PNG file. Must be one of:
            [True, False].
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    if cfv_data:
        cfv_used = copy.deepcopy(cfv_data)
    else:
        cfv_used = copy.deepcopy(cf.calc_funcs_ver)
    
    path_output = (f"../data_final/era5_orog/" +
                   f"{plot_funcs_ver}_{cfv_used}_comb_era5-orog.png")
    
    if Path(path_output).exists():
        msg_exist = ("WARNING: plot file already exists (and was " +
                     f"not overwritten): {path_output}")
        logging.warning(msg_exist)
        print(msg_exist)
        
    figrows = 2
    figcols = 3
    figwidth = figwidth_standard * 2
    figheight = (figwidth * figrows/figcols)
    fig, axes = plt.subplots(figrows, figcols, figsize=(figwidth, figheight), 
                             subplot_kw = {"projection": ccrs.PlateCarree()}
                            )
    
    create_orog_static_plot(
        param_orog="lse", region="wa", extents=None, vmin=None, vmax=None, 
        ax=axes[0][0], cfv_data=cfv_data, output=False
    )
    create_orog_static_plot(
        param_orog="lse", region="ca", extents=None, vmin=None, vmax=None, 
        ax=axes[0][1], cfv_data=cfv_data, output=False
    )
    create_orog_static_plot(
        param_orog="lse", region="sa", extents=None, vmin=None, vmax=None, 
        ax=axes[0][2], cfv_data=cfv_data, output=False
    )
    create_orog_static_plot(
        param_orog="ssgo", region="wa", extents=None, vmin=None, vmax=None, 
        ax=axes[1][0], cfv_data=cfv_data, output=False
    )
    create_orog_static_plot(
        param_orog="ssgo", region="ca", extents=None, vmin=None, vmax=None, 
        ax=axes[1][1], cfv_data=cfv_data, output=False
    )
    create_orog_static_plot(
        param_orog="ssgo", region="sa", extents=None, vmin=None, vmax=None, 
        ax=axes[1][2], cfv_data=cfv_data, output=False
    )
    
    # for idx in range(0, figrows * figcols):
    #     row = math.floor(idx / figcols)
    #     col = idx % figcols
    #     ax_title = axes[row][col].get_title()
    #     axes[row][col].set_title(chr(ord('`')+(idx+1)) + ") " + ax_title)
    
    fig.tight_layout()
    
    if (output == True) & (Path(path_output).exists() == False):
        path_output_dir = "/".join(path_output.split("/")[:-1])
        Path(path_output_dir).mkdir(parents=True, exist_ok=True)
        plt.savefig(path_output, metadata=get_plot_metadata(time_exec, func_cur, 
                                                            args_cur, args_cur_values))
        msg_create = f"CREATED: plot file: {path_output}"
        logging.info(msg_create)
        print(msg_create)
    
    plt.show()
    fig.clear()
    plt.close(fig)
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def create_all_possible_calc_plot_files(
    region, period_start, period_end, period_months, period_hours, 
    glass_source_pref, extents=None, cfv_data=None
):
    
    """
    For the given inputs, run all possible high-level plot functions over a period.
    
    Arguments:
        region (str): Region to perform calculation over.
            Must be one of ["ca", "sa", "wa"].
        period_start (str): Start of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_end (str): End of period to perform calculation over.
            Must be of form "%b-%Y" eg. "Jul-1990".
            Must be between "Jan-1981" and "Dec-2021".
        period_months (str or list): Months subset of period to perform calculation over.
            Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a subset
            list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period_hours (str or list): Hours subset of period to perform calculation over.
            Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period_start=period_start, period_end=period_end, 
                  period_months=period_months, period_hours=period_hours, 
                  glass_source_pref=glass_source_pref, 
                  extents=extents, cfv_data=cfv_data)
    
    vars_and_dvars = copy.deepcopy(cf.vars_and_dvars_era5_all)
    # Wind speeds require hourly data so we remove these from global analysis
    # (which only uses monthly averaged by hour of day data due to memory limits)
    if region == "global":
        for var_or_dvar in ["ws10", "ws100", "dws10", "dws100"]:
            vars_and_dvars.remove(var_or_dvar)
    
    for var_or_dvar in vars_and_dvars:
        plot_mdp_clim_stats_given_var_or_dvar(
            region=region, period_start=period_start, period_end=period_end, 
            period_months=period_months, glass_source_pref=glass_source_pref, 
            var_or_dvar=var_or_dvar, extents=extents, cfv_data=cfv_data, output=True
        )
        
        for hours in cf.hours_to_plot_valid:
            plot_hourly_means_given_var_or_dvar(
                region=region, period_start=period_start, period_end=period_end, 
                period_months=period_months, glass_source_pref=glass_source_pref, 
                var_or_dvar=var_or_dvar, hours_to_plot=hours,
                extents=extents, cfv_data=cfv_data, output=True
            )
        
        for months in cf.months_to_plot_valid:
            plot_monthly_means_given_var_or_dvar(
                region=region, period_start=period_start, period_end=period_end, 
                period_hours=period_hours, glass_source_pref=glass_source_pref, 
                var_or_dvar=var_or_dvar, months_to_plot=months,
                extents=extents, cfv_data=cfv_data, output=True
            )
    
    for var_or_dvar_layer in cf.var_or_dvar_layers:
        for var_or_dvar_type in cf.var_or_dvar_types:
            plot_means_given_layer_and_type(
                region=region, period_start=period_start, period_end=period_end, 
                period_months=period_months, period_hours=period_hours,
                glass_source_pref=glass_source_pref, 
                var_or_dvar_layer=var_or_dvar_layer, var_or_dvar_type=var_or_dvar_type, 
                extents=extents, cfv_data=cfv_data, output=True
            )
    
    if region == "global":
        cf.remove_handlers_if_directly_executed(func_1up)
        return None
    
    plot_wsd_clim(
        region=region, period_start=period_start, period_end=period_end, 
        period_months=period_months, period_hours=period_hours,
        glass_source_pref=glass_source_pref,  
        extents=extents, cfv_data=cfv_data, output=True
    )
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def create_all_possible_diff_plot_files(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, period1_hours, period2_hours, 
    glass_source_pref, perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    extents=None, cfv_data=None
):
    
    """
    For the given inputs, run all possible high-level difference plot functions 
    for the difference in results between two periods.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
        perc (bool): Whether to plot the difference in values as a percentage of the
            value (magnitude if negative) in period1. This is used for the comp plots
            in the plot_funcs script. Must be one of: [True, False].
        mask_perc_quantile (int): If perc is True, specify the quantile of values
            (magnitude if negative) from period1 to mask for the difference plot.
            This is used because percentage differences may be particularly high
            for values which had a low magnitude as a base in period 1.
        extents (list): Longitudinal and latitudinal extents to display in plot.
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end,
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_months=period1_months, period2_months=period2_months, 
                  period1_hours=period1_hours, period2_hours=period2_hours, 
                  glass_source_pref=glass_source_pref, 
                  perc=perc, mask_perc_quantile=mask_perc_quantile, extents=extents, 
                  cfv_data=cfv_data)
    
    vars_and_dvars = copy.deepcopy(cf.vars_and_dvars_era5_all)
    # Wind speeds require hourly data so we remove these from global analysis
    # (which only uses monthly averaged by hour of day data due to memory limits)
    if region == "global":
        for var_or_dvar in ["ws10", "ws100", "dws10", "dws100"]:
            vars_and_dvars.remove(var_or_dvar)
    
    for var_or_dvar in vars_and_dvars:
        plot_diff_mdp_clim_stats_given_var_or_dvar(
            region=region, period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end,
            period1_months=period1_months, period2_months=period2_months, 
            glass_source_pref=glass_source_pref, 
            var_or_dvar=var_or_dvar, perc=perc, mask_perc_quantile=mask_perc_quantile, 
            extents=extents, cfv_data=cfv_data, output=True
        )
        
        for hours in cf.hours_to_plot_valid:
            plot_diff_hourly_means_given_var_or_dvar(
                region=region, period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_months=period1_months, period2_months=period2_months, 
                glass_source_pref=glass_source_pref, 
                var_or_dvar=var_or_dvar, hours_to_plot=hours, 
                perc=perc, mask_perc_quantile=mask_perc_quantile,
                extents=extents, cfv_data=cfv_data, output=True
            )
        
        for months in cf.months_to_plot_valid:
            plot_diff_monthly_means_given_var_or_dvar(
                region=region, period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_hours=period1_hours, period2_hours=period2_hours, 
                glass_source_pref=glass_source_pref, 
                var_or_dvar=var_or_dvar, months_to_plot=months, 
                perc=perc, mask_perc_quantile=mask_perc_quantile,
                extents=extents, cfv_data=cfv_data, output=True
            )
    
    for var_or_dvar_layer in cf.var_or_dvar_layers:
        for var_or_dvar_type in cf.var_or_dvar_types:
            plot_diff_means_given_layer_and_type(
                region=region, period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_months=period1_months, period2_months=period2_months, 
                period1_hours=period1_hours, period2_hours=period2_hours, 
                glass_source_pref=glass_source_pref, 
                var_or_dvar_layer=var_or_dvar_layer, var_or_dvar_type=var_or_dvar_type, 
                perc=perc, mask_perc_quantile=mask_perc_quantile,
                extents=extents, cfv_data=cfv_data, output=True
            )
    
    if region == "global":
        cf.remove_handlers_if_directly_executed(func_1up)
        return None
    
    plot_diff_wsd_clim(
        region=region, period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months,
        period1_hours=period1_hours, period2_hours=period2_hours, 
        glass_source_pref=glass_source_pref,  
        perc=perc, mask_perc_quantile=mask_perc_quantile, 
        extents=extents, cfv_data=cfv_data, output=True
    )
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:


def create_all_possible_comp_plot_files(
    region, period1_start, period1_end, period2_start, period2_end, 
    period1_months, period2_months, period1_hours, period2_hours, 
    glass_source_pref, perc=False, mask_perc_quantile=mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None
):
    
    """
    For the given inputs, run all possible high-level comparison plot functions 
    displaying the results over each periods as well as the difference in results 
    between the two periods.
    
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
        period1_months (str or list): Month subset of first period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period2_months (str or list): Month subset of second period to perform calculation 
            over. Must be a str and one of: ["all", "djf", "mam", "jja", "son"], or a 
            subset list of: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12] with at least one item.
        period1_hours (str or list): Hours subset of first period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        period2_hours (str or list): Hours subset of second period to perform calculation 
            over. Must be a str and one of: ["0-5"/"night", "6-11"/"morning", 
            "12-17"/"afternoon", "18-23"/"evening", "all"], or a subset
            list of: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 
            15, 16, 17, 18, 19, 20, 21, 22, 23] with at least one item.
            Should be expressed in local time corresponding to selected region.
        glass_source_pref (str): Preferred glass data source to use when analysis is 
            over a period which is completely contained within both the available
            AVHRR and MODIS datasets. Must be one of: ["avhrr", "modis"].
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
            Must be a 4 element list in [W, E, S, N] format with longitudes
            between -180 to 180 and latitudes between -90 to 90.
        cfv_data (str): calc_funcs_ver of pre-existing data to use in plotting.
    """
    
    time_exec = datetime.today()
    func_cur = inspect.stack()[0][3]
    func_1up = inspect.stack()[1][3]
    frame_cur = inspect.currentframe()
    args_cur, _, _, args_cur_values = inspect.getargvalues(frame_cur)
    cf.create_log_if_directly_executed(time_exec, func_cur, func_1up, 
                                       args_cur, args_cur_values)
    
    cf.check_args_for_none(func_cur, args_cur, args_cur_values)
    cf.check_args(region=region, period1_start=period1_start, period1_end=period1_end,
                  period2_start=period2_start, period2_end=period2_end, 
                  period1_months=period1_months, period2_months=period2_months, 
                  period1_hours=period1_hours, period2_hours=period2_hours, 
                  glass_source_pref=glass_source_pref, 
                  perc=perc, mask_perc_quantile=mask_perc_quantile, 
                  mask_period1=mask_period1, mask_period2=mask_period2, 
                  extents=extents, cfv_data=cfv_data)
    
    vars_and_dvars = copy.deepcopy(cf.vars_and_dvars_era5_all)
    # Wind speeds require hourly data so we remove these from global analysis
    # (which only uses monthly averaged by hour of day data due to memory limits)
    if region == "global":
        for var_or_dvar in ["ws10", "ws100", "dws10", "dws100"]:
            vars_and_dvars.remove(var_or_dvar)
    
    for var_or_dvar in vars_and_dvars:
        plot_comp_mdp_clim_stats_given_var_or_dvar(
            region=region, period1_start=period1_start, period1_end=period1_end, 
            period2_start=period2_start, period2_end=period2_end,
            period1_months=period1_months, period2_months=period2_months, 
            glass_source_pref=glass_source_pref, 
            var_or_dvar=var_or_dvar, perc=perc, mask_perc_quantile=mask_perc_quantile, 
            mask_period1=mask_period1, mask_period2=mask_period2,
            extents=extents, cfv_data=cfv_data, output=True
        )
        
        for hours in cf.hours_to_plot_valid:
            plot_comp_hourly_means_given_var_or_dvar(
                region=region, period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_months=period1_months, period2_months=period2_months, 
                glass_source_pref=glass_source_pref, 
                var_or_dvar=var_or_dvar, hours_to_plot=hours, 
                perc=perc, mask_perc_quantile=mask_perc_quantile,
                mask_period1=mask_period1, mask_period2=mask_period2,
                extents=extents, cfv_data=cfv_data, output=True
            )
        
        for months in cf.months_to_plot_valid:
            plot_comp_monthly_means_given_var_or_dvar(
                region=region, period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_hours=period1_hours, period2_hours=period2_hours, 
                glass_source_pref=glass_source_pref, 
                var_or_dvar=var_or_dvar, months_to_plot=months, 
                perc=perc, mask_perc_quantile=mask_perc_quantile,
                mask_period1=mask_period1, mask_period2=mask_period2,
                extents=extents, cfv_data=cfv_data, output=True
            )
    
    for var_or_dvar_layer in cf.var_or_dvar_layers:
        for var_or_dvar_type in cf.var_or_dvar_types:
            plot_comp_means_given_layer_and_type(
                region=region, period1_start=period1_start, period1_end=period1_end, 
                period2_start=period2_start, period2_end=period2_end, 
                period1_months=period1_months, period2_months=period2_months, 
                period1_hours=period1_hours, period2_hours=period2_hours, 
                glass_source_pref=glass_source_pref, 
                var_or_dvar_layer=var_or_dvar_layer, var_or_dvar_type=var_or_dvar_type, 
                perc=perc, mask_perc_quantile=mask_perc_quantile,
                mask_period1=mask_period1, mask_period2=mask_period2,
                extents=extents, cfv_data=cfv_data, output=True
            )
    
    if region == "global":
        cf.remove_handlers_if_directly_executed(func_1up)
        return None
    
    plot_comp_wsd_clim(
        region=region, period1_start=period1_start, period1_end=period1_end, 
        period2_start=period2_start, period2_end=period2_end, 
        period1_months=period1_months, period2_months=period2_months,
        period1_hours=period1_hours, period2_hours=period2_hours, 
        glass_source_pref=glass_source_pref,  
        perc=perc, mask_perc_quantile=mask_perc_quantile, 
        mask_period1=mask_period1, mask_period2=mask_period2,
        extents=extents, cfv_data=cfv_data, output=True
    )
    
    cf.remove_handlers_if_directly_executed(func_1up)


# In[ ]:





# # Scrap

# In[ ]:




