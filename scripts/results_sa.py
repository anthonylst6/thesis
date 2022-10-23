#!/usr/bin/env python
# coding: utf-8

# # Setup
# 
# Note: the following code requires up to 30 GB of RAM over 2 hours if using 4 CPUs

# ## Import libraries for analysis

# In[ ]:


import importlib
from glob import glob
from datetime import datetime


# ## Import calc_funcs module for use in analysis

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


# ## Import plot_funcs module for use in analysis

# In[ ]:


# Choose plot_funcs_ver to use.
pfv = "latest"

assert (isinstance(pfv, str) & (len(pfv) == 5) & (pfv[:3] == "pfv") & 
        pfv[3].isnumeric() & pfv[4].isalpha() & pfv[4].islower()) | (pfv == "latest"), \
    ("pfv must be 'latest' or of form 'pfvXY' where X is a single digit number " +
     "and Y is a lowercase alphabet character. eg. pfv1e")

if pfv == "latest":
    plot_funcs_scripts = glob("plot_funcs_*.py")
    plot_funcs_scripts.sort()
    plot_funcs_module = plot_funcs_scripts[-1][:-3]
    
else:
    plot_funcs_module = "plot_funcs_" + pfv[2:]
    
pf = importlib.import_module(plot_funcs_module)

print(f"Using: {plot_funcs_module}")


# ## Analysis settings using selected periods from results below
# 
# Specify these all here in the beginning so errors from changes are minimised

# In[ ]:


# General analysis settings
region = "sa"
year_start = 1983
year_end = 2019
window_size = 5
months_wet = [1, 2, 3, 4, 5, 6]
months_dry = [7, 8, 9, 10, 11, 12]

# Periods with similar background atmospheric oscillations
period1_mid_sim = "Jul-2005"
period1_start_sim = "Jan-2003"
period1_end_sim = "Dec-2007"
period2_mid_sim = "Jul-2018"
period2_start_sim = "Jan-2016"
period2_end_sim = "Dec-2020"
if (datetime.strptime(period1_start_sim, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_sim = "modis"
else:
    glass_source_pref_sim = "avhrr"
    
# Periods with dissimilar background atmospheric oscillations
period1_mid_dis = "Mar-1985"
period1_start_dis = "Sep-1982"
period1_end_dis = "Aug-1987"
period2_mid_dis = "Jul-2010"
period2_start_dis = "Jan-2008"
period2_end_dis = "Dec-2012"
if (datetime.strptime(period1_start_dis, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_dis = "modis"
else:
    glass_source_pref_dis = "avhrr"

# Single months with very dissimilar background atmospheric oscillations
month1_dis = "Nov-1996"
month2_dis = "Sep-2015"
if (datetime.strptime(month1_dis, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_dis_month = "modis"
else:
    glass_source_pref_dis_month = "avhrr"


# In[ ]:





# # Selection of periods for analysis

# ## Years with extensive or concentrated leaf area index change
# 
# 5-year rolling averages are used to average out LAI changes from short-term climate fluctuations

# In[ ]:


# Create glass rolling avg files to identify years with significant LAI change
cf.calc_glass_rolling_avg_of_annual_diff(
    region=region, year_start=year_start, year_end=year_end, months_subset="all", 
    window_size=window_size, glass_source_pref="avhrr"
)
cf.calc_glass_rolling_avg_of_annual_diff(
    region=region, year_start=year_start, year_end=year_end, months_subset="all", 
    window_size=window_size, glass_source_pref="modis"
)


# In[ ]:


# Plot the glass rolling avg files to identify years with significant LAI change 
# using AVHRR data
pf.create_glass_rolling_plot(
    region=region, year_start=year_start, year_end=year_end, months_subset="all", 
    window_size=window_size, param_glass_mean="mlai", glass_source_pref="avhrr", 
    extents=None, vmin=None, vmax=None, cfv_data=None, output=True
)


# In[ ]:


# Plot the glass rolling avg files to identify years with significant LAI change 
# using MODIS data
pf.create_glass_rolling_plot(
    region=region, year_start=year_start, year_end=year_end, months_subset="all", 
    window_size=window_size, param_glass_mean="mlai", glass_source_pref="modis", 
    extents=None, vmin=None, vmax=None, cfv_data=None, output=True
)


# In[ ]:


# Check the difference in MLAI for periods with similar background
# atmospheric oscillations
pf.create_individual_comp_plot(
    calc_func=cf.calc_glass_mean_clim, region=region, 
    period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, 
    months_subset="all", arg_extra="mlai", glass_source_pref=glass_source_pref_sim, 
    var_or_dvar=None, perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, 
    vmin_periods=None, vmax_periods=None, vmin_diff=None, vmax_diff=None, 
    ax_period1=None, ax_period2=None, ax_diff=None, cfv_data=None, output=True
)


# In[ ]:


# Check the difference in MLAI for periods with dissimilar background
# atmospheric oscillations
pf.create_individual_comp_plot(
    calc_func=cf.calc_glass_mean_clim, region=region, 
    period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, 
    months_subset="all", arg_extra="mlai", glass_source_pref=glass_source_pref_dis, 
    var_or_dvar=None, perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, 
    vmin_periods=None, vmax_periods=None, vmin_diff=None, vmax_diff=None, 
    ax_period1=None, ax_period2=None, ax_diff=None, cfv_data=None, output=True
)


# ## Periods with similar / dissimilar background atmospheric oscillations
# 
# Select periods based on:
# - similar/dissimilar pattern in monthly values over each period
# - having similar/dissimilar 5-year rolling averages for relevant indices
# - for WA this is mainly the AAOI and DMI
# - for CA and SA this is mainly the AMOI, PDOI and ONI

# In[ ]:


# Create climate indices plot to help select periods with similar background 
# atmospheric oscillations
pf.create_climate_indices_plot(
    year_start=year_start, year_end=year_end, window_size=window_size, 
    period1_mid=period1_mid_sim, period2_mid=period2_mid_sim, 
    month1_mark=None, month2_mark=None, 
    cfv_data=None, output=True
)


# In[ ]:


# Create climate indices plot to help select periods with dissimilar background 
# atmospheric oscillations
pf.create_climate_indices_plot(
    year_start=year_start, year_end=year_end, window_size=window_size, 
    period1_mid=period1_mid_dis, period2_mid=period2_mid_dis, 
    month1_mark=month1_dis, month2_mark=month2_dis, 
    cfv_data=None, output=True
)


# In[ ]:





# # Create all output files

# ## Create all output data files over each season
# 
# - this outputs all possible files for the difference in results between periods
# - it does so by invoking the create_all_possible_diff_data_files function
# - this function in turn calls upon the create_all_possible_calc_data_files function
# - so intermediate data files for each period are also implicitly outputted

# In[ ]:


# Create all possible diff data files for periods with similar background 
# atmospheric oscillations (year round, wet season only, and dry season only)
cf.create_all_possible_diff_data_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim,
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all"
)
cf.create_all_possible_diff_data_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim,
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet
)
cf.create_all_possible_diff_data_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim,
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry
)


# In[ ]:


# Create all possible diff data files for periods with dissimilar background 
# atmospheric oscillations (year round, wet season only, and dry season only)
cf.create_all_possible_diff_data_files(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis,
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all"
)
cf.create_all_possible_diff_data_files(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis,
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet
)
cf.create_all_possible_diff_data_files(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis,
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry
)


# In[ ]:


# Create all possible diff data files for months with dissimilar background 
# atmospheric oscillations
cf.create_all_possible_diff_data_files(
    region=region, period1_start=month1_dis, period1_end=month1_dis,
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all"
)


# ## Create all static plot files

# In[ ]:


# Create climate indices plot with no labels
pf.create_climate_indices_plot(
    year_start=year_start, year_end=year_end, window_size=window_size, 
    period1_mid=None, period2_mid=None, month1_mark=None, month2_mark=None, 
    cfv_data=None, output=True
)


# In[ ]:


# Plot for land surface elevation
pf.create_orog_static_plot(
    param_orog="lse", region=region, extents=None, vmin=None, vmax=None, 
    ax=None, cfv_data=None, output=True
)


# In[ ]:


# Plot for slope of sub-gridscale orography
pf.create_orog_static_plot(
    param_orog="ssgo", region=region, extents=None, vmin=None, vmax=None, 
    ax=None, cfv_data=None, output=True
)


# ## Create all output comparison plot files over each season

# In[ ]:


# Create all possible comp plot files for periods with similar background 
# atmospheric oscillations (year round, wet season only, and dry season only)
pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, perc=False, 
    mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None
)
pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, perc=False, 
    mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None
)
pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, perc=False, 
    mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None
)


# In[ ]:


# Create all possible comp plot files for periods with dissimilar background 
# atmospheric oscillations (year round, wet season only, and dry season only)
pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, perc=False, 
    mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None
)
pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, perc=False, 
    mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None
)
pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, perc=False, 
    mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None
)


# In[ ]:


# Create all possible comp plot files for months with dissimilar background 
# atmospheric oscillations
pf.create_all_possible_comp_plot_files(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, perc=False, 
    mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None
)


# In[ ]:





# # Overview for similar periods

# ## Similar periods (all months)

# ### Similar periods (all months): MDP stats 

# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="wv100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="ws100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="u100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="v100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="mslp", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="t2", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Similar periods (all months): MDP values

# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Similar periods (all months): WSD

# In[ ]:


pf.plot_comp_wsd_clim(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset="all", 
    glass_source_pref=glass_source_pref_sim, 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ## Similar periods (wet season)

# ### Similar periods (wet season): MDP stats

# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="wv100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="ws100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="u100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="v100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="mslp", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="t2", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Similar periods (wet season): MDP values

# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Similar periods (wet season): WSD

# In[ ]:


pf.plot_comp_wsd_clim(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_sim, 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ## Similar periods (dry season)

# ### Similar periods (dry season): MDP stats

# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="wv100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="ws100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="u100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="v100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="mslp", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="t2", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Similar periods (dry season): MDP values

# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dwv100", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dmslp", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dt2", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="dnac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, var_or_dvar="nac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Similar periods (dry season): WSD

# In[ ]:


pf.plot_comp_wsd_clim(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_sim, 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:





# # Overview for dissimilar periods

# ## Dissimilar periods (all months)

# ### Dissimilar periods (all months): MDP stats 

# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="wv100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="ws100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="u100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="v100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="mslp", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="t2", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Dissimilar periods (all months): MDP values

# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Dissimilar periods (all months): WSD

# In[ ]:


pf.plot_comp_wsd_clim(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis, 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ## Dissimilar periods (wet season)

# ### Dissimilar periods (wet season): MDP stats

# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="wv100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="ws100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="u100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="v100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="mslp", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="t2", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Dissimilar periods (wet season): MDP values

# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Dissimilar periods (wet season): WSD

# In[ ]:


pf.plot_comp_wsd_clim(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_wet, 
    glass_source_pref=glass_source_pref_dis, 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ## Dissimilar periods (dry season)

# ### Dissimilar periods (dry season): MDP stats

# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="wv100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="ws100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="u100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="v100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="mslp", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="t2", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Dissimilar periods (dry season): MDP values

# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dwv100", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dmslp", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dt2", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="dnac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, var_or_dvar="nac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Dissimilar periods (dry season): WSD

# In[ ]:


pf.plot_comp_wsd_clim(
    region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
    period2_start=period2_start_dis, period2_end=period2_end_dis, months_subset=months_dry, 
    glass_source_pref=glass_source_pref_dis, 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:





# # Overview for dissimilar months

# ## Dissimilar months

# ### Dissimilar months: MDP stats 

# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="wv100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="ws100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="u100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="v100", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="mslp", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="t2", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_stats_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="nac", 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Dissimilar months: MDP values

# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dwv100", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dwv100", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dwv100", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dwv100", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dmslp", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dmslp", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dmslp", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dmslp", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dt2", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dt2", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dt2", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dt2", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dnac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dnac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dnac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="dnac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="nac", time="0-5",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="nac", time="6-11",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="nac", time="12-17",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:


pf.plot_comp_mdp_clim_values_given_var_or_dvar(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, var_or_dvar="nac", time="18-23",
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# ### Dissimilar months: WSD

# In[ ]:


pf.plot_comp_wsd_clim(
    region=region, period1_start=month1_dis, period1_end=month1_dis, 
    period2_start=month2_dis, period2_end=month2_dis, months_subset="all", 
    glass_source_pref=glass_source_pref_dis_month, 
    perc=False, mask_perc_quantile=pf.mask_perc_quantile_default, 
    mask_period1=None, mask_period2=None, extents=None, cfv_data=None, output=True
)


# In[ ]:





# # Selected analysis and insights

# In[ ]:





# # Scrap

# ## WA settings

# In[ ]:


# General analysis settings
region = "wa"
year_start = 1983
year_end = 2019
window_size = 5
months_wet = "mam"
months_dry = "djf"

# Periods with similar background atmospheric oscillations
period1_mid_sim = "Dec-1999"
period1_start_sim = "Jun-1997"
period1_end_sim = "May-2002"
period2_mid_sim = "Mar-2013"
period2_start_sim = "Sep-2010"
period2_end_sim = "Aug-2015"
if (datetime.strptime(period1_start_sim, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_sim = "modis"
else:
    glass_source_pref_sim = "avhrr"

# Periods with dissimilar background atmospheric oscillations
period1_mid_dis = "Sep-1990"
period1_start_dis = "Mar-1988"
period1_end_dis = "Feb-1993"
period2_mid_dis = "Jul-2016"
period2_start_dis = "Jan-2014"
period2_end_dis = "Dec-2018"
if (datetime.strptime(period1_start_dis, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_dis = "modis"
else:
    glass_source_pref_dis = "avhrr"

# Single months with very dissimilar background atmospheric oscillations
month1_dis = "Nov-1996"
month2_dis = "Sep-2015"
if (datetime.strptime(month1_dis, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_dis_month = "modis"
else:
    glass_source_pref_dis_month = "avhrr"


# ## CA settings

# In[ ]:


# General analysis settings
region = "ca"
year_start = 1983
year_end = 2019
window_size = 5
months_wet = [5, 6, 7, 8, 9, 10]
months_dry = [11, 12, 1, 2, 3, 4]

# Periods with similar background atmospheric oscillations
period1_mid_sim = "Jul-1983"
period1_start_sim = "Jan-1981"
period1_end_sim = "Dec-1985"
period2_mid_sim = "Jul-1994"
period2_start_sim = "Jan-1992"
period2_end_sim = "Dec-1996"
if (datetime.strptime(period1_start_sim, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_sim = "modis"
else:
    glass_source_pref_sim = "avhrr"
    
# Periods with dissimilar background atmospheric oscillations
period1_mid_dis = "Mar-1985"
period1_start_dis = "Sep-1982"
period1_end_dis = "Aug-1987"
period2_mid_dis = "Jul-2010"
period2_start_dis = "Jan-2008"
period2_end_dis = "Dec-2012"
if (datetime.strptime(period1_start_dis, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_dis = "modis"
else:
    glass_source_pref_dis = "avhrr"

# Single months with very dissimilar background atmospheric oscillations
month1_dis = "Nov-1996"
month2_dis = "Sep-2015"
if (datetime.strptime(month1_dis, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_dis_month = "modis"
else:
    glass_source_pref_dis_month = "avhrr"


# ## SA settings

# In[ ]:


# General analysis settings
region = "sa"
year_start = 1983
year_end = 2019
window_size = 5
months_wet = [1, 2, 3, 4, 5, 6]
months_dry = [7, 8, 9, 10, 11, 12]

# Periods with similar background atmospheric oscillations
period1_mid_sim = "Jul-2005"
period1_start_sim = "Jan-2003"
period1_end_sim = "Dec-2007"
period2_mid_sim = "Jul-2018"
period2_start_sim = "Jan-2016"
period2_end_sim = "Dec-2020"
if (datetime.strptime(period1_start_sim, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_sim = "modis"
else:
    glass_source_pref_sim = "avhrr"
    
# Periods with dissimilar background atmospheric oscillations
period1_mid_dis = "Mar-1985"
period1_start_dis = "Sep-1982"
period1_end_dis = "Aug-1987"
period2_mid_dis = "Jul-2010"
period2_start_dis = "Jan-2008"
period2_end_dis = "Dec-2012"
if (datetime.strptime(period1_start_dis, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_dis = "modis"
else:
    glass_source_pref_dis = "avhrr"

# Single months with very dissimilar background atmospheric oscillations
month1_dis = "Nov-1996"
month2_dis = "Sep-2015"
if (datetime.strptime(month1_dis, "%b-%Y") >= 
    datetime.strptime(cf.modis_earliest, "%b-%Y")):
    glass_source_pref_dis_month = "modis"
else:
    glass_source_pref_dis_month = "avhrr"


# In[ ]:




