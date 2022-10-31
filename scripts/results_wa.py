#!/usr/bin/env python
# coding: utf-8

# # Setup
# 
# Note: the following code requires up to 70 GB of RAM over 18 hours if using 6 CPUs
# 
# It can also be run with less RAM but this will require manual restarting of code everytime RAM limit is reached

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
region = "wa"
year_start = 1983
year_end = 2019
window_size = 5
months_wet = "jja"
months_dry = "djf"
hours_light = [8, 9, 10, 11, 12, 13, 14, 15]
hours_night = [20, 21, 22, 23, 0, 1, 2, 3]
perc = False
mask_perc_quantile = pf.mask_perc_quantile_default
mask_period1 = None
mask_period2 = None
extents = None
cfv_data = None

# Entire periods covered by AVHRR or MODIS data and for which there is FAPAR data
period_start_avhrr = "Jan-1982"
period_end_avhrr = "Dec-2018"
period_start_modis = "Jan-2001"
period_end_modis = "Dec-2020"

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

# # Periods with dissimilar background atmospheric oscillations
# period1_mid_dis = "Sep-1990"
# period1_start_dis = "Mar-1988"
# period1_end_dis = "Feb-1993"
# period2_mid_dis = "Jul-2016"
# period2_start_dis = "Jan-2014"
# period2_end_dis = "Dec-2018"
# if (datetime.strptime(period1_start_dis, "%b-%Y") >= 
#     datetime.strptime(cf.modis_earliest, "%b-%Y")):
#     glass_source_pref_dis = "modis"
# else:
#     glass_source_pref_dis = "avhrr"

# # Single months with very dissimilar background atmospheric oscillations
# month1_dis = "Nov-1996"
# month2_dis = "Sep-2015"
# if (datetime.strptime(month1_dis, "%b-%Y") >= 
#     datetime.strptime(cf.modis_earliest, "%b-%Y")):
#     glass_source_pref_dis_month = "modis"
# else:
#     glass_source_pref_dis_month = "avhrr"


# In[ ]:


# # General analysis settings
# region = "ca"
# year_start = 1983
# year_end = 2019
# window_size = 5
# months_wet = [5, 6, 7, 8, 9, 10]
# months_dry = [11, 12, 1, 2, 3, 4]
# hours_light = [8, 9, 10, 11, 12, 13, 14, 15]
# hours_night = [20, 21, 22, 23, 0, 1, 2, 3]
# perc = False
# mask_perc_quantile = pf.mask_perc_quantile_default
# mask_period1 = None
# mask_period2 = None
# extents = None
# cfv_data = None

# # Entire periods covered by AVHRR or MODIS data and for which there is FAPAR data
# period_start_avhrr = "Jan-1982"
# period_end_avhrr = "Dec-2018"
# period_start_modis = "Jan-2001"
# period_end_modis = "Dec-2020"

# # Periods with similar background atmospheric oscillations
# period1_mid_sim = "Jul-2005"
# period1_start_sim = "Jan-2003"
# period1_end_sim = "Dec-2007"
# period2_mid_sim = "Jul-2018"
# period2_start_sim = "Jan-2016"
# period2_end_sim = "Dec-2020"
# if (datetime.strptime(period1_start_sim, "%b-%Y") >= 
#     datetime.strptime(cf.modis_earliest, "%b-%Y")):
#     glass_source_pref_sim = "modis"
# else:
#     glass_source_pref_sim = "avhrr"
    
# # # Periods with dissimilar background atmospheric oscillations
# # period1_mid_dis = "Mar-1985"
# # period1_start_dis = "Sep-1982"
# # period1_end_dis = "Aug-1987"
# # period2_mid_dis = "Jul-2010"
# # period2_start_dis = "Jan-2008"
# # period2_end_dis = "Dec-2012"
# # if (datetime.strptime(period1_start_dis, "%b-%Y") >= 
# #     datetime.strptime(cf.modis_earliest, "%b-%Y")):
# #     glass_source_pref_dis = "modis"
# # else:
# #     glass_source_pref_dis = "avhrr"

# # # Single months with very dissimilar background atmospheric oscillations
# # month1_dis = "Nov-1996"
# # month2_dis = "Sep-2015"
# # if (datetime.strptime(month1_dis, "%b-%Y") >= 
# #     datetime.strptime(cf.modis_earliest, "%b-%Y")):
# #     glass_source_pref_dis_month = "modis"
# # else:
# #     glass_source_pref_dis_month = "avhrr"


# In[ ]:


# # General analysis settings
# region = "sa"
# year_start = 1983
# year_end = 2019
# window_size = 5
# months_wet = [1, 2, 3, 4, 5, 6]
# months_dry = [7, 8, 9, 10, 11, 12]
# hours_light = [8, 9, 10, 11, 12, 13, 14, 15]
# hours_night = [20, 21, 22, 23, 0, 1, 2, 3]
# perc = False
# mask_perc_quantile = pf.mask_perc_quantile_default
# mask_period1 = None
# mask_period2 = None
# extents = None
# cfv_data = None

# # Entire periods covered by AVHRR or MODIS data and for which there is FAPAR data
# period_start_avhrr = "Jan-1982"
# period_end_avhrr = "Dec-2018"
# period_start_modis = "Jan-2001"
# period_end_modis = "Dec-2020"

# # Periods with similar background atmospheric oscillations
# period1_mid_sim = "Jul-2005"
# period1_start_sim = "Jan-2003"
# period1_end_sim = "Dec-2007"
# period2_mid_sim = "Jul-2018"
# period2_start_sim = "Jan-2016"
# period2_end_sim = "Dec-2020"
# if (datetime.strptime(period1_start_sim, "%b-%Y") >= 
#     datetime.strptime(cf.modis_earliest, "%b-%Y")):
#     glass_source_pref_sim = "modis"
# else:
#     glass_source_pref_sim = "avhrr"
    
# # # Periods with dissimilar background atmospheric oscillations
# # period1_mid_dis = "Mar-1985"
# # period1_start_dis = "Sep-1982"
# # period1_end_dis = "Aug-1987"
# # period2_mid_dis = "Jul-2010"
# # period2_start_dis = "Jan-2008"
# # period2_end_dis = "Dec-2012"
# # if (datetime.strptime(period1_start_dis, "%b-%Y") >= 
# #     datetime.strptime(cf.modis_earliest, "%b-%Y")):
# #     glass_source_pref_dis = "modis"
# # else:
# #     glass_source_pref_dis = "avhrr"

# # # Single months with very dissimilar background atmospheric oscillations
# # month1_dis = "Nov-1996"
# # month2_dis = "Sep-2015"
# # if (datetime.strptime(month1_dis, "%b-%Y") >= 
# #     datetime.strptime(cf.modis_earliest, "%b-%Y")):
# #     glass_source_pref_dis_month = "modis"
# # else:
# #     glass_source_pref_dis_month = "avhrr"


# In[ ]:





# # Selection of periods for analysis

# ## Years with extensive or concentrated leaf area index change
# 
# 5-year rolling averages are used to average out LAI changes from short-term climate fluctuations

# In[ ]:


# Plot the glass rolling avg files to identify years with significant LAI change 
# using AVHRR data
pf.create_glass_rolling_plot(
    region=region, year_start=year_start, year_end=year_end, period_months="all", 
    window_size=window_size, param_glass_mean="mlai", glass_source_pref="avhrr", 
    extents=extents, vmin=None, vmax=None, cfv_data=cfv_data, output=True
)


# In[ ]:


# Plot the glass rolling avg files to identify years with significant LAI change 
# using MODIS data
pf.create_glass_rolling_plot(
    region=region, year_start=year_start, year_end=year_end, period_months="all", 
    window_size=window_size, param_glass_mean="mlai", glass_source_pref="modis", 
    extents=extents, vmin=None, vmax=None, cfv_data=cfv_data, output=True
)


# In[ ]:


# Check the difference in MLAI for periods with similar background
# atmospheric oscillations
pf.create_individual_comp_plot(
    calc_func=cf.calc_glass_mean_clim, region=region, 
    period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, 
    period1_months="all", period2_months="all", period1_hours=None, period2_hours=None,
    arg_extra="mlai", glass_source_pref=glass_source_pref_sim, 
    var_or_dvar=None, perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, extents=extents, 
    vmin_periods=None, vmax_periods=None, vmin_diff=None, vmax_diff=None, 
    ax_period1=None, ax_period2=None, ax_diff=None, cfv_data=cfv_data, output=True
)


# In[ ]:


# # Check the difference in MLAI for periods with dissimilar background
# # atmospheric oscillations
# pf.create_individual_comp_plot(
#     calc_func=cf.calc_glass_mean_clim, region=region, 
#     period1_start=period1_start_dis, period1_end=period1_end_dis, 
#     period2_start=period2_start_dis, period2_end=period2_end_dis, 
#     period1_months="all", period2_months="all", period1_hours=None, period2_hours=None,
#     arg_extra="mlai", glass_source_pref=glass_source_pref_dis, 
#     var_or_dvar=None, perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, extents=extents, 
#     vmin_periods=None, vmax_periods=None, vmin_diff=None, vmax_diff=None, 
#     ax_period1=None, ax_period2=None, ax_diff=None, cfv_data=cfv_data, output=True
# )


# In[ ]:





# ## Periods with similar / dissimilar background atmospheric oscillations
# 
# Criteria for selection of periods:
# 
# - 5-year rolling averages for relevant climate indices were similar / dissimilar
# - Change in leaf area index between the periods was extensive
# - Periods cover a similar / dissimilar amount of time spent in La Nina / El Nino (where relevant) and Negative / Positive Indian Ocean Dipole (where relevant)
# - Monthly values for relevant climate indices over each period display a similar / dissimilar time evolution pattern
# - For WA, relevant indices are AAOI and DMI
# - For CA and SA, relevant indices are AMOI, PDOI and ONI

# In[ ]:


# Create climate indices plot to help select periods with similar background 
# atmospheric oscillations
pf.create_climate_indices_plot(
    year_start=year_start, year_end=year_end, window_size=window_size, 
    period1_mid=period1_mid_sim, period2_mid=period2_mid_sim, 
    month1_mark=None, month2_mark=None, 
    cfv_data=cfv_data, output=True
)


# In[ ]:


# # Create climate indices plot to help select periods with dissimilar background 
# # atmospheric oscillations
# pf.create_climate_indices_plot(
#     year_start=year_start, year_end=year_end, window_size=window_size, 
#     period1_mid=period1_mid_dis, period2_mid=period2_mid_dis, 
#     # month1_mark=month1_dis, month2_mark=month2_dis, 
#     cfv_data=cfv_data, output=True
# )


# In[ ]:





# # Create all output plot files

# ## Static plots

# In[ ]:


# Create climate indices plot with no labels
pf.create_climate_indices_plot(
    year_start=year_start, year_end=year_end, window_size=window_size, 
    period1_mid=None, period2_mid=None, month1_mark=None, month2_mark=None, 
    cfv_data=cfv_data, output=True
)


# In[ ]:


# Combined plot for orography of each region
pf.create_comb_orog_static_plot(cfv_data=cfv_data, output=True)


# In[ ]:


# Plot for land surface elevation
pf.create_orog_static_plot(
    param_orog="lse", region=region, extents=extents, vmin=None, vmax=None, 
    ax=None, cfv_data=cfv_data, output=True
)


# In[ ]:


# Plot for slope of sub-gridscale orography
pf.create_orog_static_plot(
    param_orog="ssgo", region=region, extents=extents, vmin=None, vmax=None, 
    ax=None, cfv_data=cfv_data, output=True
)


# In[ ]:





# ## Entire modis period

# ### Entire modis period, wet vs dry months

# #### Entire modis period, wet vs dry months, all hours

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period_start_modis, period1_end=period_end_modis, 
    period2_start=period_start_modis, period2_end=period_end_modis, 
    period1_months=months_wet, period2_months=months_dry, 
    period1_hours="all", period2_hours="all",
    glass_source_pref="modis", 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# #### Entire modis period, wet vs dry months, day hours

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period_start_modis, period1_end=period_end_modis, 
    period2_start=period_start_modis, period2_end=period_end_modis, 
    period1_months=months_wet, period2_months=months_dry, 
    period1_hours=hours_light, period2_hours=hours_light,
    glass_source_pref="modis", 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# #### Entire modis period, wet vs dry months, night hours

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period_start_modis, period1_end=period_end_modis, 
    period2_start=period_start_modis, period2_end=period_end_modis, 
    period1_months=months_wet, period2_months=months_dry, 
    period1_hours=hours_night, period2_hours=hours_night,
    glass_source_pref="modis", 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# ### Entire modis period, day vs night hours

# #### Entire modis period, day vs night hours, all months

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period_start_modis, period1_end=period_end_modis, 
    period2_start=period_start_modis, period2_end=period_end_modis, 
    period1_months="all", period2_months="all", 
    period1_hours=hours_light, period2_hours=hours_night,
    glass_source_pref="modis", 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# #### Entire modis period, day vs night hours, wet months

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period_start_modis, period1_end=period_end_modis, 
    period2_start=period_start_modis, period2_end=period_end_modis, 
    period1_months=months_wet, period2_months=months_wet, 
    period1_hours=hours_light, period2_hours=hours_night,
    glass_source_pref="modis", 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# #### Entire modis period, day vs night hours, dry months

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period_start_modis, period1_end=period_end_modis, 
    period2_start=period_start_modis, period2_end=period_end_modis, 
    period1_months=months_dry, period2_months=months_dry, 
    period1_hours=hours_light, period2_hours=hours_night,
    glass_source_pref="modis", 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# In[ ]:





# ## Similar periods

# ### Similar periods, wet vs dry months

# #### Entire modis period, wet vs dry months, all hours

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, 
    period1_months=months_wet, period2_months=months_dry, 
    period1_hours="all", period2_hours="all",
    glass_source_pref=glass_source_pref_sim, 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# #### Similar periods, wet vs dry months, day hours

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, 
    period1_months=months_wet, period2_months=months_dry, 
    period1_hours=hours_light, period2_hours=hours_light,
    glass_source_pref=glass_source_pref_sim, 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# #### Similar periods, wet vs dry months, night hours

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, 
    period1_months=months_wet, period2_months=months_dry, 
    period1_hours=hours_night, period2_hours=hours_night,
    glass_source_pref=glass_source_pref_sim, 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# ### Similar periods, day vs night hours

# #### Similar periods, day vs night hours, all months

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, 
    period1_months="all", period2_months="all", 
    period1_hours=hours_light, period2_hours=hours_night,
    glass_source_pref=glass_source_pref_sim, 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# #### Similar periods, day vs night hours, wet months

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, 
    period1_months=months_wet, period2_months=months_wet, 
    period1_hours=hours_light, period2_hours=hours_night,
    glass_source_pref=glass_source_pref_sim, 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# #### Similar periods, day vs night hours, dry months

# In[ ]:


pf.create_all_possible_comp_plot_files(
    region=region, period1_start=period1_start_sim, period1_end=period1_end_sim, 
    period2_start=period2_start_sim, period2_end=period2_end_sim, 
    period1_months=months_dry, period2_months=months_dry, 
    period1_hours=hours_light, period2_hours=hours_night,
    glass_source_pref=glass_source_pref_sim, 
    perc=perc, mask_perc_quantile=mask_perc_quantile, 
    mask_period1=mask_period1, mask_period2=mask_period2, 
    extents=extents, cfv_data=cfv_data
)


# In[ ]:





# ## Entire avhrr period

# ### Entire avhrr period, wet vs dry months

# #### Entire avhrr period, wet vs dry months, all hours

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period_start_avhrr, period1_end=period_end_avhrr, 
#     period2_start=period_start_avhrr, period2_end=period_end_avhrr, 
#     period1_months=months_wet, period2_months=months_dry, 
#     period1_hours="all", period2_hours="all",
#     glass_source_pref="avhrr", 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# #### Entire avhrr period, wet vs dry months, day hours

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period_start_avhrr, period1_end=period_end_avhrr, 
#     period2_start=period_start_avhrr, period2_end=period_end_avhrr, 
#     period1_months=months_wet, period2_months=months_dry, 
#     period1_hours=hours_light, period2_hours=hours_light,
#     glass_source_pref="avhrr", 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# #### Entire avhrr period, wet vs dry months, night hours

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period_start_avhrr, period1_end=period_end_avhrr, 
#     period2_start=period_start_avhrr, period2_end=period_end_avhrr, 
#     period1_months=months_wet, period2_months=months_dry, 
#     period1_hours=hours_night, period2_hours=hours_night,
#     glass_source_pref="avhrr", 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# ### Entire avhrr period, day vs night hours

# #### Entire avhrr period, day vs night hours, all months

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period_start_avhrr, period1_end=period_end_avhrr, 
#     period2_start=period_start_avhrr, period2_end=period_end_avhrr, 
#     period1_months="all", period2_months="all", 
#     period1_hours=hours_light, period2_hours=hours_night,
#     glass_source_pref="avhrr", 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# #### Entire avhrr period, day vs night hours, wet months

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period_start_avhrr, period1_end=period_end_avhrr, 
#     period2_start=period_start_avhrr, period2_end=period_end_avhrr, 
#     period1_months=months_wet, period2_months=months_wet, 
#     period1_hours=hours_light, period2_hours=hours_night,
#     glass_source_pref="avhrr", 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# #### Entire avhrr period, day vs night hours, dry months

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period_start_avhrr, period1_end=period_end_avhrr, 
#     period2_start=period_start_avhrr, period2_end=period_end_avhrr, 
#     period1_months=months_dry, period2_months=months_dry, 
#     period1_hours=hours_light, period2_hours=hours_night,
#     glass_source_pref="avhrr", 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# In[ ]:





# ## Dissimilar periods

# ### Dissimilar periods, wet vs dry months

# #### Dissimilar periods, wet vs dry months, all hours

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
#     period2_start=period2_start_dis, period2_end=period2_end_dis, 
#     period1_months=months_wet, period2_months=months_dry, 
#     period1_hours="all", period2_hours="all",
#     glass_source_pref=glass_source_pref_dis, 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# #### Dissimilar periods, wet vs dry months, day hours

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
#     period2_start=period2_start_dis, period2_end=period2_end_dis, 
#     period1_months=months_wet, period2_months=months_dry, 
#     period1_hours=hours_light, period2_hours=hours_light,
#     glass_source_pref=glass_source_pref_dis, 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# #### Dissimilar periods, wet vs dry months, night hours

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
#     period2_start=period2_start_dis, period2_end=period2_end_dis, 
#     period1_months=months_wet, period2_months=months_dry, 
#     period1_hours=hours_night, period2_hours=hours_night,
#     glass_source_pref=glass_source_pref_dis, 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# ### Dissimilar periods, day vs night hours

# #### Dissimilar periods, day vs night hours, all months

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
#     period2_start=period2_start_dis, period2_end=period2_end_dis, 
#     period1_months="all", period2_months="all", 
#     period1_hours=hours_light, period2_hours=hours_night,
#     glass_source_pref=glass_source_pref_dis, 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# #### Dissimilar periods, day vs night hours, wet months

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
#     period2_start=period2_start_dis, period2_end=period2_end_dis, 
#     period1_months=months_wet, period2_months=months_wet, 
#     period1_hours=hours_light, period2_hours=hours_night,
#     glass_source_pref=glass_source_pref_dis, 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )


# #### Dissimilar periods, day vs night hours, dry months

# In[ ]:


# pf.create_all_possible_comp_plot_files(
#     region=region, period1_start=period1_start_dis, period1_end=period1_end_dis, 
#     period2_start=period2_start_dis, period2_end=period2_end_dis, 
#     period1_months=months_dry, period2_months=months_dry, 
#     period1_hours=hours_light, period2_hours=hours_night,
#     glass_source_pref=glass_source_pref_dis, 
#     perc=perc, mask_perc_quantile=mask_perc_quantile, 
#     mask_period1=mask_period1, mask_period2=mask_period2, 
#     extents=extents, cfv_data=cfv_data
# )
