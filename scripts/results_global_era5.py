#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# The following code requires up to XXX GB of RAM over YYY hours if using ZZZ CPUs
# Up to AAA GB of storage will be required for the outputs


# In[ ]:


# Import libraries
import importlib
from glob import glob


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


# In[ ]:


# Increase resolution of output figures since global images are very large
pf.plt.rcParams['savefig.dpi'] = 400


# In[ ]:


# Create all possible diff data files for period of LAI loss (using AVHRR data)
# Do this for year round only (since seasons are not well defined on global scale)
cf.create_all_possible_diff_data_files("global", "Jan-1981", "Dec-1985", "Jan-1992", "Dec-1996", "all")


# In[ ]:


# Create all possible diff data files for period of LAI gain (using MODIS data)
# Do this for year round only (since seasons are not well defined on global scale)
cf.create_all_possible_diff_data_files("global", "Jan-2007", "Dec-2011", "Jan-2017", "Dec-2021", "all")


# In[ ]:


# Create the relevant ERA5 analysis plots

