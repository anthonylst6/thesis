#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# The following code requires up to 100 GB of RAM over 2 hours if using 8 CPUs
# Up to 45 GB of storage will be required for the outputs


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


# Create glass rolling avg files to identify years with significant LAI change
cf.calc_glass_rolling_avg_of_annual_diff("global", 1983, 2019, "all", 5, "avhrr")
cf.calc_glass_rolling_avg_of_annual_diff("global", 1983, 2019, "all", 5, "modis")


# In[ ]:


# Plot the glass rolling avg files to identify years with significant LAI change
pf.create_glass_rolling_plot("global", 1983, 2019, "all", 5, "mlai", glass_source_pref="avhrr", output=True)
pf.create_glass_rolling_plot("global", 1983, 2019, "all", 5, "mlai", glass_source_pref="modis", output=True)


# In[ ]:




