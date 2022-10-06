import importlib
pf = importlib.import_module("plot_funcs_v1d")
pf.plt.rcParams['savefig.dpi'] = 600
pf.create_glass_rolling_plot("global", 1983, 2019, "all", 5, "mlai", glass_source_pref="modis", output=True)
# pf.create_glass_rolling_plot("global", 1983, 2019, "all", 5, "mlai", glass_source_pref="avhrr", output=True)
