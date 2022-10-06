import importlib
cf = importlib.import_module("test_profiler_calc_funcs_v1n")
# cf.priority = "memory"
cf.calc_glass_rolling_avg_of_annual_diff("wa", 1983, 2019, "all", 5, "modis")
# cf.calc_glass_rolling_avg_of_annual_diff("global", 1983, 2019, "all", 5, "avhrr")
