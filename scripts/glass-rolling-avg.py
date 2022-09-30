import importlib
cf = importlib.import_module("calc_funcs_v1n")
cf.calc_glass_rolling_avg_of_annual_diff("sa", 1983, 2019, "all", 5, "avhrr")
cf.calc_glass_rolling_avg_of_annual_diff("sa", 1983, 2019, "all", 5, "modis")
cf.calc_glass_rolling_avg_of_annual_diff("ca", 1983, 2019, "all", 5, "avhrr")
cf.calc_glass_rolling_avg_of_annual_diff("ca", 1983, 2019, "all", 5, "modis")
cf.calc_glass_rolling_avg_of_annual_diff("wa", 1983, 2019, "all", 5, "avhrr")
cf.calc_glass_rolling_avg_of_annual_diff("wa", 1983, 2019, "all", 5, "modis")
