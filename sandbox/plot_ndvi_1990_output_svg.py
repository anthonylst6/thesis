import metview as mv

# Importing : /Thesis/data_raw/noaacdrndviv5-global/AVHRR-Land_v005_AVH13C1_NOAA-11_19900710_c20170615012302.nc

avhrr_2d_land_v005_avh13c1_noaa_2d_11_19900710_c20170615012302_2e_nc = mv.read("/home/anthony/metview//Thesis/data_raw/noaacdrndviv5-global/AVHRR-Land_v005_AVH13C1_NOAA-11_19900710_c20170615012302.nc")

# Importing : /Thesis/data_raw/noaacdrndviv5-global/Coastlines

coastlines = mv.mcoast(
    map_coastline_thickness  = 2,
    map_boundaries           = "on",
    map_cities               = "on",
    map_boundaries_colour    = "black",
    map_boundaries_thickness = 2
    )

# Importing : /Thesis/data_raw/noaacdrndviv5-global/Geographical View

geographical_view = mv.geoview(
    map_area_definition = "corners",
    area                = [7,-91,17,-81],
    coastlines          = coastlines
    )

# Importing : /Thesis/data_raw/noaacdrndviv5-global/NetCDF Visualiser

netcdf_visualiser = mv.netcdf_visualiser(
    netcdf_plot_type          = "geo_matrix",
    netcdf_latitude_variable  = "latitude",
    netcdf_longitude_variable = "longitude",
    netcdf_value_variable     = "NDVI",
    netcdf_data               = avhrr_2d_land_v005_avh13c1_noaa_2d_11_19900710_c20170615012302_2e_nc
    )

# Importing : /Thesis/data_raw/noaacdrndviv5-global/Contouring

contouring = mv.mcont(
    legend                         = "on",
    contour_level_selection_type   = "level_list",
    contour_level_list             = [-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1],
    contour_shade                  = "on",
    contour_shade_technique        = "grid_shading",
    contour_shade_max_level_colour = "evergreen",
    contour_shade_min_level_colour = "white"
    )

mv.setoutput(mv.svg_output(
    output_name="/home/anthony/metview/Thesis/sandbox/sandbox_output/ndvi1990",
    output_name_first_page_number="off"),
    output_svg_meta="test metadata",
    output_svg_desc="test desc",
    output_title="test title",
    output_debug="on")

mv.plot(geographical_view, netcdf_visualiser, contouring)
