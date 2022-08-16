import metview as mv

# Importing : /sandbox/sandbox_output/Coastlines

coastlines = mv.mcoast(
    map_coastline_land_shade        = "on",
    map_coastline_land_shade_colour = "cream",
    map_boundaries                  = "on"
    )

# Importing : /sandbox/sandbox_output/Geographical View

geographical_view = mv.geoview(
    map_area_definition = "corners",
    area                = [7,-91,17,-81],
    coastlines          = coastlines
    )

# Importing : /sandbox/sandbox_output/wsmax.nc

wsmax_2e_nc = mv.read("/home/anthony/metview//sandbox/sandbox_output/wsmax.nc")

# Importing : /sandbox/sandbox_output/NetCDF Visualiser

netcdf_visualiser = mv.netcdf_visualiser(
    netcdf_plot_type            = "geo_matrix_vectors",
    netcdf_latitude_variable    = "latitude",
    netcdf_longitude_variable   = "longitude",
    netcdf_x_component_variable = "umax",
    netcdf_y_component_variable = "vmax",
    netcdf_data                 = wsmax_2e_nc
    )

mv.plot(geographical_view, netcdf_visualiser)