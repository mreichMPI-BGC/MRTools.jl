import Rasters
import Downloads, Shapefile, Extents
#@rlibrary terra
#terra = rimport("terra")
#using .terra: rast, ext, project

export countryRaster

hrev(x::Array) = reverse(x; dims=2)

function countryRaster(; res=0.5)
    shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
    shapefile_name = "country_borders.shp"
    isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)
    countries = Shapefile.Handle(shapefile_name).shapes
    world = Rasters.rasterize(countries; to = Extents.Extent(X = (-180, 180), Y = (-90, 90)), res=res, missingval=0, fill=1, shape=:line)
    #return world.data.data
end






#heatmap(proj(rand(180,360)), interpolate=false)
 
# reso = 0.5
# shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
# shapefile_name = "country_borders.shp"
# isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)
# countries = Shapefile.Handle(shapefile_name).shapes
# world = Rasters.rasterize(countries; crs=EPSG(4326), to = Extents.Extent(X = (-180, 180), Y = (-90, 90)), res=reso, missingval=0, fill=1, shape=:line)
# crs = ProjString("+proj=robin +datum=WGS84")
# wsg84_world = Rasters.resample(world, reso; crs)
# Makie.heatmap(wsg84_world)


# rast = rand((X(Projected(1:100;crs=EPSG(4326))), Y(Projected(1:100;crs=EPSG(4326)))))