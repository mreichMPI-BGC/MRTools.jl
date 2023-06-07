using RCall
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


function proj(grid, outcrs="+proj=robin +datum=WGS84")# +units=m +no_defs")
    #grid = permutedims(grid, (2,1))
    jras = R"""
        require(terra)
        rast($(permutedims(grid, (2,1))), crs = "+proj=longlat +datum=WGS84 +no_defs", extent=ext(rast())) %>% project($outcrs) %>% as.matrix(wide=T) 
    """ |> rcopy
    return permutedims(jras, (2,1))
end


 dem05() = R"""load(file="D:\\markusr\\daten\\geodata\\globalDEM0.5deg.rda"); globalDEM0.5deg %>% as.matrix()""" |> 
         rcopy  .|> Float32 |> rotr90



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