using RCall

function loess2Dsurf_R(x,y,z)
    df=DataFrame(;x,y,z)
    R"""
    d <- $df
    model= loess(z ~ x + y, d, span=0.2)
    xs = seq(min(d$x), max(d$x), length.out=30)
    ys = seq(min(d$y), max(d$y), length.out=30)
    gr = expand.grid(xs, ys); names(gr)  <- c("x", "y")
    res = list(xs, ys, predict(model, gr))
    """ |> rcopy

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