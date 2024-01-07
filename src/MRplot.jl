#using CairoMakie
using GLMakie, AlgebraOfGraphics
using ColorSchemes, Colors

# colors are sampled in [0,1], so define a function that linearly transforms x ∈ [lo, hi] to [0, 1].
to_unitrange(x, lo, hi) = (x - lo) / (hi - lo)

# cols is for simplicity here assumed to be a `ColorScheme` (e.g. `Makie.ColorSchemes.RdBu`)
# Example in combi with sharpen_colorscheme: chain ColorSchemes.BrBg sharpen_colorscheme() create_cmap(-5, 10)
function create_cmap(cols, lo, hi; categorical=false, stretch=0.0, mid=0.0, rev=false, kwargs...)
    if rev cols=reverse(cols) end
    absmax = maximum(abs, (lo, hi));
    lo_m, hi_m = to_unitrange.([lo, hi], -absmax + mid, absmax + mid);  # map lo, hi ∈ [-absmax, absmax] onto [0,1] to sample their corresponding colors
    colsvals = range(0, 1; length=length(cols))  # values on [0,1] where each color in cols is defined
    filter_colsvals = filter(∈(lo_m..hi_m), Base.unique([lo_m; colsvals; hi_m]))  # filter colsvals, keep only values in [lo_m, hi_m] + the endpoints lo_m and hi_m.
    mini, maxi = extrema(filter_colsvals)
    #println(filter_colsvals)
    if mini > 1e-3
      lamb=stretch * (0.5/(0.5 - mini)) + (1-stretch)
      filter_colsvalsNew = map(x -> x<0.5 ? 0.5 - lamb * (0.5 - x) : x , filter_colsvals) 
    elseif 1.0-maxi > 1e-3
      lamb=stretch * (0.5/(maxi - 0.5)) + (1-stretch)
      filter_colsvalsNew = map(x -> x>0.5 ? 0.5 + lamb * (x - 0.5) : x , filter_colsvals) 
    else 
      filter_colsvalsNew = filter_colsvals
    end
    #println(filter_colsvals)
    mini, maxi = extrema(filter_colsvals)
    newcols = get(cols, filter_colsvalsNew);  # select colors in filtered range; interpolate new low and hi colors.
    new_colsvals = to_unitrange.(filter_colsvals, lo_m, hi_m)  # values on [0,1] where the new colors are defined
    cmap = cgrad(newcols, new_colsvals; categorical)  # for continous cmap, set categorical=false
end

# Example in combi with create_cmap: @chain ColorSchemes.BrBg sharpen_colorscheme() create_cmap(-5, 10)
function sharpen_colorscheme(cs::ColorScheme; rev=false, fraction=0.25, repl=colorant"grey80", relwidth=0.01)
  isa(repl, String) && (repl = parse(Colorant, repl))
  cols = cs.colors
  #if reverse reverse!(cols); end
  numCols = trunc(Int, 1 / relwidth)
  numCols += iseven(numCols)
  c = cgrad(rev ? cols : reverse(cols), range(0, 1, numCols), categorical=true)

  cutMiddleAndReplace(c.colors[:], fraction, repl=repl)

end


function scale_z3d_ax(ax=current_axis(); zrel=0.3, yrel=1.0)
  plt=ax.scene.plots[1] ## Returns e.g. Surface{Tuple{StepRangeLen{Float32, Float64, Float64, Int64}, StepRangeLen{Float32, Float64, Float64, Int64}, Matrix{Float32}}}
  if !(typeof(plt) <: Combined{Makie.axis3d}) 
    ranges=map(i-> (filter(isfinite,plt[i][]) |> extrema |> collect |> diff)[1], 1:3)
  else
    ranges=map(i-> (filter(isfinite,plt[1][][i]) |> collect |> diff)[1], 1:3)
    #error("Not impemented!")
  end

  #println(ranges)
  relscale=ranges[1] ./ ranges
  relscale[2] *= yrel
  relscale[3] *= zrel
  #relscale = relscale ./ maximum(relscale)
  GLMakie.scale!(ax.scene,  0.5relscale...)
end

HDFigure() = Figure(; resolution=(1920, 1080))

one2one(; color=:grey, linestyle=:dash, kw...) =  mapping([[0]], [[1]]) * visual(ABLines; color, linestyle, kw...)
zeroline = mapping( [[0]]) * visual(HLines, color=:red, linestyle=:dash)

function evalplot(df, model::Symbol, obs::Symbol; kw...)
  data(df) * mapping(model, obs; kw...) * (visual(Scatter; kw... ) + linear()) + one2one()
end

function evalplot_AoG(mod::AbstractArray{T}, obs::AbstractArray{T}) where T<:Number
  df = DataFrame(mod=convert(Vector, mod |> vec) , obs=convert(Vector, obs |> vec))
  @transform!(df, :res = :obs - :mod)
  len = nrow(df)
  alpha = 1.0 / (log10(len) - 1)^3

  one2one =  mapping([[0]], [[1]]) * visual(ABLines, color=:red, linestyle=:dash)
  zeroline = mapping( [[0]]) * visual(HLines, color=:red, linestyle=:dash)
  p = data(df) * mapping(:mod, :obs)
  p*= visual(; color=(:black, alpha)) + linear()
  p+= one2one
  #AlgebraOfGraphics.draw(p; axis=(width=800, height=600))
  resPlot= data(df) * mapping(:mod,  :res) * (visual(; color=(:black, alpha)) + smooth(degree=1, span=0.75)) +   zeroline

  resolution = (800, 600)
  fig = Figure(; resolution)
  
  draw!(fig[1,2], resPlot)
  draw!(fig[1,1], p)
  fig


  #AlgebraOfGraphics.draw(resPlot; axis=(width=800, height=600))
end

function vari_facets(df, x, facets; verbose=false, reso = (1920,1080), kw...)
  function assign_facet(v)
      fidx = findall(facets) do fa
          v in String.(fa)
      end
      fidx = isnothing(fidx) ? missing : fidx
      return fidx
  end
  vars2plot = [x; vcat(facets...)]
  dflong = DataFrames.stack(@select(df, {vars2plot}) , Not(x))
  #dflong = DataFrames.stack(df , Not(x))
  
  @transform! dflong :facet =  assign_facet(:variable)
  verbose && print(dflong)

  plts = map(enumerate(facets)) do (i, fa)
      p = data(@subset(dflong |> dropmissing, i in :facet)) *
       #( mapping(x, :value, marker=:variable, color=:variable) * (visual(Scatter)) +
        mapping(x, :value, color=:variable) * visual(Lines; kw...)
       #)
   end
  #return plts

  #https://discourse.julialang.org/t/plotting-hierachically-organized-variables/104978
  # This comprehension replaces all of your function till `fig=...`
  # plts = [data(df) * mapping(:time, grp .=> :value, color = dims(1) => renamer(grp)) * visual(Lines)
  # Works for good data!! need to check if it works for corner cases, too (e.g. missings)
#for grp in facets]
   fig=Figure(resolution=reso)
  
   for i in 1:length(plts)
      g=draw!(fig[i,1], plts[i], axis=(; xlabelvisible=i == length(plts)))
      legend!(fig[i,2],g)
   end
  fig
end

"""
    taylor_plot(sdmod_rel::Array{Number}, correl::Array{Number}; bias=0.0, sdObs=1.0, plotBias=false, plotVarianceErr=false)

Create an (extended) Taylorplot.

Given the model's standard deviation relative to the observations `sdmod_rel`, the correlation with the
observations `correl`, and optionally the `bias`, an extended Taylor plot is created. Extended means, that apart
from the classical one (zero centered residuals, no bias) a point including the bias (indicated with an arrows
towards this point) and also a point under the scenario "= no variance error" (again with arrow) is plotted, 
when `plotBias` and `plotVarianceErr` are `true`. `sdObs` is the standard deviation of the observations (it just scales the
axes essentially.)

Please note that `sdmod_rel` and `correl` need to be Arrays (otherwise `arrows!` does not work), so for only one point
wrap it in [].

# Examples
```julia-repl
julia> taylor_plot([0.8,1.2], [0.8, 0.5], bias=0.2, plotBias=true, plotVarianceErr=true)

```
"""
function taylor_plot(sdmod_rel, correl; 
    bias=0.0, sdObs=1.0, plotBias=false, plotVarianceErr=false, labels=string.(1:length(correl)), markersize=20,  kw...)
    
    sdmax = max(maximum(sdmod_rel), 1.0)
    
    ## Function to transform coordinates
    correl_sd2taylor(sd, correl) = (x=correl * sd, y=sqrt(sd^2-(correl*sd)^2))

    ## Create Taylorplot gridlines
    correls = reduce(vcat,[-1, -0.99, -0.95, -0.9:0.1:0.9, 0.95, 0.99, 1.0])
    sds = 0:0.1:sdmax |> collect

    fig, ax, sc = scatter(0,0; color=(:black, 0), kw...);
    for c in correls
        lines!(ax, c .* sds .* sdObs * 1.02, sqrt.(sds.^2-(c .* sds).^2).*sdObs * 1.02; linestyle= :dot, color=:grey) 
    end

    for s in sds
        lines!(ax, correls .* s .* sdObs, sqrt.(s.^2 .- (correls .* s).^2).*sdObs; linestyle= s==1 ? :dash : :dot, color=:grey) 
    end

    ## Creat isolines of root mean squared difference (RMSD) and labels
    RMSD = 0.2:0.2:2*sdmax+0.2 |> collect
    xvec =  -1.0:0.01:1 |> collect

    rmsd_lab_x=Vector{Float64}()
    rmsd_lab_y=Vector{Float64}()
    rmsd_lab_text=Vector{String}()

    for r in RMSD[1:end-1]
         
        x2 = @. (xvec * r + 1) * sdObs
        y2 = @. sqrt(1-xvec^2) * r * sdObs
        
       i_plot = [i for i in 1:length(x2) if sqrt(x2[i]^2+y2[i]^2)<=1.005*maximum(sds)*sdObs]

      # This also works, but not if the array sent back by the compr is empty (collect does not work)  
      # xp, yp = zip([(xi,yi) for (xi, yi) in zip(x2, y2) if sqrt(xi^2+yi^2)<=1.005*maximum(sds)*sdObs]...) |> collect
       
        lines!(ax, x2[i_plot], y2[i_plot], color=(:blue, 0.5))

        i_lab=[i for i in 1:length(x2) if (y2[i]>0.3^2*sdObs^2/abs(x2[i])) & (x2[i]>0.0) & (sqrt(x2[i]^2+y2[i]^2)<=1.005*maximum(sds)*sdObs)]
        if length(i_lab)>0
            push!(rmsd_lab_x, x2[i_lab[1]])
            push!(rmsd_lab_y, y2[i_lab[1]]) 
            push!(rmsd_lab_text, "$(round(r, digits=2))")
         end
       
    end
    
    #### Plotting data in Taylor space

    x,y = collect.(zip(correl_sd2taylor.(sdmod_rel, correl)...) |> collect) .* sdObs

    # First the arrows, so that the point is always on top
    if plotBias
        RMSD = @. sqrt(sdObs^2 + (sdmod_rel*sdObs)^2 - 2*sdObs^2*sdmod_rel*correl)
        RMSE = @. sqrt(bias^2+RMSD^2)
        RMSDangle = @. asin(y/RMSD)
        RMSEangle = @. asin(bias/RMSE)
        alpha = @. ifelse(x > sdObs, π - RMSDangle - RMSEangle, RMSDangle - RMSEangle)
        xb = @. sdObs - cos(alpha)* RMSE
        yb = @. sin(alpha) * RMSE
        arrows!(x, y, xb.-x, yb.-y; color=1:length(x), linewidth=3)
    end

    if plotVarianceErr
        x_no_var_err, y_no_var_err = collect.(zip(correl_sd2taylor.(1.0, correl)...) |> collect) .* sdObs
        arrows!(x,y, x_no_var_err.-x, y_no_var_err.-y; color=1:length(x), linewidth=3)
    end

    # if (plotNoVarErr) {
    #     modPointsNoVarErr <- correl_sd2Taylor(1, correl) * sdObs
    #     p  <- p + geom_point(data=modPointsNoVarErr, mapping=aes(x,y), color="red", size=2, alpha=0.5)

    ## Plot data points (original Taylorplot)
    #colors=[:red, :blue, :green, :brown, :yellow, :orange]
    
    allScat = map(eachindex(x)) do i
        scatter!(ax, x[i], y[i]; strokewidth=0.5, glowwidth = 5.0, glowcolor=(:black,0.5), label=labels[i], markersize)
    end
    ## Plot point of "perfect model" 
    scatter!(ax, sdObs, 0.0, markersize=20, color=(:blue,0.5))
    
    Legend(fig[1, 2], allScat, labels)
   
     
    ### Labels for coordinate system
    correlTicks = @chain reduce(vcat, [-0.9:0.1:0.9, 0.95, 0.99]) setdiff(_, [0])
    text!(correlTicks*sdmax*sdObs*1.04, sqrt.(sdmax^2 .- (correlTicks*sdmax).^2)*sdObs*1.04; 
        text=["$(round(c, digits=2))" for c in correlTicks], align=(:center,:center), fontsize=20)
    #scatter!(correlTicks*sdmax*sdObs*1.04, sqrt.(sdmax^2 .- (correlTicks*sdmax).^2)*sdObs*1.04)
    scatter!(rmsd_lab_x,rmsd_lab_y; markersize=(60,40), color=:white, strokewidth=0.5)
    text!(rmsd_lab_x,rmsd_lab_y; text=rmsd_lab_text, align=(:center,:center), fontsize=20)

    ## Clipping plot
    xlims!(ax, min(0, 1.1*minimum(x)), max(sdObs,1.1*maximum(x)))

    (fig, ax)



    # taylor_grid = @pipe Base.product(correls, sds) |> DataFrame |> rename(_, [:correls, :sds])

    # iso_rmsd = @pipe Base.product(0.2:0.2:2*sdmax+0.2 |> collect, -1.0:0.01:1 |> collect) |> 
    #     DataFrame |> rename(_, [:RMSD, :xvec])

    # @chain iso_rmsd begin
    #     @rtransform! begin
    #         :x2 = (:xvec * :RMSD + 1)* sdObs
    #         :y2 = sqrt(1-:xvec^2)*:RMSD*sdObs
    #     end
    #     @rsubset! sqrt(:x2^2+:y2^2)<=1.005*maximum(taylor_grid.sds)*sdObs
    # end
    
    # @rtransform! taylor_grid :xv = :correls*:sds
    # @rtransform! taylor_grid begin
    #     :yv = sqrt(:sds^2-:xv^2)*sdObs
    #     :xv=:xv*sdObs 
    # end
    # @transform! taylor_grid begin
    # :sdsCat = categorical(:sds)
    # :correlCat = categorical(:correls)
    # end

    # #plt = data(taylor_grid) * visual(Lines, linestyle=:dash) *  mapping(:xv, :yv, group=(:correlCat, :sdsCat))
    # #plt += data(taylor_grid) * visual(Lines, linestyle=:dash) *  mapping(:xv, :yv, group=:sdsCat)
    # #draw(plt; axis=(width = 1225, height = 825))
    # fig=Figure()


end

qplot(df::DataFrame, x, y, geom=Scatter; color=:blue, kw...) = data(df) * mapping(x,y) * visual(geom; color )
qplot(x, y, geom=Scatter) = (df=DataFrame(;x, y); qplot(df, :x, :y, geom ))