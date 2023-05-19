#using CairoMakie
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

one2one(; color=:grey, linestyle=:dash, kw...) =  mapping([0], [1]) * visual(ABLines; color, linestyle, kw...)

function evalplot(df, model::Symbol, obs::Symbol; kw...)
  data(df) * mapping(model, obs) * (visual(Scatter; kw... ) + linear()) + one2one()
end