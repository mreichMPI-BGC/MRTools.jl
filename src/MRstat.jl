using StatsBase, Statistics, Loess, GLM, LsqFit, ForwardDiff
using SplitApplyCombine, TensorCast
using StaticModules: @with
using ImageFiltering
using Zygote: @adjoint
using RCall: @R_str, rcopy
### Loss functions
function qloss(yp, y, prob=0.5f0)
   # mean(((y .< yp) .- prob) .* (yp .- y))
   dev = yp .- y 
   mean( max.(-prob .* dev, (1-prob) .* dev) )
end

qlossfun(prob) = (yp, y) -> qloss(yp, y, prob)

@adjoint function sort(x)
    p = sortperm(x)
    x[p], x̄ -> (x̄[invperm(p)],)
  end

function trimmedLoss(yp,y, keepFrac=1.0f0)
    absRes = abs2.(yp .- y) |> vec
    trimVal = quantile(sort(absRes), keepFrac, sorted=true) 
    return sum(ifelse.(absRes .> trimVal,  0.f0 , absRes ))/(length(absRes)*(keepFrac))
    #s = sum(absRes)/length(absRes)   # using this and commenting out the two above works (no surprise)    
end

trimmedLossfun(keepFrac) = (yp, y) -> trimmedLoss(yp, y, keepFrac)

function save_quantile(x, p=[0.0,0.025,0.25,0.5,0.75, 0.975,1.0]; kw...)
    xok = filter(oknum, x)
    #xok = x |> skipmissing
    isempty(xok) && return fill(NaN, length(p))

    return quantile(xok, p; kw...)
end

function trendOfVec(y)
    xy = DataFrame(; x = eachindex(y) |> collect, y=y |> collect)
    lreg = lm(@formula(y~x), xy)
    coef(lreg)[2]
end


function trend_plot!(x, y, trend=mr_linear, ax=current_axis(); color=:black,  kw...)
    df = trend(x,y; kw...)
    @with df begin
        band!(ax, x, lower, upper; color=(color, 0.3), kw...)
        lines!(ax, x, prediction; color=color, kw...)
    end
end

function mr_linear(x, y)
    xy =  DataFrame(; x, y) |> dropmissing |> disallowmissing
    lreg = lm(@formula(y~x), xy)
    newx = DataFrame(x = range(extrema(xy.x)..., 50));
    pr = GLM.predict(lreg, newx, interval = :confidence, level = 0.95)
    return hcat(newx, pr)
    
end

function mr_loess(x, y; kw...)
    
    xy = DataFrame(; x, y) |> dropmissing |> disallowmissing
    newx = DataFrame(x = range(extrema(x)..., 100));
    model = loess(xy.x, xy.y; span=0.3, degree=1, kw...)

    prediction=Loess.predict(model, newx.x)

    y_hat = Loess.predict(model, xy.x)
    res = xy.y - y_hat

    oneLoessFit(x, y_boot) = (boot_model=Loess.loess(x, y_boot; span=0.3, degree=1, kw...); return Loess.predict(boot_model, newx.x))

    allFits = (oneLoessFit(xy.x, y_hat .+ sample(res, length(res), replace=true)) for i in 1:250) |> collect |> invert
    #@cast reshap[i][j] := allFits[j][i]
    #allFits = reduce(hcat, allFits)
    #allFits = reinterpret(allFits,  )
    qunts=[quantile(x, [0.025, 0.975]) for x in allFits] |> invert

    return DataFrame(; x=newx.x, prediction, lower=qunts[1], upper=qunts[2])
    
end

# Calculate confidence or prediction intervals from LsqFit (from https://github.com/JuliaNLSolvers/LsqFit.jl/pull/180/files/8ff8e30bbc34f9e302a71e8445fc140259fe35c8)
function lsqFitBands(model::Function, x, fit::LsqFit.LsqFitResult, alpha=0.05; purpose=:confidence)
    g = p -> ForwardDiff.gradient(p -> model(x, p), fit.param)
    c = g(fit.param)' * estimate_covar(fit) * g(fit.param)
    if purpose === :prediction
        c = c + 1
    end
    dist = TDist(dof(fit))
    critical_values = quantile(dist, 1 - alpha/2)
    return sqrt(c*rss(fit)/dof(fit))*critical_values
end

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

function trimReplace!(x; trimFrac=0.0, mode="NaN", kwargs...)
    if trimFrac > 0.0
        qunts = quantile(filter(!isnan, x), (trimFrac / 2.0, 1.0 - trimFrac / 2.0))
        #print(qunts)
        if mode == "NaN"
            map!(x -> (x < qunts[1]) | (x > qunts[2]) ? NaN : x, x, x)
        else
            map!(x -> max.(min.(x, qunts[2]), qunts[1]), x, x)
        end
    end
end
  
function trimReplace(x; trimFrac=0.0, kwargs...)
    xx = copy(x)
    trimReplace!(xx; trimFrac=trimFrac, kwargs...)
    #println(count(isnan, xx))
    xx

end
  
function loess4Miss(x, y; span=0.3, useR=false, kwargs...)
    if count(!isnan, y) < 3
        return (similar(y, Missing))
    end
    if useR
    else

        i = 1
        while isnan(y[i])
            i += 1
        end
        y[1:i] .= y[i]

        i = 0
        while isnan(y[end-i])
            i += 1
        end
        y[end-i:end] .= y[end-i]

        df = dropmissing(DataFrame(x=x, y=replace(y, NaN => missing), copycols=false), :y)
        model = loess(df.x, df.y; span=span, degree=1)
        Loess.predict(model, x)
    end
end


regBivar(x,y) = regBivar(x,y,lm)
#regBivar(x,y, fun) = @c fun([fill(1, length(x)) x], y) coef() NamedTuple(zip(string(fun) .* ["_inter", "_slope"] .|> Symbol, _ ))
regBivar(x,y, fun) = @chain begin
    DataFrame(x=x,y=y;copycols=false)
    fun(@formula(y ~ 1 + x), _)
    coef
    NamedTuple(zip(string(fun) .* ["_inter", "_slope"] .|> Symbol, _ ))
end

mapwindow2(fun, x,y, win)= mapwindow(x-> fun(getindex.(x,1), getindex.(x,2)), zip(x,y)|>collect, win)

function harmonics(period::Number, nwave::Number, fun=sin)
    function h(x)
        return fun(2π * nwave * x / period)
    end
    return h
end

function generate_harmonics(v, nwaves)
     mx = maximum(filter(oknum, v))
     res_sin = map(n->harmonics(mx, n, sin).(v), nwaves)
     res_cos = map(n->harmonics(mx, n, cos).(v), nwaves)

     names = [string.("cos_", nwaves)..., string.("sin_", nwaves)...]

     @chain reduce(hcat, append!(res_cos, res_sin)) DataFrame(names)
end