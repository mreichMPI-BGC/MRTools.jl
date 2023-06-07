#using Missings

catchCondition(f; condition=isempty,  default=missing) = x -> condition(x) ? default : f(x)
catchAny(f; condition=isempty,  default=missing) = x -> condition(x) ? default : try f(x) catch e  default end
catchMissing(f; catchFun=catchCondition) = x -> catchFun(f)(skipmissing(x))

function extend_dims(A,which_dim) #from https://stackoverflow.com/questions/42312319/how-do-i-add-a-dimension-to-an-array-opposite-of-squeeze
    s = [size(A)...]
    insert!(s, which_dim, 1)
    return reshape(A, s...)
end

# Flatten a nested array (from discourse)
flat(arr) = mapreduce(x -> x == [] || x[1] === x ? x : flat(x), vcat, arr, init=[])


oknum(x::Union{Number, Missing}) = !ismissing(x) && isfinite(x)
oknum(x::String) = !ismissing(x) && !isequal(x, "")

nok2miss(x::Union{Number, Missing}) = oknum(x) ? x : missing
nok2miss!(df::DataFrame) = @transform! df {} = nok2miss({names(df, Number)})
miss2nan!(df::DataFrame) =  @transform!(df, {} = coalesce({names(df, Union{Missing, Number})}, NaN)) |> disallowmissing!

ok_frac(x::AbstractArray) = sum(oknum.(x))/length(x)

function expand_grid(; kws...)
    names, vals = keys(kws), values(kws)
    return DataFrame(NamedTuple{names}(t) for t in Iterators.product(vals...))
end


function cutMiddleAndReplace(colvec, cutfract; repl="grey80")
    n = length(colvec)
    iseven = n % 2 == 0
    halfwidth = trunc(Int, n * cutfract * 0.5)
    cent = (n+1) ÷ 2
    outVec = vcat(colvec[1:(cent - halfwidth + iseven - 1)], colvec[(cent + halfwidth + 1):end])
    insert!(outVec, cent - halfwidth + iseven, repl)
    return outVec
end

function nest(df::DataFrame, group_cols; nest_cols = nothing)
    nest_cols = isnothing(nest_cols) ? names(df[!, Not(group_cols)]) : nest_cols
    combine(groupby(df, group_cols),
          nest_cols => ((args...) -> Ref(DataFrame([nest_cols[i] => a for (i,a) ∈ enumerate(args)]))) => :nested)
end

function unnest_dfcol(df::DataFrame, col; appendName=true)
    #[DataFrame(gr = reduce(vcat, fill.(df_nested.gr, nrow.(df_nested.X_DataFrame)))) reduce(vcat, df_nested.X_DataFrame)]
    df_gr = select(df, Not(col))
    gr_expand = @chain [reduce(vcat, fill(DataFrame(g), nrow(n[col])))   for (g,n)  in zip(eachrow(df_gr), eachrow(df))] reduce(vcat, _)
    unnested = reduce(vcat, df[:,col] )
    appendName && rename!(x-> string(col) * "_" * x, unnested)

    #return hcat(gr_expand, unnested)
    return [gr_expand unnested]
end

function unnest(df::DataFrame, col; appendName=true)
    dfNew = select(df, col => AsTable) 
    appendName && rename!(x-> string(col) * "_" * x, dfNew)
    return hcat(dfNew, select(df,Not(col)), makeunique=true)
end