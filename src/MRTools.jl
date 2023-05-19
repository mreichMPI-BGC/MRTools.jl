module MRTools

using Reexport
@reexport using RCall
# @reexport using Bootstrap
# @reexport using CategoricalArrays
# @reexport using Clustering
@reexport using CSV
@reexport using DataFrames
@reexport using NCDatasets
# @reexport using Distances
@reexport using Distributions
@reexport using GLM
# @reexport using HypothesisTests
# @reexport using KernelDensity
@reexport using Loess
# @reexport using MixedModels
# @reexport using MultivariateStats
# @reexport using ShiftedArrays
@reexport using Statistics
@reexport using StatsBase
# @reexport using TimeSeries

@reexport using ColorSchemes, Colors
@reexport using LsqFit, ForwardDiff
@reexport using SplitApplyCombine, TensorCast
@reexport using StaticModules: @with
@reexport using ImageFiltering
@reexport using Zygote: @adjoint
@reexport using Chain: @chain
#@reexport using RCall: @R_str, rcopy


# Write your package code here.
foreach(include,  ["Misc.jl", "MR_geo.jl", "MRplot.jl", "MR_io.jl", "MRstat.jl"])

for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end