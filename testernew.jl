using Plots
using JLD2

include("estimate.jl")
include("generate.jl")

ρs = 0:0.1:15
filenames = string.("./results/data/run_5_4_26/results_", replace.(string.(ρs), "." => "_"), ".jld2")
collect = []
for i in filenames
    push!(collect, load_object(i)[1])
end

plotter = zeros(Float64, length(ρs) - 1, 3)
for i in eachindex(ρs)[2:end]
    plotter[i, :] = quantile(log(collect[i] / ρs[i]), [0.055, 0.5, 0.945])
end

plot(ρs[2:end],  plotter[:, 1], fillrange = plotter[:, 3], 
    xlabel = "ρ", ylabel = "log % error",
    fillalpha = 0.15, c = "#29a0c8", linecolor = false,
    label = false, grid = false)
plot!(ρs[2:end], plotter[:, 2],
    c = :black, label = false)