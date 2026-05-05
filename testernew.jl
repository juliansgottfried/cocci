using Plots
using Measures
using JLD2
import Statistics

include("estimate.jl")
include("generate.jl")

ρs = 0:0.1:15
filenames = string.("./results/data/run_5_5_26/results_", replace.(string.(ρs), "." => "_"), ".jld2")
collect = []
for i in filenames
    push!(collect, load_object(i)[1])
end

plotter = zeros(Float64, length(ρs) - 1, 3)
for i in eachindex(ρs)[2:end]
    plotter[i - 1, :] = Statistics.quantile(log2.(collect[i] / ρs[i]), [0.055, 0.5, 0.945])
end
plotter[plotter .== -Inf] .= minimum(plotter[plotter .!= -Inf])

plot(ρs[2:end], plotter[:, 1], fillrange = plotter[:, 3], 
    xlabel = "recombination rate ρ", ylabel = "log₂ ρ hat / ρ",
    title = "Estimation error",
    xlim = [0, 15], fillalpha = 0.15, c = "#29a0c8",
    linecolor = false, label = false, grid = false)
hline!([0], c = :red, alpha = 0.7, label = false)
plot!(ρs[2:end], plotter[:, 2],
    c = :black, linewidth = 1.2, label = false)
