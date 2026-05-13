using Plots
using JLD2
import Statistics

include("estimate.jl")
include("generate.jl")

dρ = 0.5
maxρ = 10
nρ = Int(div(maxρ, dρ)) + 1
# collect = [generate.getobjlocal("data", "5_11_26", dρ, idx1, idx2) for idx1 in 1:nρ, idx2 in 1:nρ]

collect = [generate.getobjlocal("prob", "5_11_26", dρ, idx1, idx2) for idx1 in 1:nρ, idx2 in 1:nρ]

collect[1, :]

idx1 = 20
idx2 = 20
ρ0 = dρ * (idx1 - 1)
ρ1 = dρ * (idx2 - 1)
histogram(collect[idx1, idx2][:, 1], bins = 20,
    color = :white, fill = :grey, label = false)
histogram(collect[idx1, idx2][:, 2], bins = 20,
    color = :white, fill = :grey, label = false)

Statistics.quantile(collect[idx1, idx2][:, 1], [0.1, 0.5, 0.9])
Statistics.quantile(collect[idx1, idx2][:, 2], [0.1, 0.5, 0.9])

plotter = zeros(Float64, nρ, nρ, 3)
for i in 1:nρ, j in 1:nρ
    ρ0 = dρ * (i - 1)
	ρ1 = dρ * (j - 1)
    plotter[i, j, :] = Statistics.quantile.(log2.(collect[i, j] .- [ρ0; ρ1]), [0.055, 0.5, 0.945])
end
plotter[plotter .== -Inf] .= minimum(plotter[plotter .!= -Inf])

plot(ρs[2:end], plotter[:, 1], fillrange = plotter[:, 3], 
    xlabel = "recombination rate ρ", ylabel = "log₂ relative error",
    xlim = [minimum(ρs), maximum(ρs)], fillalpha = 0.15, c = "#29a0c8",
    linecolor = false, label = false, grid = false)
hline!([0], c = :red, alpha = 0.7, label = false)
plot!(ρs[2:end], plotter[:, 2],
    c = :black, linewidth = 1.2, label = false)

