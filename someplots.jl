using Plots
using JLD2
import Statistics

include("estimate.jl")
include("generate.jl")

dρ = 0.1
maxρ = 10 - dρ
nρ = floor(Int, maxρ / dρ) + 1

data0 = [load_object(generate.getfilenamelocal("data", "5_14_26", true, ρ)) for ρ in 0:dρ:maxρ]
data1 = [load_object(generate.getfilenamelocal("data", "5_14_26", false, ρ)) for ρ in 0:dρ:maxρ]

ρs = 0:dρ:maxρ

aic = zeros(Float64, nρ, 2, 3)
ratio = zeros(Float64, nρ, 4, 3)
qs = [0.055, 0.5, 0.945]
for i in 1:nρ
    # i = 1
    aic[i, 1, :] = Statistics.quantile(data0[i][:, 2] - data0[i][:, 4], qs)
    aic[i, 2, :] = Statistics.quantile(data1[i][:, 2] - data1[i][:, 4], qs)
    if ρs[i] == 0 
        ratio[i, 1, :] = [Inf; Inf; Inf]
        ratio[i, 2, :] = [Inf; Inf; Inf]
        ratio[i, 3, :] = [Inf; Inf; Inf]
        ratio[i, 4, :] = [Inf; Inf; Inf]
    else
        ratio[i, 1, :] = Statistics.quantile(log2.(data0[i][:, 1] ./ ρs[i]), qs)
        ratio[i, 2, :] = Statistics.quantile(log2.(data0[i][:, 3] ./ ρs[i]), qs)
        ratio[i, 3, :] = Statistics.quantile(log2.(data1[i][:, 1] ./ ρs[i]), qs)
        ratio[i, 4, :] = Statistics.quantile(log2.(data1[i][:, 3] ./ ρs[i]), qs)
    end
end
ratio[ratio .== -Inf] .= minimum(ratio[ratio .!= -Inf])
ratio[ratio .== Inf] .= maximum(ratio[ratio .!= Inf])

aicplot = function(i)
    plot(ρs, aic[:, i, 1], fillrange = aic[:, i, 3], 
        xlabel = "recombination rate ρ", ylabel = "Δ AIC",
        title = i == 1 ? "constant data" : "varying data",
        xlim = [0, maxρ], fillalpha = 0.15, c = "#29a0c8",
        linecolor = false, label = false, grid = false)
    hline!([1], c = :red, alpha = 0.7, label = false)
    plot!(ρs, aic[:, i, 2], c = :black, linewidth = 1.2, label = false)
end

ratioplot = function(i, j)
    k = 2(i - 1) + j
    part1 = i == 1 ? "constant data" : "varying data"
    part2 = j == 1 ? "constant model" : "varying model"
    plot(ρs, ratio[:, k, 1], fillrange = ratio[:, k, 3], 
        xlabel = "recombination rate ρ", ylabel = "log2 relative error",
        title = "$part1, $part2",
        xlim = [0, maxρ], fillalpha = 0.15, c = "#29a0c8",
        linecolor = false, label = false, grid = false)
    hline!([0], c = :red, alpha = 0.7, label = false)
    plot!(ρs, log2.(maxρ ./ ρs), color = :red, label = false)
    plot!(ρs, ratio[:, k, 2], c = :black, linewidth = 1.2, label = false)
end

aicplot(1)
aicplot(2)

ratioplot(1, 1)
ratioplot(2, 1)
ratioplot(1, 2)
ratioplot(2, 2)
