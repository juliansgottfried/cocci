using Plots
using JLD2
import Statistics

include("estimate.jl")
include("generate.jl")

dρ = 0.1
maxρ = 10
nρ = Int(div(maxρ, dρ)) + 1

data0 = [load_object(generate.getfilenamelocal("data", "5_13_26", true, ρ)) for ρ in 0:dρ:maxρ]
data1 = [load_object(generate.getfilenamelocal("data", "5_13_26", false, ρ)) for ρ in 0:dρ:maxρ]

aic = zeros(Float64, nρ, 2, 3)
ratio = zeros(Float64, nρ, 4, 3)
for i in 1:nρ
    aic[i, 1, :] = Statistics.quantile.(data0[i][:, 2] - data0[i][:, 4], 
        [0.055, 0.5, 0.945])
    aic[i, 2, :] = Statistics.quantile.(data1[i][:, 4] - data0[i][:, 2], 
        [0.055, 0.5, 0.945])
    ρ = 0:dρ:maxρ
    ratio[i, 1, :] = Statistics.quantile.(data0[i][:, 1] ./ ρ[i], [0.055, 0.5, 0.945])
    ratio[i, 2, :] = Statistics.quantile.(data0[i][:, 3] ./ ρ[i], [0.055, 0.5, 0.945])
    ratio[i, 3, :] = Statistics.quantile.(data1[i][:, 1] ./ ρ[i], [0.055, 0.5, 0.945])
    ratio[i, 4, :] = Statistics.quantile.(data1[i][:, 3] ./ ρ[i], [0.055, 0.5, 0.945])
end
ratio[ratio .== -Inf] .= minimum(ratio[ratio .!= -Inf])
ratio[ratio .== Inf] .= maximum(ratio[ratio .!= Inf])

aicplot = function(i)
    plot(ρs, aic[:, i, 1], fillrange = aic[:, i, 3], 
        xlabel = "recombination rate ρ", ylabel = "Δ AIC",
        xlim = [0, maxρ], fillalpha = 0.15, c = "#29a0c8",
        linecolor = false, label = false, grid = false)
    hline!([1], c = :red, alpha = 0.7, label = false)
    plot!(ρs, aic[:, i, 2], c = :black, linewidth = 1.2, label = false)
end

ratioplot = function(i)
    plot(ρs, ratio[:, i, 1], fillrange = ratio[:, i, 3], 
        xlabel = "recombination rate ρ", ylabel = "relative error",
        xlim = [0, maxρ], fillalpha = 0.15, c = "#29a0c8",
        linecolor = false, label = false, grid = false)
    hline!([1], c = :red, alpha = 0.7, label = false)
    plot!(ρs, ratio[:, i, 2], c = :black, linewidth = 1.2, label = false)
end

aicplot(1)
aicplot(2)

ratioplot(1)
ratioplot(2)
ratioplot(3)
ratioplot(4)