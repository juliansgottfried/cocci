using Plots, JLD2, StatsBase, Statistics
import Distributions

include("estimate.jl")
include("generate.jl")

dρ = 0.1
maxρ = 20 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ

J = 500

# Retrodiction plots

loadit(filename) = if isfile(filename) return load_object(filename) end

data0 = [loadit(generate.getfilenamelocal("data", "5_17_26_d", true, ρ)) for ρ in 0:dρ:maxρ]
data1 = [loadit(generate.getfilenamelocal("data", "5_17_26_d", false, ρ)) for ρ in 0:dρ:maxρ]
skip = isnothing.(data0) .| isnothing.(data1)

aic = zeros(Float64, nρ, 2, 3)
raw = zeros(Float64, nρ, 4, 3)
qs = [0.055, 0.5, 0.945]
for i in 1:nρ
    if skip[i] 
        aic[i, :, :] .= Inf
        raw[i, :, :] .= Inf
    else
        aic[i, 1, :] = quantile(2(data0[i][:, 2] - data0[i][:, 4]), qs)
        aic[i, 2, :] = quantile(2(data1[i][:, 4] - data1[i][:, 2]), qs)
        # aic[i, 1, 2] = Statistics.mean(2(data0[i][:, 2] - data0[i][:, 4]))
        # aic[i, 2, 2] = Statistics.mean(2(data1[i][:, 4] - data1[i][:, 2]))

        raw[i, 1, :] = quantile(data0[i][:, 1], qs)
        raw[i, 2, :] = quantile(data0[i][:, 3], qs)
        raw[i, 3, :] = quantile(data1[i][:, 1], qs)
        raw[i, 4, :] = quantile(data1[i][:, 3], qs)
    end
end

aicplot = function(i, lbound, rbound)
    ploterror = aic[:, i, :]
    ploterror[ploterror .< lbound] .= lbound
    ploterror[ploterror .> rbound] .= rbound
    plot(ρs, ploterror[:, 1], fillrange = ploterror[:, 3], 
        xlabel = "recombination rate", ylabel = "Δ AIC",
        title = i == 1 ? "constant data" : "varying data",
        xlim = [0, maxρ], ylim = [lbound, rbound],
        fillalpha = 0.15, c = "#29a0c8",
        linecolor = false, label = false, grid = false)
    hline!([2], c = :red, alpha = 0.7, label = false)
    plot!(ρs, aic[:, i, 2], c = :black, linewidth = 1.2, label = false)
end

rawplot = function(i, j)
    k = 2(i - 1) + j
    part1 = i == 1 ? "constant data" : "varying data"
    part2 = j == 1 ? "constant model" : "varying model"
    plot(ρs, raw[:, k, 1], fillrange = raw[:, k, 3], 
        xlabel = "recombination rate", ylabel = "estimate",
        title = "$part1, $part2",
        xlim = [0, maxρ], fillalpha = 0.15, c = "#29a0c8",
        linecolor = false,
        label = false, grid = false)
    plot!(ρs, maxρ * ones(nρ),
        linestyle = :dash, color = "#29a0c8", label = false)
    plot!(ρs, zeros(nρ), 
        linestyle = :dash, color = "#29a0c8", label = false)
    plot!(ρs, raw[:, k, 2], c = :black, 
        linewidth = 1.2, label = false)
    plot!([0, maxρ], [0, maxρ], 
        color = :red, label = false)
end

aicplot(1, -100, 100)
aicplot(2, -100, 100)

rawplot(1, 1)
rawplot(1, 2)
rawplot(2, 1)
rawplot(2, 2)


# For Bayes

getp = function(data, i, ρs, J)
    counts = data1[i][:, 3]
    counts = round.(sort(counts), sigdigits = 3)
    p = [sum(counts .== ρ) / J for ρ in ρs]
end

p0 = getp(data1, 1, ρs, J)

pc0 = zeros(nρ)
for i in 2:nρ
    pc0 .+= getp(data1, i, ρs, J) ./ (nρ - 1)
end

ll0 = mean(ρhat[:, 7])
llc0 = mean(ρhat[:, 8])


# P(0 | rho_hat distr) =
# P(rho_hat dist | 0) * L(0)
# divided by
# P(rho_hat dist | 0) * L(0) + P(rho_hat dist | >0) * L(>0)

ll0 = mean(ρhat[:, 7])
llc0 = mean(ρhat[:, 8])

# what is P(rho_hat dist | 0)
# multinomial distribution! perfect

# and P(rho_hat dist | >0) ?
# just the average of all those

# easy
# do tomorrow!