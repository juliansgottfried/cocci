using Plots, JLD2, StatsBase

include("estimate.jl")
include("generate.jl")

dρ = 1
maxρ = 10 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ
J = 500

# Retrodiction plots

loadit(filename) = if isfile(filename) return load_object(filename) end

data = [load_object(generate.getfilenamegridlocal("data", "5_18_26_e", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]
skip = isnothing.(data)

data[1,1]

val = zeros(Float64, nρ, nρ, 2, 3)
qs = [0.055, 0.5, 0.945]
for i in 1:nρ, j in 1:nρ
    if skip[i, j] 
        val[i, j, :, :] .= Inf
    else
        tmpdat = data[i, j][:, 1:2]
        val[i, j, 1, :] = quantile(tmpdat[:, 1], qs)
        val[i, j, 2, :] = quantile(tmpdat[:, 2], qs)
    end
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

medianρ0 = val[:, :, 1, 2]
medianρ1 = val[:, :, 2, 2]'

plot(ρs, medianρ0,
    xlabel = "recombination rate", ylabel = "estimate",
    title = "ρ0",
    color = :black, label = false)
plot(ρs, medianρ1, 
    xlabel = "recombination rate", ylabel = "estimate",
    title = "ρ1",
    color = :black, label = false)
