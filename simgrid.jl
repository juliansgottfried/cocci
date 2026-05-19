using Plots, JLD2, StatsBase

include("estimate.jl")
include("generate.jl")

dρ = 1
maxρ = 10 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ

# Retrodiction plots

loadit(filename) = if isfile(filename) return load_object(filename) end

data = [load_object(generate.getfilenamegridlocal("data", "5_19_26_a", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]
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

plotit = function(i)
    title = i == 1 ? "intercept" : "slope"
    dat = i == 1 ? val[:, :, i, 2] : val[:, :, i, 2]'
    plot(ρs, dat,
        xlim = [0, maxρ], ylim = [0, maxρ],
        xlabel = "recombination rate", ylabel = "estimate",
        alpha = 0.5, color = :black, linewidth = 0.5,
        title = title, label = false, grid = false)
    plot!([0, maxρ], [0, maxρ], 
        color = :red, label = false)
end

plotit(1)
plotit(2)
