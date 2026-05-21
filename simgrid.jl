using Plots, JLD2, StatsBase
import Distributions

include("estimate.jl")
include("generate.jl")

dρ = 0.5
maxρ = 20 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ
J = 500

loadit(filename) = if isfile(filename) return load_object(filename) end

data = [loadit(generate.getfilenamegridlocal("data", "5_19_26_a", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]
skip = isnothing.(data)

val = zeros(Float64, nρ, nρ, 2, 3)
qs = [0.055, 0.5, 0.945]
for i in 1:nρ, j in 1:nρ
    if skip[i, j] 
        val[i, j, :, :] .= -Inf
    else
        tmpdat = data[i, j][:, 1:2]
        val[i, j, 1, :] = quantile(tmpdat[:, 1], qs)
        val[i, j, 2, :] = quantile(tmpdat[:, 2], qs)
    end
end

plotit = function(i)
    xlab = i == 1 ? "ρ0" : "ρ1"
    dat = val[:, :, i, 2]
    if i == 2 dat = dat' end
    plot(ρs, dat,
        xlim = [0, maxρ], ylim = [0, maxρ],
        xlabel = xlab, ylabel = "estimate",
        alpha = 1, color = :black, linewidth = 0.4,
        label = false, grid = false)
    plot!([0, maxρ], [0, maxρ], 
        color = :red, label = false)
end

plotheat = function(i, cutoff)
    title = i == 1 ? "ρ0 estimate" : "ρ1 estimate"
    diff = zeros(Float64, nρ, nρ)
    compρs = collect(ρs)
    compρs[1] = dρ
    if i == 1 
        for j in 1:nρ
            diff[j, :] = (val[j, :, i, 2] .- compρs[j]) ./ compρs[j]
        end
    else
        for j in 1:nρ
            diff[:, j] = (val[:, j, i, 2] .- compρs[j]) ./ compρs[j]
        end
    end
    if i != 1 diff = diff' end

    diff[.!isfinite.(diff)] .= 0

    diff[diff .> cutoff] .= cutoff
    diff[diff .< -cutoff] .= -cutoff


    heatmap(ρs, ρs, diff,
            xlabel = "ρ0", ylabel = "ρ1",
            xlim = [0, maxρ], ylim = [0, maxρ],
            clim = (-cutoff, cutoff),
            label = "sq", c = cgrad(:diff, rev = false), 
            title = title,
            aspect_ratio = 1, colorbar_title = "relative error")
end

plotit(1)
plotit(2)

cutoff = 2
plotheat(1, cutoff)
plotheat(2, cutoff)
