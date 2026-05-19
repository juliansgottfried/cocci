using Plots, JLD2, StatsBase

include("estimate.jl")
include("generate.jl")

dρ = 0.5
maxρ = 20 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ
J = 500

# Retrodiction plots

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
    title = i == 1 ? "ρ0" : "ρ1"
    dat = val[:, :, i, 2]
    if i == 2 dat = dat' end
    plot(ρs, dat,
        xlim = [0, maxρ], ylim = [0, maxρ],
        xlabel = "recombination rate", ylabel = "estimate",
        alpha = 0.5, color = :black, linewidth = 0.5,
        title = title, label = false, grid = false)
    plot!([0, maxρ], [0, maxρ], 
        color = :red, label = false)
end

makeheat = function(i)
    diff = zeros(Float64, nρ, nρ)
    compρs = collect(ρs)
    compρs[1] = dρ
    if i == 1 
        for j in 1:nρ
            diff[j, :] = log2.((val[j, :, i, 2] .- compρs[j]) .^ 2 ./ compρs[j])
        end
    else
        for j in 1:nρ
            diff[:, j] = log2.((val[:, j, i, 2] .- compρs[j]) .^ 2 ./ compρs[j])
        end
    end
    if i != 1 diff = diff' end
    diff
end

plotheat = function(input, i)
    title = i == 1 ? "ρ0 estimate" : "ρ1 estimate"
    heatmap(ρs, ρs, input,
            xlabel = "ρ0", ylabel = "ρ1",
            xlim = [0, maxρ], ylim = [0, maxρ],
            clim = (0, topval),
            label = "sq", c = :thermal, title = title,
            aspect_ratio = 1, colorbar_title = "relative error")
end

plotit(1)
plotit(2)

heat0 = makeheat(1)
heat1 = makeheat(2)
topval = max(maximum(heat0), maximum(heat1))

plotheat(heat0, 1)
plotheat(heat1, 2)

# if rho0 is large, it estimates an equally large rho1
# it says that the highest likelihood is for
# both very large

# rho0 is accurately estimated for all values
# rho1 is accurate estimated if rho0 is low
# but if rho0 is high, it thinks that rho1 is also how


# i mean it looks like the intercept of the rho1 estimate
# is in fact the rho0 estimate
# and then it goes up
# for some reason the estimator for rho1 is actually rh0+rho1