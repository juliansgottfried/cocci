using Plots, JLD2, StatsBase

include("estimate.jl")
include("generate.jl")

dρ = 0.5
maxρ = 20 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ

loadit(filename) = if isfile(filename) return load_object(filename) end
data = [loadit(generate.getfilenamegridlocal("data", "5_19_26_c", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]
skip = isnothing.(data)

unifprior = ones(Float64, nρ) ./ nρ

getprior = function(λ)
    expprior = λ .* exp.(-λ .* ρs)
    expprior ./= sum(expprior)
    expprior
end

plot(ρs, getprior(1), xlabel = "ρ", ylabel = "density",
    label = false, color = :black, grid = false)

jointprior = function(prior)
    [prior[i] * prior[j] for i in 1:nρ, j in 1:nρ]
end

getstats = function(i, j, prior)
    tmpdat = data[i, j]
    agg = log.(jointprior(prior))
    for k in 1:J
        agg .+= tmpdat[k, :, :]
    end
    agg .-= sum(agg)
    agg
end

plotheat = function(i, j, prior)
    heatmap(ρs, ρs, getstats(i, j, prior)',
        xlabel = "ρ0", ylabel = "ρ1",
        xlim = [0, maxρ], ylim = [0, maxρ],
        c = cgrad(:diff, rev = true), 
        title = "ρ0: $ρ0, ρ1: $ρ1",
        aspect_ratio = 1, colorbar_title = "log P")
end

plotp0 = [sum(getstats(i, 1, unifprior)[:, 1]) for i in 1:nρ]
# plotp0 = [sum(getstats(i, 1, getprior(1))[:, 1]) for i in 1:nρ]
plot(ρs, plotp0,
    xlabel = "ρ0", ylabel = "P(ρ1 = 0)",
    label = false, grid = false, color = :blacl)

plotheat(1, 1, unifprior)
# plotheat(1, 1, getprior(1))
