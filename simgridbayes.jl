using Plots, JLD2, StatsBase

include("estimate.jl")
include("generate.jl")

logsumexp = function(x)
    xmax = maximum(x)
    xmax + log(sum(exp.(x .- xmax)))
end

dρ = 0.5
maxρ = 20 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ
J = 500

loadit(filename) = if isfile(filename) return load_object(filename) end
data = [loadit(generate.getfilenamegridlocal("data", "5_19_26_c", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]
skip = isnothing.(data)

lunifprior = log.(ones(Float64, nρ) ./ nρ)

lgetprior = function(λ)
    expprior = log(λ) .- λ .* ρs
    expprior .- logsumexp(expprior)
end

# plot(ρs, lgetprior(10), xlabel = "ρ", ylabel = "density",
#     label = false, color = :black, grid = false)

jointprior = function(prior0, prior1)
    [prior0[i] + prior1[j] for i in 1:nρ, j in 1:nρ]
end

getstats = function(i, j, ndata, real)
    if real tmpdat = ρhat
    else tmpdat = data[i, j]
    end
    agg = jointprior(prior0, prior1)
    for k in 1:ndata
        agg .+= tmpdat[k, :, :]
    end
    agg .- logsumexp(agg)
end

plotheat = function(i, j, cutoff, ndata, real)
    ρ0 = ρs[i]
    ρ1 = ρs[j]
    vals = getstats(i, j, ndata, real)
    vals[vals .< cutoff] .= cutoff
    heatmap(ρs, ρs, vals',
        xlabel = "ρ0", ylabel = "ρ1",
        xlim = [0, maxρ], ylim = [0, maxρ],
        c = cgrad(:default, rev = false), 
        title = "ρ0: $ρ0, ρ1: $ρ1",
        aspect_ratio = 1, colorbar_title = "log P")
end

plotheatall = function(idx, ndata)
    best = zeros(Float64, nρ, nρ)
    for i in 1:nρ, j in 1:nρ
        if skip[i, j] best[i, j] = Inf 
        else
            vals = getstats(i, j, ndata, false)
            best[i, j] = ρs[argmax(vals)[idx]]
        end
    end
    heatmap(ρs, ρs, best',
        xlabel = "ρ0", ylabel = "ρ1",
        xlim = [0, maxρ], ylim = [0, maxρ],
        c = cgrad(:default, rev = false),
        clim = (0, 20), grid = false,
        aspect_ratio = 1, colorbar_title = "hat")
end

plotthresh = function(cutoff, use, ndata)
    p0 = zeros(Float64, nρ, nρ)
    for i in 1:nρ, j in 1:nρ
        if skip[i, j] p0[i, j] = Inf 
        else
            vals = getstats(i, j, ndata, false)
            p0[i, j] = logsumexp(vals[:, 1])
        end
    end
    if use p0 = p0 .< cutoff end
    heatmap(ρs, ρs, p0',
        xlabel = "ρ0", ylabel = "ρ1",
        xlim = [0, maxρ], ylim = [0, maxρ],
        c = cgrad(:default, rev = true),
        aspect_ratio = 1, colorbar_title = "log p0")
end

plotbic = function(ndata)
    k = 1
    bic = zeros(Float64, nρ, nρ)
    for i in 1:nρ, j in 1:nρ
        if skip[i, j] bic[i, j] = Inf 
        else
            vals = getstats(i, j, ndata, false)
            bic0 = k * log(ndata) - 2maximum(vals[:, 1])
            bic1 = k * log(ndata) - 2maximum(vals[1, :])
            bic[i, j] = bic1 < bic0
        end
    end
    heatmap(ρs, ρs, bic',
        xlabel = "ρ0", ylabel = "ρ1",
        xlim = [0, maxρ], ylim = [0, maxρ],
        c = cgrad(:default, rev = true),
        title = "support for varying model",
        aspect_ratio = 1, colorbar_title = "")
end

# plots

ndata = 500
λ0 = 1
λ1 = 10
thresh = -900
prior0 = lunifprior
prior1 = lgetprior(λ1)

plotbic(ndata)

plotthresh(thresh, true, ndata)
plotthresh(thresh, false, ndata)

# plotheat(2, 15, thresh, ndata, false)
# plotheat(9, 1, thresh, ndata, false)
plotheatall(1, ndata)
plotheatall(2, ndata)


# data
iseffect = logsumexp(getstats(1, 1, ndata, true)[:, 1]) > thresh
plotheat(1, 1, thresh, ndata, true)
