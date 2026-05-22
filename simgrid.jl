using Plots, JLD2, StatsBase

include("estimate.jl")
include("generate.jl")

logsumexp = function(x)
    xmax = maximum(x)
    xmax + log(sum(exp.(x .- xmax)))
end

plotheat = function(input, lower, upper, title, colorbar_title, binary)
    if binary input[input .!= 0] .= maxρ end
    heatmap(ρs, ρs, input',
            xlabel = "ρ0", ylabel = "ρ1",
            xlim = [0, maxρ], ylim = [0, maxρ],
            clim = (lower, upper),
            color = cgrad(:thermal, rev = false), 
            title = title,
            aspect_ratio = 1, colorbar_title = colorbar_title)
end

plotq = function(xvar, idx)
    if xvar == 1 eg = store[:, idx, :, xvar]
    else eg = store[idx, :, :, xvar] end
    xlab = xvar == 1 ? "ρ0" : "ρ1"
    title = xvar == 1 ? "ρ1 = $(ρs[idx])" : "ρ0 = $(ρs[idx])"
    plot(ρs, eg[:, 1], fillrange = eg[:, 3], 
        xlabel = xlab, ylabel = "estimate",
        title = title,
        xlim = [0, maxρ], fillalpha = 0.15, color = "#29a0c8",
        linecolor = false, label = false, grid = false)
    plot!(ρs, maxρ * ones(nρ),linestyle = :dash, color = "#29a0c8", label = false)
    plot!(ρs, zeros(nρ), linestyle = :dash, color = "#29a0c8", label = false)
    plot!(ρs, eg[:, 2], color = :black, linewidth = 1.2, label = false)
    plot!([0, maxρ], [0, maxρ], color = :red, label = false)
end

dρ = 0.5
maxρ = 19.5 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ

J = 500

gatherdata = [load_object(generate.getfilenamegridlocal("data", "5_20_26_a", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]
for i in 1:nρ, j in 1:nρ
    tmpgather = gatherdata[i, j][:, 1:2]
    tmpgather[tmpgather .> maxρ] .= maxρ
    gatherdata[i, j][:, 1:2] = tmpgather
end

store = zeros(Float64, nρ, nρ, 3, 2)
qs = [0.055, 0.5, 0.945]
for i in 1:nρ, j in 1:nρ store[i, j, :, :] = reduce(hcat, [quantile(gatherdata[i, j][:, k], qs) for k in 1:2]) end

estimates = zeros(Float64, nρ, nρ, 2)
for i in 1:nρ, j in 1:nρ estimates[i, j, :] = store[i, j, 2, :] end

plotq(2, 1)
plotheat(estimates[:, :, 1], 0, maxρ, "estimating ρ0", "", false)
plotheat(estimates[:, :, 2], 0, maxρ, "estimating ρ1", "", false)

getprior = function(λ)
    prior = log(λ) .- λ .* ρs
    prior .- logsumexp(prior)
end

getcounts = function(input)
    counts = zeros(Float64, nρ, nρ)
    idx = Int.(div.(input[:, 1:2], dρ)) .+ 1
    for i in 1:J
        counts[idx[i, 1], idx[i, 2]] += 1
    end
    counts
end

guass = function(σ, i, j)
    (1 / (sqrt(2pi) * σ)) * exp(-(i^2 + j^2) / 2σ^2)
end

σ = 0.8
gausskern = [guass(σ, i, j) for i in -1:1, j in -1:1]

getprobs = function(i, j)
    counts = getcounts(gatherdata[i, j])
    counts[counts .== 0] .= 0.1
    smoothed = copy(counts)
    for k in 1:nρ, l in 1:nρ
        bottom = max(1, k - 1)
        top = min(nρ, k + 1)
        left = max(1, l - 1)
        right = min(nρ, l + 1)
        sect = counts[bottom:top, left:right]
        kern = gausskern[(bottom:top) .- k .+ 2, (left:right) .- l .+ 2]
        smoothed[k, l] = sum(sect .* kern)
    end
    log.(smoothed ./ sum(smoothed))
end

# plotheat(getprobs(1, 39), -Inf, Inf, "", "", false)

allprobs = [getprobs(i, j) for i in 1:nρ, j in 1:nρ]


λ0 = 30
λ1 = 2
prior0 = getprior(λ0)
prior1 = getprior(λ1)

posterior = function(input)
    tmpcounts = getcounts(input)
    terms = [prior0[i] + prior1[j] + sum(allprobs[i, j] .* tmpcounts) for i in 1:nρ, j in 1:nρ]
    # terms[terms .== 0] .= maximum(terms[terms .< 0]) + log(10)
    terms .-= logsumexp(terms)
    terms
end

posts = [posterior(gatherdata[i, j]) for i in 1:nρ, j in 1:nρ]

# pick0 = 10
# pick1 = 20
# plotheat(posts[pick0, pick1], -Inf, Inf, "ρ0: $(ρs[pick0]), ρ1: $(ρs[pick1])", "log P", false)

p0 = [logsumexp(posts[i, j][:, 1]) for i in 1:nρ, j in 1:nρ]
plotheat(p0, -Inf, Inf, "P(ρ1 = 0)", "log", false)

bottomn = 2
alpha = 0.01
cutoff = Int(floor(alpha * bottomn * nρ)) + 1
thresh = sort(reduce(vcat, p0[:, 1:bottomn]))[cutoff]
issig = p0 .< thresh
plotheat(issig, 0, 1, "significance", "", false)

# real
plotheat(posterior(ρhat), -Inf, Inf, "", "", false)
logsumexp(posterior(ρhat)[:, 1]) < thresh
