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

getcounts1d = function(input)
    counts = zeros(Float64, nρ, 2)
    idx = Int.(div.(input[:, 1:2], dρ)) .+ 1
    for j in 1:J
        counts[idx[j, 1], 1] += 1
        counts[idx[j, 2], 2] += 1
    end
    counts
end

guass1d = function(σ, i)
    (1 / (sqrt(2pi) * σ)) * exp(-i^2 / 2σ^2)
end

σ = 1
gausskern1d = guass1d.(σ, -1:1)

getprobs1d = function(i, j)
    counts = getcounts1d(gatherdata[i, j])
    counts[counts .== 0] .= 0.1
    smoothed = copy(counts)
    for k in 1:nρ
        lbound = max(1, k - 1)
        rbound = min(nρ, k + 1)
        sect = counts[lbound:rbound, :]
        kern = gausskern1d[(lbound:rbound) .- k .+ 2]
        smoothed[k, :] = sum(reduce(hcat, [sect[:, l] .* kern for l in 1:2]), dims = 1)
    end
    log.(reduce(hcat, [smoothed[:, l] ./ sum(smoothed[:, l]) for l in 1:2]))
end

allprobs = [getprobs1d(i, j) for i in 1:nρ, j in 1:nρ]

λ0 = 0.0001
λ1 = 0.0001
prior0 = getprior(λ0)
prior1 = getprior(λ1)

posterior1d = function(input)
    tmpcounts = getcounts1d(input)
    terms = [prior0[i] + prior1[j] + sum(allprobs[i, j] .* tmpcounts) for i in 1:nρ, j in 1:nρ]
    terms .-= logsumexp(terms)
    terms[terms .== 0] .= maximum(terms[terms .< 0]) + log(1000)
    terms
end

posts = [posterior1d(gatherdata[i, j]) for i in 1:nρ, j in 1:nρ]

pick0 = 1
pick1 = 20
# plotheat(posts[pick0, pick1], -1000, Inf, "ρ0: $(ρs[pick0]), ρ1: $(ρs[pick1])", "log P", false)

p0 = [logsumexp(posts[i, j][:, 1]) for i in 1:nρ, j in 1:nρ]
plotheat(p0, -Inf, Inf, "P(ρ1 = 0)", "log", false)

alpha = -20000
sum(p0[:, 3:end] .< alpha) / sum(p0 .< alpha)

issig = p0 .< alpha
plotheat(issig, 0, 1, "significance", "", false)

# real
exp(sum(posterior(ρhat)[:, 1])) < alpha