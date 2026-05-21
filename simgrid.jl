using Plots, JLD2, StatsBase

include("estimate.jl")
include("generate.jl")

logsumexp = function(x)
    xmax = maximum(x)
    xmax + log(sum(exp.(x .- xmax)))
end

plotheat = function(input, lower, upper, title, colorbar_title)
    heatmap(ρs, ρs, input',
            xlabel = "ρ0", ylabel = "ρ1",
            xlim = [0, maxρ], ylim = [0, maxρ],
            clim = (lower, upper),
            c = cgrad(:diff, rev = false), 
            title = title,
            aspect_ratio = 1, colorbar_title = colorbar_title)
end

dρ = 0.5
maxρ = 20 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ

J = 500

gatherdata = [load_object(generate.getfilenamegridlocal("data", "5_20_26_a", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]

store = zeros(Float64, nρ, nρ, 3, 2)
qs = [0.055, 0.5, 0.945]
for i in 1:nρ, j in 1:nρ store[i, j, :, :] = reduce(hcat, [quantile(gatherdata[i, j][:, k], qs) for k in 1:2]) end

estimates = zeros(Float64, nρ, nρ, 2)
for i in 1:nρ, j in 1:nρ estimates[i, j, :] = store[i, j, 2, :] end

# plot example

plotheat(estimates[:, :, 1], 0, maxρ, "estimating ρ0", "")
plotheat(estimates[:, :, 2], 0, maxρ, "estimating ρ1", "")



lgetprior = function(λ)
    prior = log(λ) .- λ .* ρs
    prior .- logsumexp(prior)
end

getcounts = function(input)
    counts = zeros(Float64, nρ, nρ)
    idx = Int.(div.(input[:, 1:2], dρ)) .+ 1
    for j in 1:J
        counts[idx[j, 1], idx[j, 2]] += 1
    end
    counts
end

# bin by 2
getprobs = function(i, j)
    counts = getcounts(gatherdata[i, j])
    counts[counts .== 0] .= 0.01
    for i in 2:(nρ - 1), j in 2:(nρ - 1)
        counts[i, j] = sum(counts[(i - 1):(i + 1), (j - 1):(j + 1)]) / 9
    end
    counts ./ sum(counts)
end

allprobs = [getprobs(i, j) for i in 1:nρ, j in 1:nρ]

λ = 0.001
prior = lgetprior(λ)

posterior = function(input)
    terms = [prior[i] + prior[j] + sum(log.(allprobs[i, j]) .* getcounts(input)) for i in 1:nρ, j in 1:nρ]
    terms .- logsumexp(terms)
end

posts = [posterior(gatherdata[i, j]) for i in 1:nρ, j in 1:nρ]

pick0 = 1
pick1 = 1
plotheat(posts[pick0, pick1], -Inf, 0, "ρ0: $(ρs[pick0]), ρ0: $(ρs[pick1])", "log P")

p0 = [sum(posts[i, j][:, 1]) for i in 1:nρ, j in 1:nρ]
plotheat(p0, -Inf, 0, "", "P(0)")

alpha = 0.05
issig = exp.(p0s) .< alpha
plotheat(issig, 0, 1, "", "significant")

# real
exp(sum(posterior(ρhat)[:, 1])) < alpha