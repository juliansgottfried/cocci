using DelimitedFiles, Plots
import StatsBase, StatsPlots, Distributions

include("estimate.jl")
include("generate.jl")
include("kernels.jl")

logsumexp = function(x)
    xmax = maximum(x)
    xmax + log(sum(exp.(x .- xmax)))
end

dρ = 0.5
maxρ = 20 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ
J = 500

n = 17

gather = [load_object(generate.getfilenamegridlocal("prob", "5_19_26_a", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]

pseudo = estimate.getpseudo(gather, n)

covariate = readdlm("rodent_data/covariate.csv", ',', Any, '\n')
alleles = readdlm("sampling_data/alleles.csv", ',', Any, '\n')

L = maximum(alleles[:, end])
chunksize = Int(floor(L / 100))

bestchunk = nothing
mostloci = 0
for i in 1:(L - chunksize)
    if mod(i, 1000) != 0 continue end
    chunkloci = (alleles[:, end] .>= i) .& (alleles[:, end] .<= i + chunksize)
    chunk = alleles[chunkloci, :]
    nloci = size(chunk)[1]
    if nloci == 0 continue end
    chunk[:, end] .-= minimum(chunk[:, end]) - 1
    nloci = size(chunk)[1]
    if nloci > mostloci
        mostloci = nloci
        bestchunk = chunk
    end
end

nsample = 50
window = 10000
nwindow = chunksize - window

# G = 365
# η = 0.95
# maxbrk = 50
# pvec = kernels.getkernels(G, η, maxbrk)
pvec = 1

S = 500
ρhat = zeros(Float64, S, 3)
for i in 1:S
    subset = falses(mostloci)
    println(i)
    while sum(subset) < nsample
        windowidx = StatsBase.sample((1:nwindow), 1)[1]
        subset = (bestchunk[:, end] .>= windowidx) .& (bestchunk[:, end] .<= windowidx + window)
    end
    sampleidx = StatsBase.sample((1:mostloci)[subset], nsample, replace = false)
    samples = bestchunk[sampleidx, :]

    configs = zeros(Int, Int(nsample * (nsample - 1) / 2), 3)
    dists = zeros(Float64, size(configs)[1])
    for i in 2:nsample
        for j in 1:(i - 1)
            c = Int((i * (i - 3)) / 2 + 1 + j)
            configs[c, :] = generate.getconfig(samples[[i; j], 1:(end - 1)])
            dists[c] = abs(samples[i, end] - samples[j, end]) / window
        end
    end


    loglik = estimate.getlgrid(n, gather, pseudo,
        dρ, maxρ, configs, dists, pvec)
    idx = argmax(loglik)
    bestρ0 = dρ * (idx[1] - 1)
	bestρ1 = dρ * (idx[2] - 1)
    ρhat[i, :] = [bestρ0; bestρ1; maximum(loglik)]
end

nbins = Int.(floor.(nρ ./ 1))
histogram(ρhat[:, 1], linecolor = :white, color = :black, 
    xlabel = "recombination rate", ylabel = "density",
    xlim = [0, maxρ], normalize = true, bins = nbins,
    title = "ρ0 estimate",
    label = false, grid = false)
histogram(ρhat[:, 2], linecolor = :white, color = :black,
    xlabel = "recombination rate", ylabel = "density",
    xlim = [0, maxρ], normalize = true, bins = nbins,
    title = "ρ1 estimate",
    label = false, grid = false)



loadit(filename) = if isfile(filename) return load_object(filename) end
data = [loadit(generate.getfilenamegridlocal("data", "5_19_26_a", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]

getcounts = function(input)
    idx = Int.(div.(input[:, 1:2], dρ)) .+ 1
    [sum(idx[:, i] .== j) for j in 1:nρ, i in 1:2]
end

getprob = function(focal)
    focal = 20
    counts = zeros(Int, nρ, 2)
    for i in 1:nρ
        if skip[i, focal] continue end
        counts .+= getcounts(data[i, focal])
    end
    prob = counts ./ sum(counts, dims = 1)
    probsmooth = prob
    for i in 2:(nρ - 1)
        probsmooth[i, :] = sum(prob[(i - 1):(i + 1), :], dims = 1) ./ 3
    end
    # probsmooth[probsmooth .== 0] .= minimum(probsmooth[probsmooth .> 0]) / 10
    probsmooth[:, 1] ./= sum(probsmooth[:, 1])
    probsmooth[:, 2] ./= sum(probsmooth[:, 2])
    plot(ρs, probsmooth)
    [Distributions.Multinomial(S, probsmooth[:, 1]); Distributions.Multinomial(S, probsmooth[:, 2])]
end

distr = [getprob(i) for i in 1:nρ]

lgetprior = function(λ)
    expprior = log(λ) .- λ .* ρs
    expprior .- logsumexp(expprior)
end

λ = 0.001
prior = lgetprior(λ)

subsample = similar(data)
for i in 1:nρ, j in 1:nρ
    if isnothing(data[i, j])
        subsample[i, j] = nothing
        continue
    end
    subsample[i, j] = data[i, j][1:S, :]
end

p0 = function(input)
    if isnothing(input) return end
    obs = getcounts(input)
    terms = [prior[i] + Distributions.logpdf(distr[i][1], obs[:, 1]) + Distributions.logpdf(distr[i][2], obs[:, 2]) for i in 1:nρ]
    terms .-= logsumexp(terms)
    # plot(ρs, terms, grid = false, label = false, color = :black)
    terms[1]
end

p0s = [p0(subsample[i, j]) for i in 1:nρ, j in 1:nρ]
p0s[isnothing.(p0s)] .= Inf
alpha = 0.01
cut = exp.(p0s) .< alpha
heatmap(ρs, ρs, cut')

exp.(p0(ρhat)) < alpha




# remove multinomial because normalization constant cancels 

# P(rho1 = 0 | rho0_hat, rho1_hat) = 
# P(rho0_hat, rho1_hat | rho1 = 0)P(rho1 = 0) /
# Σ_r P(rho0_hat, rho1_hat | rho1 = r)P(rho1 = r)



loadit(filename) = if isfile(filename) return load_object(filename) end
data = [loadit(generate.getfilenamegridlocal("data", "5_19_26_a", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]

counts2d = function(input, N)
    counts = zeros(Float64, nρ, nρ)
    idx = Int.(div.(input[:, 1:2], dρ)) .+ 1
    for j in 1:N
        counts[idx[j, 1], idx[j, 2]] += 1
    end
    counts
end

getprobs = function(focal)
    counts = zeros(Float64, nρ, nρ)
    for i in 1:nρ
        if skip[i, focal] continue end
        counts .+= counts2d(data[i, focal], J)
    end
    counts[counts .== 0] .= minimum(counts[counts .> 0]) / 10
    counts ./= sum(counts)
end

lgetprior = function(λ)
    expprior = log(λ) .- λ .* ρs
    expprior .- logsumexp(expprior)
end

allprobs = [getprobs(i) for i in 1:nρ]

λ = 0.001
prior = lgetprior(λ)

subsample = similar(data)
for i in 1:nρ, j in 1:nρ
    if isnothing(data[i, j])
        subsample[i, j] = nothing
        continue
    end
    subsample[i, j] = data[i, j][1:S, :]
end

p0 = function(input, N)
    if isnothing(input) return end
    obs = counts2d(input, N)
    terms = [prior[i] + sum(log.(allprobs[i]) .* obs) for i in 1:nρ]
    terms .-= logsumexp(terms)
    terms[1]
end

p0s = [p0(subsample[i, j], S) for i in 1:nρ, j in 1:nρ]
p0s[isnothing.(p0s)] .= Inf
alpha = 0.001
cut = exp.(p0s) .< alpha
heatmap(ρs, ρs, cut')
heatmap(ρs, ρs, p0s')

exp(p0(ρhat, J)) < alpha
