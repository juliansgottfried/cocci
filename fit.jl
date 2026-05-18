using DelimitedFiles, Plots
import StatsBase, StatsPlots, Distributions

include("estimate.jl")
include("generate.jl")

logsumexp = function(x)
    xmax = maximum(x)
    xmax + log(sum(exp.(x .- xmax)))
end

dρ = 0.1
maxρ = 100 - dρ

n = 17

collect0 = [load_object(generate.getfilenamelocal("prob", "5_18_26_b", true, ρ)) for ρ in 0:dρ:maxρ]
collect1 = [load_object(generate.getfilenamelocal("prob", "5_18_26_b", false, ρ)) for ρ in 0:dρ:maxρ]

pseudo0 = estimate.getpseudo(collect0, n)
pseudo1 = estimate.getpseudo(collect1, n)
pseudo = min(pseudo0, pseudo1)

covariate = readdlm("rodent_data/covariate.csv", ',', Any, '\n')
alleles = readdlm("sampling_data/alleles.csv", ',', Any, '\n')

nloci = size(alleles)[1]
nsample = 25
window = 5000
nwindow = maximum(alleles[:, end]) - window

S = 500
ρhat = zeros(Float64, S, 8)
for i in 1:S
    subset = falses(nloci)
    println(i)
    while sum(subset) < nsample
        windowidx = StatsBase.sample((1:nwindow), 1)[1]
        subset = (alleles[:, end] .>= windowidx) .& (alleles[:, end] .<= windowidx + window)
    end
    sampleidx = StatsBase.sample((1:nloci)[subset], nsample, replace = false)
    samples = alleles[sampleidx, :]

    configs = zeros(Int, Int(nsample * (nsample - 1) / 2), 3)
    dists = zeros(Float64, size(configs)[1])
    for i in 2:nsample
        for j in 1:(i - 1)
            c = Int((i * (i - 3)) / 2 + 1 + j)
            configs[c, :] = generate.getconfig(samples[[i; j], 1:(end - 1)])
            dists[c] = abs(samples[i, end] - samples[j, end]) / window
        end
    end
    loglik0, loglik1 = estimate.getl(n, collect0, collect1, pseudo, pseudo,
            dρ, maxρ, configs, dists)

    ρhat[i, 5:6] = [loglik0[1]; logsumexp(loglik0[2:end])]
    ρhat[i, 7:8] = [loglik1[1]; logsumexp(loglik1[2:end])]
    
    idx0 = argmax(loglik0)
    idx1 = argmax(loglik1)
    lik0 = maximum(loglik0)
    lik1 = maximum(loglik1)
    bestρ0 = dρ * (idx0 - 1)
    bestρ1 = dρ * (idx1 - 1)
    ρhat[i, 1:4] = [bestρ0; lik0; bestρ1; lik1]
end

aic = -2(ρhat[:, 2] .- ρhat[:, 4])
StatsPlots.density(aic, color = :black, 
    xlabel = "Δ AIC", ylabel = "density",
    title = "constant - varying",
    label = false, grid = false)
vline!([2], c = :red, alpha = 0.7, label = false)

rate1 = sum((aic .> 0)) / J
diff = quantile(aic, 0.5)
pass = diff > 2

nbins = Int(floor((maxρ + 1) / 1))
histogram(ρhat[:, 1], linecolor = :white, color = :black, 
    xlabel = "recombination rate", ylabel = "count",
    bins = nbins, label = false, grid = false)
histogram(ρhat[:, 3], linecolor = :white, color = :black,
    xlabel = "recombination rate", ylabel = "count",
    bins = nbins, label = false, grid = false)


# p-value

simdρ = 0.1
simmaxρ = 20 - simdρ
simnρ = length(0:simdρ:simmaxρ)
ρs = 0:simdρ:simmaxρ
J = 500
data1 = [load_object(generate.getfilenamelocal("data", "5_17_26_d", false, ρ)) for ρ in ρs]

getp = function(input, ρs, J)
    counts = round.(sort(input), sigdigits = 3)
    p = [sum(counts .== ρ) / J for ρ in ρs]
end

counts0 = data1[1][:, 3]
p0 = getp(counts0, ρs, J)
# p0[p0 .== 0] .= minimum(p0[p0 .> 0]) / 10
# p0 ./= sum(p0)

p0smooth = copy(p0)
for i in 2:(nρ - 1)
    p0smooth[i] = sum(p0[(i - 1):(i + 1)]) / 3
end
p0smooth[p0smooth .== 0] .= minimum(p0smooth[p0smooth .> 0]) / 10
p0smooth ./= sum(p0smooth)

pc0 = zeros(nρ)
for i in 2:nρ
    countsc0 = data1[i][:, 3]
    pc0 .+= getp(countsc0, ρs, J) ./ (nρ - 1)
end
pc0[pc0 .== 0] .= minimum(pc0[pc0 .> 0]) / 10
pc0 ./= sum(pc0)

ll0 = StatsBase.mean(ρhat[:, 7])
llc0 = StatsBase.mean(ρhat[:, 8])

obs = round.(ρhat[:, 3], sigdigits = 3)
obs = obs[obs .<= simmaxρ]
cutJ = length(obs)
table = Int.(round.(cutJ * getp(obs, ρs, cutJ), sigdigits = 3))

distr0 = Distributions.Multinomial(cutJ, p0smooth)
distrc0 = Distributions.Multinomial(cutJ, pc0)

lpdf0 = Distributions.logpdf(distr0, table)
lpdfc0 = Distributions.logpdf(distrc0, table)

lpvalue = lpdf0 + ll0 - (lpdf0 + ll0 + log(1 + exp(lpdfc0 + llc0 - lpdf0 - ll0)))
pvalue = exp(lpvalue)


# now - do the kernel fitting
# no p-value, bc can't fit all that many