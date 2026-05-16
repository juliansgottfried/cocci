using DelimitedFiles, Plots
import StatsBase, StatsPlots

include("estimate.jl")
include("generate.jl")

dρ = 1
maxρ = 100 - dρ

n = 17

collect0 = [load_object(generate.getfilenamelocal("prob", "5_16_26_b", true, ρ)) for ρ in 0:dρ:maxρ]
collect1 = [load_object(generate.getfilenamelocal("prob", "5_16_26_b", false, ρ)) for ρ in 0:dρ:maxρ]

pseudo0 = estimate.getpseudo(collect0, n)
pseudo1 = estimate.getpseudo(collect1, n)

covariate = readdlm("rodent_data/covariate.csv", ',', Any, '\n')
alleles = readdlm("sampling_data/alleles.csv", ',', Any, '\n')

nloci = size(alleles)[1]
nsample = 25
window = 500
nwindow = maximum(alleles[:, end]) - window

S = 1000
ρhat = zeros(Float64, S, 4)
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
    loglik0, loglik1 = estimate.getl(n, collect0, collect1, pseudo0, pseudo1,
            dρ, maxρ, configs, dists)
    idx0 = argmax(loglik0)
    idx1 = argmax(loglik1)
    lik0 = maximum(loglik0)
    lik1 = maximum(loglik1)
    bestρ0 = dρ * (idx0 - 1)
    bestρ1 = dρ * (idx1 - 1)
    ρhat[i, :] = [bestρ0; lik0; bestρ1; lik1]
end

StatsPlots.density(2(ρhat[:, 4] .- ρhat[:, 2]), color = :black, label = false, grid = false)
vline!([2], c = :red, alpha = 0.7, label = false)

histogram(ρhat[:, 3], linecolor = :white, color = :black,
     bins = Int(floor(maxρ + 1) / 4), label = false, grid = false)

#= boxplot(ρhat[:, [1;3]],
    color = :white, whisker_width = 0.2, outliers = false,
    xticks = false, label = false, grid = false)
histogram(ρhat[:, 1], linecolor = :white, color = :black, 
    bins = Int(floor(maxρ + 1) / 4), label = false, grid = false) =#
    

# significance testing
nullρhat = generate.repeated(collect0, collect1, pseudo0, pseudo1,
    n, 17, 0.3, 1000, dρ, maxρ,
	0, 5, covariate, 1)

StatsBase.mean(2(nullρhat[:, 4] .- nullρhat[:, 2]))

StatsPlots.density(2(nullρhat[:, 4] .- nullρhat[:, 2]), color = :black, label = false, grid = false)
vline!([-2], c = :red, alpha = 0.7, label = false)

histogram(nullρhat[:, 3], linecolor = :white, color = :black,
     bins = Int(floor(maxρ + 1) / 4), label = false, grid = false)
