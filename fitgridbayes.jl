using DelimitedFiles, Plots
import StatsBase, StatsPlots, Distributions

include("estimate.jl")
include("generate.jl")
include("kernels.jl")

dρ = 0.5
maxρ = 20 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ

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
window = 1000
nwindow = chunksize - window

# G = 365
# η = 0.95
# maxbrk = 50
# pvec = kernels.getkernels(G, η, maxbrk)
pvec = 1

S = 100
ρhat = zeros(Float64, S, nρ, nρ)
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
    ρhat[i, :, :] = estimate.getlgrid(n, gather, pseudo,
            dρ, maxρ, configs, dists, pvec)
end
