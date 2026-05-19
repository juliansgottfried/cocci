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

n = 17

gather = [load_object(generate.getfilenamegridlocal("prob", "5_19_26_a", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]

pseudo = estimate.getpseudo(gather, n)

covariate = readdlm("rodent_data/covariate.csv", ',', Any, '\n')
alleles = readdlm("sampling_data/alleles.csv", ',', Any, '\n')

chunksize = Int(floor(maximum(alleles[:, end]) / 20))
chunkidx = 2580020
chunkloci = (alleles[:, end] .>= chunkidx) .& (alleles[:, end] .<= chunkidx + chunksize)
chunk = alleles[chunkloci, :]
chunk[:, end] .-= minimum(chunk[:, end]) - 1

nloci = size(chunk)[1]
nsample = 50
window = 400000
nwindow = chunksize - window

# G = 365
# η = 0.95
# maxbrk = 50
# pvec = kernels.getkernels(G, η, maxbrk)
pvec = 1

S = 100
ρhat = zeros(Float64, S, 3)
for i in 1:S
    subset = falses(nloci)
    println(i)
    while sum(subset) < nsample
        windowidx = StatsBase.sample((1:nwindow), 1)[1]
        subset = (chunk[:, end] .>= windowidx) .& (chunk[:, end] .<= windowidx + window)
    end
    sampleidx = StatsBase.sample((1:nloci)[subset], nsample, replace = false)
    samples = chunk[sampleidx, :]

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


# p-value

fauxdρ = 0.1
fauxmaxρ = 20 - fauxdρ
fauxnρ = length(0:fauxdρ:fauxmaxρ)
fauxρs = 0:fauxdρ:fauxmaxρ
J = 500
faux = [load_object(generate.getfilenamelocal("data", "5_18_26_d", false, ρ)) for ρ in fauxρs]

getp = function(input, ρs, J)
    counts = round.(sort(input), sigdigits = 3)
    [sum(counts .== ρ) / J for ρ in ρs]
end

counts0 = faux[1][:, 3]
p0 = getp(counts0, fauxρs, J)
p0smooth = copy(p0)
for i in 2:(fauxnρ - 1)
    p0smooth[i] = sum(p0[(i - 1):(i + 1)]) / 3
end
p0smooth[p0smooth .== 0] .= minimum(p0smooth[p0smooth .> 0]) / 10
p0smooth ./= sum(p0smooth)

pc0 = zeros(fauxnρ)
for i in 2:fauxnρ
    countsc0 = data1[i][:, 3]
    pc0 .+= getp(countsc0, fauxρs, J) ./ (fauxnρ - 1)
end
pc0[pc0 .== 0] .= minimum(pc0[pc0 .> 0]) / 10
pc0 ./= sum(pc0)

ll0 = sum(ρhat[:, 7]) / J
llc0 = sum(ρhat[:, 8]) / J

obs = round.(ρhat[:, 3], sigdigits = 3)
obs = obs[obs .<= simmaxρ]
cutJ = length(obs)
table = Int.(round.(cutJ * getp(obs, fauxρs, cutJ), sigdigits = 3))

distr0 = Distributions.Multinomial(cutJ, p0smooth)
distrc0 = Distributions.Multinomial(cutJ, pc0)

lpdf0 = Distributions.logpdf(distr0, table)
lpdfc0 = Distributions.logpdf(distrc0, table)

x = lpdfc0 + llc0 - lpdf0 - ll0
pvalue = 1 / exp(x)
