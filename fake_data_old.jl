using Distributions
using StatsBase

N = 1000
L = 10000
G = 10000

r = 5 / N
μ = 0.001 / N

mutations = Binomial(L, μ)

genomes = falses(L, N)
prev = similar(genomes)

for i in 1:G
    prev .= genomes
    for j in 1:N
        if rand() < r
            # println("break")
            parents = StatsBase.sample(1:N, 2)
            junction = ceil(Int, rand() * L)
            child = [prev[1:(junction - 1), parents[1]] ; 
                prev[junction:L, parents[2]]]
        else 
            child = prev[:, StatsBase.sample(1:N)]
        end
        
        mutate = StatsBase.sample(1:L, rand(mutations), replace = false)
        child[mutate] .|= trues(size(mutate)[1])
        
        genomes[:, j] = child
    end
end

getconfig = function(input)
    config = sum(eachcol(input))
    push!(config, sum(input[1, :] .& input[2, :]))
end

data = genomes[:, StatsBase.sample(1:N, n, replace = false)]
loci = (1:L)[0 .< sum(eachcol(data)) .< n]
subsetdata = data[loci, :]

k = size(loci)[1]
configs = zeros(Int, Int(k * (k - 1) / 2), 3)
dists = zeros(Float64, size(configs)[1])
for i in 1:k
    for j in 1:(i - 1)
        c = Int((i * (i - 3)) / 2 + 1 + j)
        configs[c, :] = getconfig(subsetdata[[i; j], :])
        dists[c] = abs(loci[i] - loci[j])/L
    end
end

configs
dists