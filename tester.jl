using Plots
using JLD2

include("estimate.jl")
include("generate.jl")

ρs = 0:0.1:15
filenames = string.("./results/results_", replace.(string.(ρs), "." => "_"), ".jld2")
collect = []
for i in filenames
    push!(collect, load_object(i)[1])
end

# from cluster run
n = 50
l1 = 100
m = 100000

minval = Inf
for i in eachindex(ρs)
    collect[i]
    for j in 1:(n - 1), k in 1:j
        tmpvec = collect[i][j, k][2]
        tmpmin = minimum(tmpvec[tmpvec .> 0])
        if tmpmin < minval
            minval = tmpmin
        end
    end
end
minval /= 10

# for ground truth
ρ = 5
θ = 5
configs, dists = generate.getconfigs(n, l1, ρ, θ)

pseudo = minval
loglik = estimate.getl(ρs, n, collect, configs, dists, pseudo)

plot(ρs, loglik, xlabel = "ρ", ylabel = "loglik",
    label = false, color = :black)
vline!([ρ], label = false)

ρhat = ρs[argmax(loglik)]
