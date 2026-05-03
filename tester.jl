using Plots
using JLD2

include("estimate.jl")
include("generate.jl")

n = 50
l1 = 50
m = 100000
ρs = 0:0.1:10

results = estimate.rhogrid(n, l1, ρs, m)
# load_object("results/results_1.jld2")

ρ = 2
θ = 100
configs, dists = generate.getconfigs(n, l1, ρ, θ)

pseudo = 0.1
loglik = estimate.getl(ρs, n, results, configs, dists, pseudo)

plot(ρs, loglik, xlabel = "ρ", ylabel = "loglik",
    label = false, color = :black)
vline!([ρ], label = false)