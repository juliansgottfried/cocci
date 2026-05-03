using Plots
using JLD2

include("estimate.jl")
include("generate.jl")

# insert code to accumulate results!
# load_object("results/results.jld2")

ρ = 2
θ = 100
configs, dists = generate.getconfigs(n, l1, ρ, θ)

pseudo = 0.1
loglik = estimate.getl(ρs, n, results, configs, dists, pseudo)

plot(ρs, loglik, xlabel = "ρ", ylabel = "loglik",
    label = false, color = :black)
vline!([ρ], label = false)