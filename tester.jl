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

# for ground truth
ρ = 10
θ = 100
configs, dists = generate.getconfigs(n, l1, ρ, θ)

pseudo = 0.001
loglik = estimate.getl(ρs, n, collect, configs, dists, pseudo)

plot(ρs, loglik, xlabel = "ρ", ylabel = "loglik",
    label = false, color = :black)
vline!([ρ], label = false)