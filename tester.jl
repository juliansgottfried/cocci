using Plots
import StatsPlots
using JLD2

include("estimate.jl")
include("generate.jl")

ρs = 0:0.1:15
filenames = string.("./results/prob/run_5_3_26/results_", replace.(string.(ρs), "." => "_"), ".jld2")
collect = []
for i in filenames
    push!(collect, load_object(i)[1])
end

# from cluster run
n = 50
l1 = 100
m = 100000

pseudo = Inf
for i in eachindex(ρs)
    collect[i]
    for j in 1:(n - 1), k in 1:j
        tmpvec = collect[i][j, k][2]
        tmpmin = minimum(tmpvec[tmpvec .> 0])
        if tmpmin < pseudo
            pseudo = tmpmin
        end
    end
end
pseudo /= 10

# for ground truth
ρ = 5
θ = 10

# single run
configs, dists = generate.getconfigs(n, l1, ρ, θ)
loglik = estimate.getl(ρs, n, collect, configs, dists, pseudo)
plot(ρs, loglik, xlabel = "ρ", ylabel = "loglik",
    label = false, color = :black)
vline!([ρ], label = false, color = :red)

# many runs
J = 100
ρhat = generate.repeated(ρs, collect, pseudo, n, l1, ρ, θ, J)
plot(1:J, ρhat/ρ,  xlabel = "sample", ylabel = "ρ_hat / ρ",
    label = false, color = :black)
hline!([1], label = false, color = :red)

StatsPlots.density(ρhat/ρ, label = false, color = :black)