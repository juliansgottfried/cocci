#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2

@everywhere n = 50
@everywhere l1 = 100
@everywhere m = 100000
@everywhere ρs = 0:0.1:15

@everywhere filenames = string.("/scratch/users/jgottf/cocci/results/prob/run_5_3_26/results_", replace.(string.(ρs), "." => "_"), ".jld2")
@everywhere collect = []
@everywhere for i in filenames
	push!(collect, load_object(i)[1])
end

@everywhere pseudo = Inf
@everywhere for i in eachindex(ρs)
	collect[i]
    for j in 1:(n - 1), k in 1:j
        tmpvec = collect[i][j, k][2]
        tmpmin = minimum(tmpvec[tmpvec .> 0])
        if tmpmin < pseudo
            pseudo = tmpmin
        end
    end
end
@everywhere pseudo /= 10

@everywhere J = 100
@everywhere θ = 10

pmap(ρs) do ρ
	ρhat = generate.repeated(ρs, collect, pseudo, n, l1, ρ, θ, J)
	save_object(string("/scratch/users/jgottf/cocci/results/data/run_5_4_26/results_", replace(string(ρ), "." => "_"), ".jld2"), [ρhat, (ρ, θ, J)])
end
