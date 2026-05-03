#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("estimate.jl")
@everywhere include("generate.jl")
@everywhere using JLD2

@everywhere n = 25
@everywhere l1 = 50
@everywhere m = 1000
@everywhere ρs = 0:0.1:10

pmap(ρs) do ρ
	results = estimate.montecarlo(n, l1, ρ, m)
    save_object(string("/scratch/users/jgottf/cocci/results/results_",replace(string(ρ), "." => "_"), ".jld2"), [results, (n, l1, m, ρ)])
end
