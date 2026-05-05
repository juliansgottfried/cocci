#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2

@everywhere n = 50
@everywhere l1 = 100
@everywhere m = 100000
@everywhere ρs = 6:0.1:15

pmap(ρs) do ρ
	results = estimate.montecarlo(n, l1, ρ, m)
	save_object(string("/scratch/users/jgottf/cocci/results/prob/run_5_3_26/results_", replace(string(ρ), "." => "_"), ".jld2"), [results, (n, l1, m, ρ)])
end
