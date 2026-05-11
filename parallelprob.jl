#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2

@everywhere n = 17
@everywhere l1 = 50
@everywhere m = 100000

@everywhere dρ = 0.5
@everywhere maxρ = 10
@everywhere nρ = Int(div(maxρ, dρ)) + 1

@everywhere dt = 0.01
@everywhere maxtime = 3
@everywhere change = 15
@everywhere covariate = generate.buildcov(dt, maxtime, change)

pmap(1:nρ^2) do i
	idx1 = div(i, nρ + 1) + 1
	idx2 = i - nρ * (idx1 - 1)
	ρ0 = dρ * (idx1 - 1)
	ρ1 = dρ * (idx2 - 1)
	results = estimate.montecarlo(n, l1, ρ0, ρ1, covariate, dt, m)
	save_object(generate.getfilename("prob", "5_11_26", ρ0, ρ1), results)
end
