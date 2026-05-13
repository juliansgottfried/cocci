#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2

@everywhere n = 17
@everywhere l1 = 25
@everywhere m = 100000

@everywhere dρ = 0.1
@everywhere maxρ = 10
@everywhere nρ = floor(Int, maxρ / dρ) + 1

@everywhere dt = 0.01
@everywhere maxtime = 3
@everywhere change = 20
@everywhere covariate = generate.buildcov(dt, maxtime, change)

pmap(1:2nρ) do i
	ρ = (0:dρ:maxρ)[mod(i - 1, nρ) + 1]
	isρ0 = i <= nρ
	if isρ0 results = estimate.montecarlo(n, l1, m, ρ, 0.0, covariate, dt)
	else results = estimate.montecarlo(n, l1, m, 0.0, ρ, covariate, dt) end
	save_object(generate.getfilename("prob", "5_13_26", isρ0, ρ), results)
end
