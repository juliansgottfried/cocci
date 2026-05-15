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
@everywhere maxρ = 10 - dρ
@everywhere nρ = floor(Int, maxρ / dρ) + 1

@everywhere dt = 0.01
@everywhere maxtime = 1
@everywhere change = 30
@everywhere covariate = generate.buildcov(dt, maxtime, change)

pmap(1:2nρ) do i
	ρ = (0:dρ:maxρ)[mod(i - 1, nρ) + 1]
	isρ0 = i <= nρ
	filename = generate.getfilename("prob", "5_15_26_b", isρ0, ρ)
	if !isfile(filename)
		println("ρ: $ρ, ρ0: $isρ0")
		results = estimate.montecarlo(n, l1, m, isρ0 * ρ, !isρ0 * ρ, covariate)
		save_object(filename, results)
	end
end

#= pmap((nρ + 1):2nρ) do i
	ρ = (0:dρ:maxρ)[mod(i - 1, nρ) + 1]
	isρ0 = i <= nρ
	filename = generate.getfilename("prob", "5_15_26", isρ0, ρ)
	if !isfile(filename)
		println("ρ: $ρ, ρ0: $isρ0")
		results = estimate.montecarlo(n, l1, m, isρ0 * ρ, !isρ0 * ρ, covariate)
		save_object(filename, results)
	end
end =#