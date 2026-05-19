#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2, DelimitedFiles

@everywhere n = 17
@everywhere l1 = 25
@everywhere m = 100000

@everywhere dρ = 1
@everywhere maxρ = 10 - dρ
@everywhere ρs = 0:dρ:maxρ
@everywhere nρ = length(ρs)

@everywhere covariate = readdlm("/scratch/users/jgottf/cocci/rodent_data/covariate.csv", ',', Any, '\n')

#= @everywhere dt = 0.01
@everywhere maxtime = 1
@everywhere change = 30
@everywhere covariate = generate.buildcov(dt, maxtime, change) =#

#= pmap(1:2nρ) do i
	ρ = ρs[mod(i - 1, nρ) + 1]
	isρ0 = i <= nρ
	filename = generate.getfilename("prob", "5_18_26_b", isρ0, ρ)
	if !isfile(filename)
		results = estimate.montecarlo(n, l1, m, isρ0 * ρ, !isρ0 * ρ, covariate)
		save_object(filename, results)
	end
end =#

#= pmap((nρ + 1):2nρ) do i
	ρ = ρs[mod(i - 1, nρ) + 1]
	isρ0 = i <= nρ
	filename = generate.getfilename("prob", "5_17_26_d", isρ0, ρ)
	if !isfile(filename)
		results = estimate.montecarlo(n, l1, m, isρ0 * ρ, !isρ0 * ρ, covariate)
		save_object(filename, results)
	end
end =#

pmap(1:nρ^2) do i
	ρ0 = ρs[Int(div(i - 1, nρ)) + 1]
	ρ1 = ρs[Int(mod(i - 1, nρ)) + 1]
	filename = generate.getfilenamegrid("prob", "5_18_26_e", ρ0, ρ1)
	if !isfile(filename)
		results = estimate.montecarlo(n, l1, m, ρ0, ρ1, covariate)
		save_object(filename, results)
	end
end
