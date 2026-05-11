#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2

@everywhere n = 25
@everywhere l1 = 50
@everywhere m = 100000
@everywhere maxρ = 30
@everywhere ρs = generate.makegrid(maxρ)

@everywhere dt = 0.01
@everywhere maxtime = 3
@everywhere change = 20
@everywhere covariate = generate.buildcov(dt, maxtime, change)

pmap(ρs) do ρ
	results = estimate.montecarlo(n, l1, ρ[1], ρ[2], covariate, dt, m)
	save_object(generate.getfilename("prob", "5_8_26", ρ), [results])
end
