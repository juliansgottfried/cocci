#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2

@everywhere n = 17
@everywhere l1 = 50

@everywhere dρ = 0.5
@everywhere maxρ = 10
@everywhere nρ = Int(div(maxρ, dρ)) + 1

@everywhere collect = [generate.getobj("prob", "5_11_26", dρ, idx1, idx2) for idx1 in 1:nρ, idx2 in 1:nρ]

@everywhere pseudo = estimate.getpseudo(collect, n)

@everywhere dt = 0.01
@everywhere maxtime = 3
@everywhere change = 15
@everywhere covariate = generate.buildcov(dt, maxtime, change)

@everywhere J = 1000
@everywhere θ = 10

pmap(1:nρ^2) do i
	idx1 = div(i - 1, nρ) + 1
	idx2 = i - nρ * (idx1 - 1)
	ρ0 = dρ * (idx1 - 1)
	ρ1 = dρ * (idx2 - 1)
	# θ = iszero(ρ0 + ρ1) ? 0.1 : ρ0 + ρ1
	ρhat = generate.repeated(collect, dρ, nρ, pseudo, n, l1, ρ0, ρ1, covariate, dt, θ, J)
	save_object(generate.getfilename("data", "5_11_26", ρ0, ρ1), ρhat)
end
