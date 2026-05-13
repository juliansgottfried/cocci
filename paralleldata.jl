#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2

@everywhere n = 17
@everywhere l1 = 25

@everywhere dρ = 0.1
@everywhere maxρ = 10
@everywhere nρ = Int(div(maxρ, dρ)) + 1

@everywhere dt = 0.01
@everywhere maxtime = 3
@everywhere change = 20
@everywhere covariate = generate.buildcov(dt, maxtime, change)

@everywhere J = 1000
@everywhere θ = 10

@everywhere collect0 = [generate.getfilename("prob", "5_13_26", true, ρ) for ρ in 0:dρ:maxρ]
@everywhere collect1 = [generate.getfilename("prob", "5_13_26", false, ρ) for ρ in 0:dρ:maxρ]

@everywhere pseudo0 = estimate.getpseudo(collect0, n)
@everywhere pseudo1 = estimate.getpseudo(collect1, n)

pmap(1:2nρ) do i
	# idx1 = div(i - 1, nρ) + 1
	# idx2 = i - nρ * (idx1 - 1)
	
	# ρ0 = dρ * (idx1 - 1)
	# ρ1 = dρ * (idx2 - 1)

	θ = iszero(ρ0 + ρ1) ? 0.1 : ρ0 + ρ1
	ρ = (0:dρ:maxρ)[mod(i - 1, nρ) + 1]
	isρ0 = i <= nρ
	filename = generate.getfilename("data", "5_11_26", ρ0, ρ1)
	if !isfile(filename)
		ρhat = generate.repeated(collect, dρ, nρ, pseudo, n, l1, ρ0, ρ1, covariate, dt, θ, J)
		save_object(filename, ρhat)
	end
end
