#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2

@everywhere n = 17
@everywhere l1 = 25
@everywhere θ = 10
@everywhere J = 1000

@everywhere dρ = 0.1
@everywhere maxρ = 10
@everywhere nρ = Int(div(maxρ, dρ)) + 1

@everywhere dt = 0.01
@everywhere maxtime = 3
@everywhere change = 20
@everywhere covariate = generate.buildcov(dt, maxtime, change)

@everywhere collect0 = [load_object(generate.getfilename("prob", "5_13_26", true, ρ)) for ρ in 0:dρ:maxρ]
@everywhere collect1 = [load_object(generate.getfilename("prob", "5_13_26", false, ρ)) for ρ in 0:dρ:maxρ]

@everywhere pseudo0 = estimate.getpseudo(collect0, n)
@everywhere pseudo1 = estimate.getpseudo(collect1, n)

pmap(1:2nρ) do i
	# idx1 = div(i - 1, nρ) + 1
	# idx2 = i - nρ * (idx1 - 1)
	# ρ0 = dρ * (idx1 - 1)
	# ρ1 = dρ * (idx2 - 1)
	ρ = (0:dρ:maxρ)[mod(i - 1, nρ) + 1]
	isρ0 = i <= nρ
	# θ = iszero(ρ) ? 0.1 : ρ
	filename = generate.getfilename("data", "5_13_26", isρ0, ρ)
	if !isfile(filename)
		ρhat = generate.repeated(collect0, collect1, pseudo0, pseudo1,
			n, l1, θ, J, dρ, nρ,
			isρ0 * ρ, !isρ0 * ρ, covariate, dt)
		save_object(filename, ρhat)
	end
end
