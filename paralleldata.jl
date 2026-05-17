#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2, DelimitedFiles

@everywhere n = 17
@everywhere l1 = 25
@everywhere θ = 10
@everywhere J = 500

@everywhere dρ = 0.1
@everywhere maxρ = 10 - dρ
@everywhere nρ = length(0:dρ:maxρ)

@everywhere covariate = readdlm("/scratch/users/jgottf/cocci/rodent_data/covariate.csv", ',', Any, '\n')

#= @everywhere dt = 0.01
@everywhere maxtime = 1
@everywhere change = 30
@everywhere covariate = generate.buildcov(dt, maxtime, change) =#

@everywhere pvec = 1

@everywhere collect0 = [load_object(generate.getfilename("prob", "5_16_26_c", true, ρ)) for ρ in 0:dρ:maxρ]
@everywhere collect1 = [load_object(generate.getfilename("prob", "5_16_26_c", false, ρ)) for ρ in 0:dρ:maxρ]

@everywhere pseudo0 = estimate.getpseudo(collect0, n)
@everywhere pseudo1 = estimate.getpseudo(collect1, n)

pmap(1:2nρ) do i
	# idx1 = div(i - 1, nρ) + 1
	# idx2 = i - nρ * (idx1 - 1)
	# ρ0 = dρ * (idx1 - 1)
	# ρ1 = dρ * (idx2 - 1)
	ρ = (0:dρ:maxρ)[mod(i - 1, nρ) + 1]
	isρ0 = i <= nρ
	# θ = ρ + dρ
	filename = generate.getfilename("data", "5_16_26_f", isρ0, ρ)
	if !isfile(filename)
		println("ρ: $ρ, ρ0: $isρ0")
		ρhat = generate.repeated(collect0, collect1, pseudo0, pseudo1,
			n, l1, θ, J, dρ, maxρ,
			isρ0 * ρ, !isρ0 * ρ, covariate, pvec)
		save_object(filename, ρhat)
	end
end
