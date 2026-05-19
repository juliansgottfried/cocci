#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2, DelimitedFiles

@everywhere n = 17
@everywhere l1 = 25
@everywhere θ = 20
@everywhere nsample = 50
@everywhere J = 500

@everywhere dρ = 0.5
@everywhere maxρ = 20 - dρ
@everywhere ρs = 0:dρ:maxρ
@everywhere nρ = length(ρs)

@everywhere covariate = readdlm("/scratch/users/jgottf/cocci/rodent_data/covariate.csv", ',', Any, '\n')

#= @everywhere dt = 0.01
@everywhere maxtime = 1
@everywhere change = 30
@everywhere covariate = generate.buildcov(dt, maxtime, change) =#

@everywhere pvec = 1

# @everywhere collect0 = [load_object(generate.getfilename("prob", "5_18_26_e", true, ρ)) for ρ in ρs]
# @everywhere collect1 = [load_object(generate.getfilename("prob", "5_18_26_e", false, ρ)) for ρ in ρs]
@everywhere collect = [load_object(generate.getfilenamegrid("prob", "5_19_26_a", ρ0, ρ1)) for ρ0 in ρs, ρ1 in ρs]

# @everywhere pseudo0 = estimate.getpseudo(collect0, n)
# @everywhere pseudo1 = estimate.getpseudo(collect1, n)
# @everywhere pseudo = min(pseudo0, pseudo1)
@everywhere pseudo = estimate.getpseudo(collect, n)

#= pmap(1:2nρ) do i
	ρ = (0:dρ:maxρ)[mod(i - 1, nρ) + 1]
	isρ0 = i <= nρ
	filename = generate.getfilename("data", "5_18_26_d", isρ0, ρ)
	if !isfile(filename)
		ρhat = generate.repeated(collect0, collect1, pseudo, pseudo,
			n, l1, θ, nsample, J, dρ, maxρ,
			isρ0 * ρ, !isρ0 * ρ, covariate, pvec)
		save_object(filename, ρhat)
	end
end =#

pmap(1:nρ^2) do i
	ρ0 = ρs[Int(div(i - 1, nρ)) + 1]
	ρ1 = ρs[Int(mod(i - 1, nρ)) + 1]
	filename = generate.getfilenamegrid("data", "5_19_26_a", ρ0, ρ1)
	if !isfile(filename)
		ρhat = generate.repeatedgrid(collect, pseudo,
			n, l1, θ, nsample, J, dρ, maxρ,
			ρ0, ρ1, covariate, pvec)
		save_object(filename, ρhat)
	end
end
