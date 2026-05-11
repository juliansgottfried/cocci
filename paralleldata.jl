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

getobj = function(dρ, idx1, idx2)
	ρ0 = dρ * (idx1 - 1)
	ρ1 = dρ * (idx2 - 1)
	filename = generate.getfilenamelocal("prob", "5_8_26", ρ0, ρ1)
	load_object(filename)[1]
end

collect = [getobj(dρ, idx1, idx2) for idx1 in 1:nρ, idx2 in 1:nρ]

@everywhere pseudo = estimate.getpseudo(collect, n)

@everywhere dt = 0.01
@everywhere maxtime = 3
@everywhere change = 20
@everywhere covariate = generate.buildcov(dt, maxtime, change)

@everywhere J = 100
@everywhere θ = 10

pmap(ρs) do ρ
	# θ = iszero(ρ[1] + ρ[2]) ? 0.1 : ρ[1] + ρ[2]
	ρhat = generate.repeated(ρs, collect, pseudo, n, l1, ρ[1], ρ[2], covariate, dt, θ, J)
	save_object(generate.getfilename("data", "5_8_26", ρ0, ρ1), [ρhat])
end