#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@everywhere include("/scratch/users/jgottf/cocci/estimate.jl")
@everywhere include("/scratch/users/jgottf/cocci/generate.jl")
@everywhere using JLD2

@everywhere n = 50
@everywhere l1 = 100
@everywhere maxρ = 30
@everywhere ρs = generate.makegrid(maxρ)
@everywhere filenames = getfilename.("prob", "5_6_26", ρs)

@everywhere collect = []
@everywhere for i in filenames
	push!(collect, load_object(i)[1])
end

@everywhere pseudo = estimate.getpseudo(collect, n)

@everywhere dt = 0.01
@everywhere maxtime = 4
@everywhere times = 0:dt:maxtime
@everywhere growth = -2
@everywhere covariate = generate.buildcov(dt, maxtime, growth)

@everywhere J = 1000
@everywhere θ = 10

pmap(ρs) do ρ
	# θ = iszero(ρ) ? 0.1 : ρ
	ρhat = generate.repeated(ρs, collect, pseudo, n, l1, ρ[1], ρ[2], covariate, dt, θ, J)
	save_object(generate.getfilename("data", "5_6_26", ρ), 
	[ρhat, (n, ρ[1], ρ[2], θ, J, dt, maxtime, growth)])
end
