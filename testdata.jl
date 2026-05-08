include("estimate.jl")
include("generate.jl")
using JLD2

n = 25
l1 = 50
maxρ = 1
ρs = generate.makegrid(maxρ)
filenames = generate.getfilename.("prob", "5_8_26", ρs)

collect = []
for i in filenames
	push!(collect, load_object(i)[1])
end

pseudo = estimate.getpseudo(collect, n)

dt = 0.01
maxtime = 4
times = 0:dt:maxtime
growth = -2
covariate = generate.buildcov(dt, maxtime, growth)

J = 10
θ = 10

pmap(ρs) do ρ
	# θ = iszero(ρ) ? 0.1 : ρ
	ρhat = generate.repeated(ρs, collect, pseudo, n, l1, ρ[1], ρ[2], covariate, dt, θ, J)
	save_object(generate.getfilename("data", "5_8_26", ρ), 
	[ρhat, (n, ρ[1], ρ[2], θ, J, dt, maxtime, growth)])
end
