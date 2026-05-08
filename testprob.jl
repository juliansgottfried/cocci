include("estimate.jl")
include("generate.jl")
using JLD2

n = 25
l1 = 50
m = 5
maxρ = 1
ρs = generate.makegrid(maxρ)

dt = 0.01
maxtime = 4
times = 0:dt:maxtime
growth = -2
covariate = generate.buildcov(dt, maxtime, growth)

map(ρs) do ρ
	results = estimate.montecarlo(n, l1, ρ[1], ρ[2], covariate, dt, m)
	save_object(generate.getfilename("prob", "5_8_26", ρ), 
		[results, (n, m, ρ[1], ρ[2], dt, maxtime, growth)])
end
