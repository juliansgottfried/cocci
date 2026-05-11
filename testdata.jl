include("estimate.jl")
include("generate.jl")
using JLD2

n = 25
l1 = 50
maxρ = 30
ρs = generate.makegrid(maxρ)
filenames = generate.getfilenamelocal.("prob", "5_8_26", ρs)

collect = []
for i in filenames
	push!(collect, load_object(i)[1])
end

pseudo = estimate.getpseudo(collect, n)

dt = 0.01
maxtime = 4
change = -2980.9579870417
covariate = generate.buildcov(dt, maxtime, change)

J = 100
θ = 10

ρ = ρs[100]
configs, dists = generate.getconfigs(n, l1, ρ[1], ρ[2], covariate, dt, θ)




getstore = function(collect, l, config, pseudo)
	store = zeros(Float64, l)
	for i in 1:l
		extract = collect[i][config[1], config[2]]
		store[i] = extract[2][extract[1] .== config[3]][1]
	end
	store .+= pseudo
	if config[1] != config[2] store .*= 2 end
	store
end

l = length(ρs)
loglik = zeros(Float64, l)
denom = zeros(Float64, l)
for i in 1:l denom[i] = estimate.getsum(collect[i], n, pseudo) end

pairedconfigs = [configs dists]
sortedpairedconfigs = sortslices(pairedconfigs, dims = 1, by = x -> (x[1], x[2], x[3]), rev = false)
sortedconfigs = Int.(sortedpairedconfigs[:, 1:3])
sorteddists = sortedpairedconfigs[:, 4]

ρmat = reduce(hcat, ρs)'

currentconfig = [0; 0; 0]
for i in 1:size(sortedconfig)[1]
	println(i)
	tmpconfig = sortedconfigs[i, :]
	if tmpconfig != currentconfig
		currentconfig = tmpconfig
		store = getstore(collect, l, currentconfig, pseudo)
	end
	
	scaled = sorteddists[i] * ρmat
	relocate = zeros(Int, l)
	for j in 1:l
		a = (scaled[j, 1] .- ρmat[:, 1]) .^ 2
		b = (scaled[j, 2] .- ρmat[:, 2]) .^ 2
		relocate[j] = (1:length(ρs))[(a .== minimum(a)) .& (b .== minimum(b))][1]
	end

	loglik .+= log.(store[relocate] ./ denom[relocate])
end


spacing = 0.1
fakeρ = 0:spacing:10
dist = 0.3
Int.(div.(fakeρ .* dist, spacing)) .+ 1
# and could do that for both indices


dρ = 0.5
maxρ = 10
0:dρ:maxρ

i = 1
ρi = dρ*(i-1)
nρ = Int(div(maxρ, dρ)) + 1


# but we need a way to store results
# well, it should in fact be in a matrix
# yes
# rho0 indices along the rows, rho1 along the columns
# so we can expect it to be a matrix


# getl is the only function to change in estimate
# repeated is the only function to change in generate

# and the parallel files



for i in 1:nρ^2
	idx1 = div(i, nρ + 1) + 1
	idx2 = i - nρ * (idx1 - 1)
	
end

# figure out how many times nρ fits into it. that's i
# then get modulo; that's j