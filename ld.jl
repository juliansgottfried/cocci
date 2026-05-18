using Plots

include("kernels.jl")

getrsq = function(ρidx, prob0, n)
    configs = []
    probs = []
    for i in 1:(n - 1), j in 1:i 
        for k in eachindex(prob0[ρidx][i, j][1])
            push!(configs, [i; j; prob0[ρidx][i, j][1][k]])
            push!(probs, [prob0[ρidx][i, j][2][k]])
        end
    end
    configs = reduce(hcat, configs)'
    probs = reduce(vcat, probs)
    probs ./= sum(probs)
    freqs = configs ./ n
    D = freqs[:, 3] .- freqs[:, 1] .* freqs[:, 2]
    rsq = D .^ 2 ./ (freqs[:, 1] .* (1 .- freqs[:, 1]) .* freqs[:, 2] .* (1 .- freqs[:, 2]))
    sum(probs .* rsq)
end

fetchidx = function(whichkernel, accum, ρ, dρ)
    Int.(div.(ρ .* accum[whichkernel, :], dρ)) .+ 1
end

prob0 = [load_object(generate.getfilenamelocal("prob", "5_15_26_d", true, ρ)) for ρ in 0:dρ:maxρ]
n = 17
dists = 0:0.01:1
rsq = [getrsq(i, prob0, n) for i in 1:nρ]

G = 365
ηs = 0.95:0.0001:0.9999
maxbrk = 50

kernelcollect = kernels.getkernels(G, ηs, maxbrk)

nη = length(ηs)

accum = ones(Float64, nη, length(dists))
for i in 1:nη
    for j in 1:(maxbrk - 1)
        accum[i, :] .-= kernelcollect[:, i][j] .* (1 .- dists) .^ j
    end
end

whichkernel = 1:20:nη
ρidx = [fetchidx(kern, accum, maxρ, dρ) for kern in whichkernel]

plotrsq = zeros(Float64, length(dists), length(whichkernel))
for i in eachindex(whichkernel)
    plotrsq[:, i] = [rsq[ρidx[i][j]] for j in eachindex(dists)]
end

plot(dists, plotrsq,
    palette = palette([:red, :black], length(whichkernel)),
    linewidth = 0.5, alpha = 1,
    xlabel = "distance",
    ylabel = "r²",
    label = false, grid = false)
