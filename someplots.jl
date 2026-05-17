using Plots, JLD2
import Statistics, Distributions

include("estimate.jl")
include("generate.jl")
include("kernels.jl")

dρ = 0.1
maxρ = 10 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ

# Retrodiction plots

loadit(filename) = if isfile(filename) return load_object(filename) end

data0 = [loadit(generate.getfilenamelocal("data", "5_16_26_g", true, ρ)) for ρ in 0:dρ:maxρ]
data1 = [loadit(generate.getfilenamelocal("data", "5_16_26_g", false, ρ)) for ρ in 0:dρ:maxρ]
skip = isnothing.(data0) .| isnothing.(data1)

aic = zeros(Float64, nρ, 2, 3)
raw = zeros(Float64, nρ, 4, 3)
qs = [0.055, 0.5, 0.945]
for i in 1:nρ
    if skip[i] 
        aic[i, :, :] .= -Inf
        raw[i, :, :] .= -Inf
    else
        aic[i, 1, :] = Statistics.quantile(2(data0[i][:, 2] - data0[i][:, 4]), qs)
        aic[i, 2, :] = Statistics.quantile(2(data1[i][:, 4] - data1[i][:, 2]), qs)
       
        raw[i, 1, :] = Statistics.quantile(data0[i][:, 1], qs)
        raw[i, 2, :] = Statistics.quantile(data0[i][:, 3], qs)
        raw[i, 3, :] = Statistics.quantile(data1[i][:, 1], qs)
        raw[i, 4, :] = Statistics.quantile(data1[i][:, 3], qs)
    end
end

aicplot = function(i)
    plot(ρs, aic[:, i, 1], fillrange = aic[:, i, 3], 
        xlabel = "recombination rate", ylabel = "Δ AIC",
        title = i == 1 ? "constant data" : "varying data",
        xlim = [0, maxρ], fillalpha = 0.15, c = "#29a0c8",
        linecolor = false, label = false, grid = false)
    hline!([2], c = :red, alpha = 0.7, label = false)
    plot!(ρs, aic[:, i, 2], c = :black, linewidth = 1.2, label = false)
end

rawplot = function(i, j)
    k = 2(i - 1) + j
    part1 = i == 1 ? "constant data" : "varying data"
    part2 = j == 1 ? "constant model" : "varying model"
    plot(ρs, raw[:, k, 1], fillrange = raw[:, k, 3], 
        xlabel = "recombination rate", ylabel = "estimate",
        title = "$part1, $part2",
        xlim = [0, maxρ], fillalpha = 0.15, c = "#29a0c8",
        linecolor = false,
        label = false, grid = false)
    plot!(ρs, maxρ * ones(nρ),
        linestyle = :dash, color = "#29a0c8", label = false)
    plot!(ρs, zeros(nρ), 
        linestyle = :dash, color = "#29a0c8", label = false)
    plot!(ρs, raw[:, k, 2], c = :black, 
        linewidth = 1.2, label = false)
    plot!([0, maxρ], [0, maxρ], 
        color = :red, label = false)
end

aicplot(1)
aicplot(2)

rawplot(1, 1)
rawplot(1, 2)
rawplot(2, 1)
rawplot(2, 2)


# LD plots

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
