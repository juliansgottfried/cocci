using Plots, JLD2
import Statistics, Distributions

include("estimate.jl")
include("generate.jl")

dρ = 0.1
maxρ = 20 - dρ
nρ = length(0:dρ:maxρ)
ρs = 0:dρ:maxρ

loadit(filename) = if isfile(filename) return load_object(filename) end

data0 = [loadit(generate.getfilenamelocal("data", "5_15_26_d", true, ρ)) for ρ in 0:dρ:maxρ]
data1 = [loadit(generate.getfilenamelocal("data", "5_15_26_d", false, ρ)) for ρ in 0:dρ:maxρ]
skip = isnothing.(data0) .| isnothing.(data1)

aic = zeros(Float64, nρ, 2, 3)
ratio = zeros(Float64, nρ, 4, 3)
qs = [0.055, 0.5, 0.945]
for i in 1:nρ
    if skip[i] 
        aic[i, :, :] .= -Inf
        ratio[i, :, :] .= -Inf
    else
        aic[i, 1, :] = Statistics.quantile(2(data0[i][:, 2] - data0[i][:, 4]), qs)
        aic[i, 2, :] = Statistics.quantile(2(data1[i][:, 4] - data1[i][:, 2]), qs)
        if ρs[i] == 0 
            ratio[i, :, :] .= Inf
        else
            tmp0 = data0[i][:, [1; 3]]
            tmp0[tmp0 .== 0] .= 0.1
            tmp1 = data1[i][:, [1; 3]]
            tmp1[tmp1 .== 0] .= 0.1
            ratio[i, 1, :] = Statistics.quantile(log2.(tmp0[:, 1] ./ ρs[i]), qs)
            ratio[i, 2, :] = Statistics.quantile(log2.(tmp0[:, 2] ./ ρs[i]), qs)
            ratio[i, 3, :] = Statistics.quantile(log2.(tmp1[:, 1] ./ ρs[i]), qs)
            ratio[i, 4, :] = Statistics.quantile(log2.(tmp1[:, 2] ./ ρs[i]), qs)
        end
    end
end
# ratio[ratio .== -Inf] .= minimum(ratio[ratio .!= -Inf])
# ratio[ratio .== Inf] .= maximum(ratio[ratio .!= Inf])

aicplot = function(i)
    plot(ρs, aic[:, i, 1], fillrange = aic[:, i, 3], 
        xlabel = "recombination rate", ylabel = "Δ AIC",
        title = i == 1 ? "constant data" : "varying data",
        xlim = [0, maxρ], fillalpha = 0.15, c = "#29a0c8",
        linecolor = false, label = false, grid = false)
    hline!([2], c = :red, alpha = 0.7, label = false)
    plot!(ρs, aic[:, i, 2], c = :black, linewidth = 1.2, label = false)
end

ratioplot = function(i, j)
    k = 2(i - 1) + j
    part1 = i == 1 ? "constant data" : "varying data"
    part2 = j == 1 ? "constant model" : "varying model"
    plot(ρs, ratio[:, k, 1], fillrange = ratio[:, k, 3], 
        xlabel = "recombination rate", ylabel = "log2 relative error",
        title = "$part1, $part2",
        xlim = [0, maxρ], fillalpha = 0.15, c = "#29a0c8",
        linecolor = false,
        label = false, grid = false)
    hline!([0], c = :red, alpha = 0.7, label = false)
    plot!(ρs, log2.(maxρ ./ ρs), 
        linestyle = :dash, color = "#29a0c8", label = false)
    plot!(ρs, log2.(dρ ./ ρs), 
        linestyle = :dash, color = "#29a0c8", label = false)
    plot!(ρs, ratio[:, k, 2], c = :black, 
        linewidth = 1.2, label = false)
end

aicplot(1)
aicplot(2)

ratioplot(1, 1)
ratioplot(1, 2)
ratioplot(2, 1)
ratioplot(2, 2)

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

rsq = [getrsq(i, prob0, n) for i in 1:nρ]

getidx(whichkernel) = Int.(div.(ρ .* accum[whichkernel, :], dρ)) .+ 1

prob0 = [load_object(generate.getfilenamelocal("prob", "5_15_26_d", true, ρ)) for ρ in 0:dρ:maxρ]
n = 17
ρ = 19.9
dists = 0:0.01:1

nkern = 500
nbrks = 49
accum = ones(Float64, nkern, length(dists))
for i in 1:nkern
    for j in 1:nbrks
        accum[i, :] .-= kernels[:, i][j] .* (1 .- dists) .^ j
    end
end

whichkernel = 1:20:nkern
ρidx = getidx.(whichkernel)

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
