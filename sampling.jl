using Plots
using Distributions
using LinearAlgebra
using Colors

function prob(m1, m2, K)
    i = max(0, K - m2):min(K, m1)
    term1 = log.(binomial.(big(m1), big.(i))) .+ log.(binomial.(big(m2), big.(K .- i)))
    term2 = log.((factorial.(big.(m1 .+ K .- 2i)))) .+ log.(factorial.(big.(m2 .- K .+ 2i))) .- log.(factorial.(big.(m1 + m2 + 1)))
    sum(exp.(term1 .+ term2))
end

size = 50
probs = zeros(size, size, 2size + 1)
for i in 1:size
    for j in 1:i
        ks = 1:(ceil(Int, (i + j - 2) / 2) + 1)
        for k in ks
            probs[i, j, k] = prob(i - 1, j - 1, k - 1)
            probs[i, j, i + j - k] = probs[i, j, k]
        end
    end
end
for k in 1:(2size + 1)
    probs[:, :, k] = probs[:, :, k] .+ probs[:, :, k]' .- diagm(diag(probs[:, :, k]))
end

# m1 between 0 and 50
# m2 between 0 and 50
# K betweem 0 and m1 + m2
# retrieve(m1, m2, K) = probs[max(m1, m2) + 1, min(m1, m2) + 1, K + 1]

# plot(0:10, freq[1:11],
#         ylim = [0, 1],
#         label = false, color = :black)

# strains = 2
# segs = 5
# distr = Multinomial(segs, (1/strains) * ones(strains))
# rand(distr)

function segfreq(G, η)
    freq = [1; zeros(size - 1)]
    # sums = zeros(G)
    for i in 1:G
        pairs = freq * freq'
        freq *= η
        for j in 1:(size - 1)
            freq[j + 1] += (1 - η) * sum(probs[:, :, j] .* pairs)
        end
        # sums[i] = sum(freq)
        freq ./= sum(freq)
    end
    # freq[2:size] ./ sum(freq[2:size])
    freq ./ sum(freq)
    # freq[2:size]
    # sums
end


G = 365
ηs = 0.95:0.0001:0.9999
n_η = length(ηs)

out = segfreq.(G, ηs)
seg_mat = reduce(hcat, out)

seg_mat_rev = seg_mat[:, n_η:-1:1]

kernels = seg_mat_rev[2:50, :]
for i in 1:500
    kernels[:, i] ./= sum(kernels[:, i])
end
kernels
plot(kernels)
plot(kernels[:, 400])

L = 0:0.01:1
accum = ones(Float64, 500, 101)
for i in 1:500
    vec = seg_mat_rev[:, i]
    for j in 0:49
        accum[i, :] .-= vec[j + 1] .* (1 .- L) .^ j
    end
end
for i in 1:500
    accum[i, :] ./= maximum(accum[i, :])
end

colors = ["#EB2F89", "#EACD55"]
l = @layout [a{0.95w} b]
p1 = plot(L, accum',
            xlim = [0, 1], ylim = [0, 1];
            palette = palette(colors, 500),
            alpha = 0.25, label = false,
            linewidth = 0.8,
            xlabel = "distance between loci", ylabel= "two-locus ρ",
            grid = false)
plot!(L, accum'[:, [1; 100; 200; 300; 400; 500]],
            xlim = [0, 1], ylim = [0, 1];
            color = :black, linewidth = 0.6,
            alpha = 0.8,
            label = false,
            grid = false)
p2 = heatmap(rand(2, 2), clims = (0.1, 50),
            framestyle = :none,
            c = palette(colors, ηs),
            cbar = true,
            lims = (-1, 0),
            cbartitle = "1000 * r")
plot(p1, p2, layout=l)




# heatmap(1:cutoff, ηs, transpose(seg_mat),
#         color = cgrad(:devon, rev = false),
#         xlabel= " # junctions", ylabel = "selfing rate η")

function plot_curve(i)
    plot(1:(size - 1), seg_mat_rev[:, 1:(i - 1)],
            xlim = [1, 8], ylim = [0, 1];
            color = :black,
            linewidth = 0.2,
            linealpha = 0.5,
            label = false,
            grid = false,
            xlabel = "# junctions", ylabel=  " density",
            guidefont = font(15, "Mono"),
            tickfont = font(10, "Mono"))
    plot!(1:(size - 1), seg_mat_rev[:, i],  
                xlim = [1, 8], ylim = [0, 1];
                color = palette(colors, 100)[i],
                alpha = 1, label = false,
                title = string("Selfing rate: ", ηs[end - i + 1]),
                titlefontcolor = string("#",hex(RGB(palette(colors, 100)[i]))),
                linewidth = 4)
end


colors = ["#EB2F89", "#EACD55"]
l = @layout [a{0.95w} b]
p1 = plot(1:(size - 1), seg_mat_rev[:, 101:end],
            xlim = [1, 8], ylim = [0, 1];
            color = :black,
            linewidth = 0.2,
            linealpha = 0.1,
            label = false)
plot!(1:(size - 1), seg_mat_rev[:, 1:100],
            xlim = [1, 8], ylim = [0, 1];
            palette = palette(colors, 100),
            alpha = 0.8, label = false,
            linewidth = 0.8,
            xlabel = "# junctions",ylabel= " density")
plot!(1:(size - 1), seg_mat_rev[:, [1; 25; 50; 75; 100]],
            xlim = [1, 8], ylim = [0, 1];
            color = :black, linewidth = 0.7,
            alpha = 0.7,
            label = false)
p2 = heatmap(rand(2, 2), clims = (0.1, 10),
            framestyle = :none,
            c = palette(colors, ηs[401:end]),
            cbar = true,
            lims = (-1, 0),
            cbartitle = "(1-η)e3")
plot(p1, p2, layout=l)

Gs = 10:1:10100
N = 100
η = Gs ./ N

out = segfreq_2.(Gs, N)
seg_mat = reduce(hcat, out)

# heatmap(1:(size - 1), Gs, transpose(seg_mat),
#         color = cgrad(:devon, rev = false),
#         xlabel= " # junctions", ylabel = "G")

colors = ["#EB2F89", "#EACD55"]
l = @layout [a{0.95w} b]
p1 = plot(1:(size - 1), seg_mat,
            xlim = [1, 6],
            palette = palette(colors, length(Gs)),
            alpha=0.5, label = false,
            linewidth=0.2,
            xlabel="# junctions",ylabel="density")
plot!(1:(size - 1), seg_mat[:, 650:50:850],
            xlim = [1, 6],
            color = :black, alpha = 0.5,
            label = false)
p2 = heatmap(rand(2, 2), clims = (0.1, 10.1),
            framestyle = :none,
            c = palette(colors, η),
            cbar = true,
            lims = (-1, 0),
            cbartitle = "η")
plot(p1, p2, layout=l)

strains = 3
junctions = 20
rand(Multinomial(junctions, ones(strains) ./ strains))
