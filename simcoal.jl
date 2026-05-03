using Plots
using DelimitedFiles
using Random
using Distributions
Random.seed!(999)

function simulate(L, K, ρ, χ)
    T = 0
    RA = 0
    seqs = Dict(i => ones(Bool, L) for i = 1:K)
    counts = K .* ones(Int, L)
    last_idx = K + 1
    
    while K > 1
        rate = (K * (K - 1 + ρ)) / 2
        T += rand(Exponential(rate))
        u = rand()
        p_CA = (K - 1) / (K - 1 + ρ)
        
        if (u < p_CA)
            idx = sample(collect(keys(seqs)), 2, replace = false)
            new_seq = seqs[idx[1]] .|| seqs[idx[2]]
            delete!(seqs,idx[1])
            delete!(seqs,idx[2])
            seqs[last_idx] = new_seq
            K -= 1
            last_idx += 1
            counts -= new_seq
            remove = counts .!= 1
            
            for key in keys(seqs) 
                seqs[key] = seqs[key] .&& remove
                
                if (sum(seqs[key]) == 0)
                    delete!(seqs,key)
                    K -= 1
                end
            end
        else
            idx = sample(collect(keys(seqs)), 1)
            new_seq = seqs[idx[1]]
            delete!(seqs,idx[1])
            breaks = [1; sort(sample(1:L, χ - 1, replace = false)); L + 1]
            recomb = zeros(χ, L)
            keep = zeros(Bool, χ)
            
            for i in 2:(χ + 1)
                start = breaks[i - 1]
                finish = breaks[i] - 1
                segment = new_seq[start:finish]
                
                if (sum(segment) > 0)
                    recomb[i - 1, start:finish] = segment
                    keep[i - 1] = 1
                end
            end
            recomb_cut = recomb[keep, :]
            n_parents = size(recomb_cut)[1]
            
            for i in 1:n_parents
                seqs[last_idx] = recomb_cut[i, :]
                last_idx += 1
            end
            K += n_parents - 1
            RA += n_parents - 1
        end
        # println(collect(keys(seqs)))
        # println(K)
    end
    [T, RA]
end


L = 100 # locus length ≥ 1
K = 50 # sample size ≥ 1
# ρ = 1 # population recombination rate ≥ 0
# χ = 2 # multifurcation number ≥ 2

# out = simulate(L, K, ρ, χ)[1]

scale(vec) = (vec .- minimum(vec)) ./ (maximum(vec) .- minimum(vec))

ρs = 1:0.1:5
χs = 2:10
len_ρ = length(ρs)
len_χ = length(χs)
Ts_mean = zeros(len_ρ, len_χ)
Ts_var = zeros(len_ρ, len_χ)
Ts_mean_scale = zeros(len_ρ, len_χ)
Ts_var_scale = zeros(len_ρ, len_χ)
N = 10000

for j in 1:len_χ
    for i in 1:len_ρ
        sims = zeros(N)
        for k in 1:N
            sims[k] = simulate(L, K, ρs[i], χs[j])[1]
        end
        Ts_mean[i, j] = mean(sims)
        Ts_var[i, j] = var(sims)
    end
    Ts_mean_scale[:, j] = scale(Ts_mean[:, j])
    Ts_var_scale[:, j] = scale(Ts_var[:, j])
end

heatmap(ρs,χs,transpose(Ts_mean),
        color=cgrad(:rainbow,rev=false),
        xlabel="ρ",ylabel="χ",
        plot_title="E[T₅₀]")
heatmap(ρs,χs,transpose(Ts_var_scale),
        color=cgrad(:rainbow,rev=false),
        xlabel="ρ",ylabel="χ",
        plot_title="Var[T₅₀]")

# writedlm("T_mean.csv", Ts_mean, ",")
