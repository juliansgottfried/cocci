module estimate

import StatsBase
import JLD2

initialize = function(n, l1)
    leavesn = falses(n, l1, 2)
    for i in 1:n leavesn[i, i, :] .= true end

    ancestor = falses(l1)
    ancestor[1:n] .= 1

    timen = zeros(Float64, l1)
    
    leavesb = falses(n, l1, 2)
    
    timeb = zeros(Float64, l1)

    (leavesn, ancestor, timen, leavesb, timeb)
end

expand = function(leavesn, ancestor, timen, leavesb, timeb, n, l1)
    newleavesn = falses(n, 2l1, 2)
    newleavesn[:, 1:l1, :] = leavesn

    newancestor = falses(2l1)
    newancestor[1:l1] = ancestor

    newtimen = zeros(Float64, 2l1)
    newtimen[1:l1] = timen

    newleavesb = falses(n, 2l1, 2)
    newleavesb[:, 1:l1, :] = leavesb

    newtimeb = zeros(Float64, 2l1)
    newtimeb[1:l1] = timeb

    (newleavesn, newancestor, newtimen, newleavesb, newtimeb, 2l1)
end

simulator = function(n, l1, ρ)
    leavesn, ancestor, timen, leavesb, timeb = initialize(n, l1)
    
    numbern = n
    numberb = 0
    time = 0.0
    K = n
    
    while K > 1
        # println(K)
        
        ρ_eff = 0.0
        hasleaves = (sum(leavesn, dims = 1) .> 0)[1, :, :]
        
        if sum(hasleaves[ancestor, 1]) > 1 && sum(hasleaves[ancestor, 2]) > 1
            ρ_eff = ρ * sum(hasleaves[ancestor, 1] .& hasleaves[ancestor, 2]) / K
        end
       
        λ = (K * (K - 1 + ρ_eff)) / 2
        time -= log(1 - rand()) / λ
        
        if rand() > ρ_eff / (K - 1 + ρ_eff)
            chosen = StatsBase.sample((1:l1)[ancestor], 2, replace = false)
            ancestor[chosen] .= 0
            numbern += 1
            numberb += 2
            
            if numbern > l1 || numberb > l1
                leavesn, ancestor, timen, leavesb, timeb, l1 = expand(leavesn, ancestor, timen, leavesb, timeb, n, l1)
            end
            
            ancestor[numbern] = 1
            timen[numbern] = time
            leavesn[:, numbern, :] = leavesn[:, chosen[1], :] .| leavesn[:, chosen[2], :]
            
            timeb[[numberb - 1; numberb]] = timen[numbern] .- timen[chosen]
            leavesb[:, [numberb - 1; numberb], :] = leavesn[:, chosen, :]
        else
            chosen = StatsBase.sample((1:l1)[ancestor .& hasleaves[:, 1] .& hasleaves[:, 2]])
            ancestor[chosen] = 0
            numbern += 2
            numberb += 1
            
            if numbern > l1 || numberb > l1
                leavesn, ancestor, timen, leavesb, timeb, l1 = expand(leavesn, ancestor, timen, leavesb, timeb, n, l1)
            end
            
            ancestor[[numbern - 1; numbern]] .= 1
            timen[[numbern - 1; numbern]] .= time
            leavesn[:, numbern - 1, 1] = leavesn[:, chosen, 1]
            leavesn[:, numbern, 2] = leavesn[:, chosen, 2]

            timeb[numberb] = time - timen[chosen]
            leavesb[:, numberb, :] = leavesn[:, chosen, :]
        end

        K = sum(ancestor)
    end

    # println(K)

    (leavesb, timeb, numberb)
end

geth! = function(blackidx, whiteidx, config3, leavesb, timeb)
    for i in blackidx, j in whiteidx
        shared = sum(leavesb[:, i, 1] .& leavesb[:, j, 2])
        config3[2][shared .== config3[1]] .+= timeb[i] * timeb[j]
    end
end

configloop = function(allconfigs, n, leavesb, timeb, numberb)
    checker = falses(numberb, n - 1, 2)
    numberleaves = sum(leavesb, dims = 1)[1, :, :]
        
    for i in 1:(n - 1)
        checker[:, i, :] = numberleaves[1:numberb, :] .== i
    end
    
    for config1 in 1:(n - 1)
        blackidx = (1:numberb)[checker[:, config1, 1]]
        for config2 in 1:(n - 1)
            whiteidx = (1:numberb)[checker[:, config2, 2]]
            geth!(blackidx, whiteidx, allconfigs[config1, config2], leavesb, timeb)
        end
    end

    allconfigs
end

buildallconfigs = function(i, j, n)
    k = max(0, i + j - n):min(i, j)
    [k, zeros(Float64, size(k)[1])]
end

sumallconfigs! = function(allconfigs, add, n, m)
    for i in 1:(n - 1), j in 1:(n - 1)
        allconfigs[i, j][2] += add[i, j][2] ./ m
    end
end

montecarlo = function(n, l1, ρ, m)
    allconfigs = [buildallconfigs(i, j, n) for i = 1:(n - 1), j = 1:(n - 1)]

    for i in 1:m
        mod(i,1000) == 0 && println("Iteration ", i)
        leavesb, timeb, numberb = simulator(n, l1, ρ)
        add = configloop(allconfigs, n, leavesb, timeb, numberb)
        sumallconfigs!(allconfigs, add, n, m)
    end

    allconfigs
end

rhogrid = function(n, l1, ρs, m)
    for i in eachindex(ρs)
        println("ρ: ", ρs[i])
        results = montecarlo(n, l1, ρs[i], m)
        JLD2.save_object(string("/scratch/users/jgottf/cocci/results/results_",replace(string(ρs[i]), "." => "_"), ".jld2"), [results, (n, l1, m, ρs[i])])
    end

    results
end

getsum = function(allconfigs, n, pseudo)
    summer = 0
    for i in 1:(n - 1), j in 1:(n - 1)
        summer += sum(allconfigs[i, j][2] .+ pseudo)
    end

    summer
end

getl = function(ρs, n, results, configs, dists, pseudo)
    l = size(ρs)[1]
    loglik = zeros(Float64, l)
    for i in 1:l
        # println("i: ",i)
        denom = getsum(results[i], n, pseudo)
        for j in 1:size(configs)[1]
            # println("j: ",j)
            ρidx = argmin(((dists[j] * ρs[i]) .- ρs) .^ 2)
            extract = results[ρidx][configs[j, 1], configs[j, 2]]
            numer = extract[2][extract[1] .== configs[j, 3]][1] + pseudo
            qc = log(numer / denom)
            if !isinf(qc) loglik[i] += qc 
            else 0 end
        end
    end

    loglik
end

end
