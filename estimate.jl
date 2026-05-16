module estimate

import StatsBase, JLD2

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

getintegral = function(covariate, dt, time)
    if time >= covariate[end, 1] 
        return sum(covariate[:, 2]) * dt + covariate[end, 2] * (time - covariate[end, 1])
    end
    timecut = Int(time ÷ dt) + 1
    accum = 0
    if timecut > 1 accum = sum(covariate[1:(timecut - 1), 2]) * dt end
    accum + covariate[timecut, 2] * mod(time, dt)
end

getvalue = function(covariate, dt, time)
    timecut = Int(time ÷ dt) + 1
    if timecut >= size(covariate)[1] return covariate[end, 2] end
    covariate[timecut, 2]
end

Λ = function(K, ρ0, ρ1, time, covariate, dt)
    0.5K * ((K - 1 + ρ0) * time + ρ1 * getintegral(covariate, dt, time))
end

targetf = function(K, ρ0, ρ1, next, covariate, dt, time)
    0.5K * (K - 1 + ρ0 + ρ1 * getvalue(covariate, dt, next)) * exp(Λ(K, ρ0, ρ1, time, covariate, dt) - Λ(K, ρ0, ρ1, next, covariate, dt))
end

M = function(λ, K, ρ1, covariate, dt, time)
    (λ + 0.5K * ρ1 * getvalue(covariate, dt, time)) / (λ * exp(-λ * time))
end

simulator = function(n, l1, ρ0, ρ1, covariate)
    leavesn, ancestor, timen, leavesb, timeb = initialize(n, l1)
    
    dt = covariate[2, 1] - covariate[1, 1]

    numbern = n
    numberb = 0
    time = 0.0
    next = 0.0
    K = n
    
    while K > 1
        # println(K)
        
        ϕ = 0.0
        hasleaves = (sum(leavesn, dims = 1) .> 0)[1, :, :]
        if sum(hasleaves[ancestor, 1]) > 1 && sum(hasleaves[ancestor, 2]) > 1
            ϕ = sum(hasleaves[ancestor, 1] .& hasleaves[ancestor, 2]) / K
        end
       
        λ = 0.5K * (K - 1 + ϕ * ρ0)
        comparison = false
        while !comparison
            next = time - log(1 - rand()) / λ
            comparison = rand() <= targetf(K, ϕ * ρ0, ϕ * ρ1, next, covariate, dt, time) / (M(λ, K, ϕ * ρ1, covariate, dt, time) * λ * exp(-λ * next))
        end
        time = next
        
        if rand() < (K - 1) / (K - 1 + ϕ * ρ0 + ϕ * ρ1 * getvalue(covariate, dt, time))
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

geth! = function(blackidx, whiteidx, thisconfig, leavesb, timeb)
    for i in blackidx, j in whiteidx
        shared = sum(leavesb[:, i, 1] .& leavesb[:, j, 2])
        thisconfig[2][shared .== thisconfig[1]] .+= timeb[i] * timeb[j]
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
        for config2 in 1:config1
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
    for i in 1:(n - 1), j in 1:i
        allconfigs[i, j][2] += add[i, j][2] ./ m
    end
end

montecarlo = function(n, l1, m, ρ0, ρ1, covariate)
    allconfigs = [buildallconfigs(i, j, n) for i = 1:(n - 1), j = 1:(n - 1)]

    for i in 1:m
        mod(i, 1000) == 0 && println("Iteration ", i)
        leavesb, timeb, numberb = simulator(n, l1, ρ0, ρ1, covariate)
        add = configloop(allconfigs, n, leavesb, timeb, numberb)
        sumallconfigs!(allconfigs, add, n, m)
    end

    allconfigs
end

getsum = function(allconfigs, n, pseudo)
    summer = 0
    for i in 1:(n - 1), j in 1:i
        tmpsum = sum(allconfigs[i, j][2] .+ pseudo)
        if i != j tmpsum *= 2 end
        summer += tmpsum
    end

    summer
end

orderconfigs = function(configs, dists)
	configs = [sort(configs[:, 1:2], dims = 2, rev = true) configs[:, 3]]
	pairedconfigs = [configs dists]
	sortedpairedconfigs = sortslices(pairedconfigs, dims = 1, by = x -> (x[1], x[2], x[3]), rev = false)
	sortedconfigs = Int.(sortedpairedconfigs[:, 1:3])
	sorteddists = sortedpairedconfigs[:, 4]
	(sortedconfigs, sorteddists)
end

getstore = function(collect, nρ, config, pseudo)
	store = zeros(Float64, nρ)
	for i in 1:nρ
		extract = collect[i][config[1], config[2]]
		store[i] = extract[2][extract[1] .== config[3]][1]
	end
	store .+= pseudo
	if config[1] != config[2] store .*= 2 end
	store
end

getlsubtask = function(n, collect, pseudo, dρ, maxρ, configs, dists)
    nρ = length(0:dρ:maxρ)
    
    loglik = zeros(Float64, nρ)
    denom = getsum.(collect, n, pseudo)
    
    sortedconfigs, sorteddists = orderconfigs(configs, dists)

    currentconfig = [0; 0; 0]
    store = zeros(Float64, nρ)
    for i in 1:size(sortedconfigs)[1]
        # println(i)
        tmpconfig = sortedconfigs[i, :]
        if tmpconfig != currentconfig
            currentconfig = tmpconfig
            store = getstore(collect, nρ, currentconfig, pseudo)
        end

        relocate = Int.(div.((0:dρ:maxρ) .* sorteddists[i], dρ)) .+ 1
        loglik .+=  log.(store[relocate] ./ denom[relocate])
    end

    loglik
end

getl = function(n, collect0, collect1, pseudo0, pseudo1,
            dρ, maxρ, configs, dists)
    loglik0 = getlsubtask(n, collect0, pseudo0, dρ, maxρ, configs, dists)
    loglik1 = getlsubtask(n, collect1, pseudo1, dρ, maxρ, configs, dists)
    (loglik0, loglik1)
end

#= getl = function(ρs, n, results, configs, dists, pseudo)
    # ρs = []
    # for i in 1:10 push!(ρs, [0.0; i/10]) end
    l = length(ρs)
    loglik = zeros(Float64, l)
    denom = zeros(Float64, l)
    for i in 1:l denom[i] = getsum(results[i], n, pseudo) end
    for i in 1:l
        printlnln(i)
        for j in 1:size(configs)[1]
            config1 = maximum(configs[j, 1:2])
            config2 = minimum(configs[j, 1:2])

            a = [(ρ[1] - dists[j] * ρs[i][1]).^2 for ρ in ρs]
            b = [(ρ[2] - dists[j] * ρs[i][2]).^2 for ρ in ρs]
            ρidx = (1:length(ρs))[(a .== minimum(a)) .& (b .== minimum(b))][1]

            extract = results[ρidx][config1, config2]
            numer = extract[2][extract[1] .== configs[j, 3]][1] + pseudo
            if config1 != config2 numer *= 2 end
            qc = log(numer / denom[ρidx])
            if !isinf(qc) loglik[i] += qc 
            else println("uh oh") end
        end
    end

    loglik
end =#

getpseudo = function(collect, n)
    pseudo = Inf
    for i in eachindex(collect)
        for j in 1:(n - 1), k in 1:j
            tmpvec = collect[i][j, k][2]
            nonzero = tmpvec[tmpvec .> 0]
            if length(nonzero) == 0 continue end
            tmpmin = minimum(nonzero)
            if tmpmin < pseudo
                pseudo = tmpmin
            end
        end
    end
    
    pseudo / 10
end

end
