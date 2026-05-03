module Main

using StatsBase: sample

initialize = function(n, l1, l2)
    # node black and white leaf lists
    leavesn = zeros(Int, l1, 2, l2)
    # leaf nodes have themselves as black and white leaves
    leavesn[1:n, :, 1] = [1:n 1:n]

    # number of leaves in each node's black and white list
    lengthn = zeros(Int, l1, 2)
    # each node contains only one black and white leaf, to start
    lengthn[1:n, 1:2] .= 1

    # node times
    timen = zeros(Float64, l1)

    # whether each node is currently ancestral
    ancestor = falses(l1)
    # each node begin as an ancestor
    ancestor[1:n] .= 1
    
    # initial number of nodes
    numbern = n

    # branch black and white leaf lists
    leavesb = zeros(Int, l1, 2, l2)

    # number of leaves in each branches' black and white list
    lengthb = zeros(Int, l1, 2)

    # branch times
    timeb = zeros(Float64, l1)

    # number of branches
    numberb = 0

    (leavesn, lengthn, timen, ancestor, numbern, leavesb, lengthb, timeb, numberb)
end

expand_l = function(leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1, l2)
    new_leavesn = zeros(Int, 2l1, 2, l2)
    new_leavesn[1:l1, :, :] = leavesn

    new_lengthn = zeros(Int, 2l1, 2)
    new_lengthn[1:l1, :] = lengthn

    new_timen = zeros(Float64, 2l1)
    new_timen[1:l1] = timen

    new_ancestor = falses(2l1)
    new_ancestor[1:l1] = ancestor

    new_leavesb = zeros(Int, 2l1, 2, l2)
    new_leavesb[1:l1, :, :] = leavesb

    new_lengthb = zeros(Int, 2l1, 2)
    new_lengthb[1:l1, :] = lengthb

    new_timeb = zeros(Float64, 2l1)
    new_timeb[1:l1] = timeb

    (new_leavesn, new_lengthn, new_timen, new_ancestor, new_leavesb, new_lengthb, new_timeb, 2l1)
end

expand_llist = function(leavesn, leavesb, l1, l2)
    new_leavesn = zeros(Int, l1, 2, 2l2)
    new_leavesn[:, :, 1:l2] = leavesn

    new_leavesb = zeros(Int, l1, 2, 2l2)
    new_leavesb[:, :, 1:l2] = leavesb
    
    (new_leavesn, new_leavesb, 2l2)
end

simulator = function(n, l1, l2, ρ)
    leavesn, lengthn, timen, ancestor, numbern, leavesb, lengthb, timeb, numberb = initialize(n, l1, l2)

    # begin time
    time = 0

    # begin ancestor count
    K = n

    # stop simulation when a single ancestor has been found
    while K > 1
        # println(K)

        ρ_eff = 0

        # only recombination if neither black nor white have found common ancestor
        if sum(lengthn[ancestor, 1]) > 1 && sum(lengthn[ancestor, 2]) > 1
            # scale recombination by how many ancestors have both black and white
            ρ_eff = ρ * sum(sum(eachcol(lengthn[ancestor, :] .!= 0)) .== 2) / K
        end

        # get event rate
        λ = (K * (K - 1 + ρ_eff)) / 2

        # get time of next event
        time -= log(1 - rand()) / λ

        if rand() > ρ_eff / (K - 1 + ρ_eff) # if CA:
            # sample two ancestors for coalescence
            chosen = sample((1:l1)[ancestor], 2, replace = false)

            # coalesced ancestors are no longer ancestral
            ancestor[chosen] .= 0

            # add a new node
            numbern += 1
            if numbern > l1
                leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1 = expand_l(leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1, l2)
            end

            # new node is ancestral
            ancestor[numbern] = 1

            # set node's time
            timen[numbern] = time

            # calculate merged black and white leaf list lengths
            new_length = sum(eachrow(lengthn[chosen, :]))
            if new_length[1] > l2 || new_length[2] > l2
                leavesn, leavesb, l2 = expand_llist(leavesn, leavesb, l1, l2)
            end

            # set leaf list lengths
            lengthn[numbern, :] = new_length

            # get black leaf list lengths for both sampled ancestors
            midpoint = lengthn[chosen[1], :]
            
            # merge lists (separately for black and white) and place in new node
            leavesn[numbern, :, :] = [leavesn[chosen[1], 1, 1:midpoint[1]] ; 
                leavesn[chosen[2], 1, 1:(l2 - midpoint[1])] ;;
                leavesn[chosen[1], 2, 1:midpoint[2]] ; 
                leavesn[chosen[2], 2, 1:(l2 - midpoint[2])]]'

            # add two new branches
            numberb += 2
            if numberb > l1
                leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1 = expand_l(leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1, l2)
            end

            # set branch times separately
            timeb[[numberb - 1; numberb]] = timen[numbern] .- timen[chosen]

            # place separate list lengths in each branch
            lengthb[[numberb - 1; numberb], :] = lengthn[chosen, :]

            # place separate node lists in each branch
            leavesb[[numberb - 1; numberb], :, :] = leavesn[chosen, :, :]
        
        else # if RE:
            # choose one ancestor to recombine, from those with both colors
            chosen = sample((1:l1)[ancestor .& (sum(eachcol(lengthn .!= 0)) .== 2)])

            # no longer ancestral
            ancestor[chosen] = 0

            # add two new nodes
            numbern += 2
            if numbern > l1
                leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1 = expand_l(leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1, l2)
            end

            # make nodes ancestral
            ancestor[[numbern - 1; numbern]] .= 1

            # set the same time for both
            timen[[numbern - 1; numbern]] .= time

            # split list lengths into the two nodes
            lengthn[numbern - 1, 1] = lengthn[chosen, 1]
            lengthn[numbern, 2] = lengthn[chosen, 2]

            # separate black and white leaf lists into each node
            leavesn[numbern - 1, 1, :] = leavesn[chosen, 1, :]
            leavesn[numbern, 2, :] = leavesn[chosen, 2, :]

            # add a single branch
            numberb += 1
            if numberb > l1
                leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1 = expand_l(leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1, l2)
            end

            # set branch time
            timeb[numberb] = time - timen[chosen]

            # set branch list length
            lengthb[numberb, :] = lengthn[chosen, :]

            # set branch leaf list
            leavesb[numberb, :, :] = leavesn[chosen, :, :]
        end

        # compute new ancestral count
        K = sum(ancestor)
    end

    # println(K)

    (leavesb, lengthb, timeb)
end

geth! = function(config12, config3, leavesb, lengthb, timeb)
    l = size(leavesb)[1]
    checker = falses(l, 2)
    for i in 1:l
        checker[i, :] = lengthb[i, :] .== config12
    end

    blackidx = (1:l)[checker[:, 1]]
    whiteidx = (1:l)[checker[:, 2]]
    
    for i in 1:size(config3[1])[1]
        for j in blackidx, k in whiteidx
            if sum(intersect(leavesb[j, 1, :], leavesb[k, 2, :]) .!= 0) == config3[1][i]
                config3[2][i] += timeb[j] * timeb[k]
            end
        end
    end
end

configloop = function(configs, n, leavesb, lengthb, timeb)
    for i in 1:(n - 1)
        for j in 1:(n - 1)
            geth!([i; j], configs[i, j], leavesb, lengthb, timeb)
        end
    end

    configs
end

buildconfigs = function(i, j, n)
    k = max(0, i + j - n):min(i, j)
    [k, zeros(Float64, size(k)[1])]
end

sumconfigs! = function(configs, add, n)
    for i in 1:(n - 1), j in 1:(n - 1)
        configs[i, j][2] += add[i, j][2]
    end
end

divideconfigs! = function(configs, n, m)
    for i in 1:(n - 1), j in 1:(n - 1)
        configs[i, j][2] /= m
    end
end

montecarlo = function(n, l1, l2, ρ, m)
    configs = [buildconfigs(i, j, n) for i = 1:(n - 1), j = 1:(n - 1)]
    for _ in 1:m
        leavesb, lengthb, timeb = Main.simulator(n, l1, l2, ρ)
        add = configloop(configs, n, leavesb, lengthb, timeb)
        sumconfigs!(configs, add, n)
    end
    divideconfigs!(configs, n, m)

    configs
end

rhogrid = function(n, l1, l2, ρs, m)
    results = []
    for ρ in ρs
        push!(results, montecarlo(n, l1, l2, ρ, m))
    end
    results
end

getsum = function(configs, n, pseudo)
    summer = 0
    for i in 1:(n - 1), j in 1:(n - 1)
        summer += sum(configs[i, j][2] .+ pseudo)
    end
    summer
end

getqc = function(results, ρidx, n, config, pseudo)
    tmpresult = results[ρidx][config[1], config[2]]
    tmpresult[2][findfirst(tmpresult[1] .== config[3])] + pseudo
end

getl = function(ρs, n, results, configs, dists, pseudo)
    l = size(ρs)[1]
    loglik = ones(Float64, l)
    denoms = similar(loglik)
    for i in 1:l denoms[i] = getsum(results[l], n, pseudo) end
    for i in 1:l
        for j in 1:size(configs)[1]
            ρidx = argmin(((dists[j] * ρs[i]) .- ρs) .^ 2)
            qc = log(getqc(results, ρidx, n, configs[j, :], pseudo) / denoms[ρidx])
            if !isinf(qc) loglik[i] += qc end
        end
    end
    loglik
end

end

using .Main
