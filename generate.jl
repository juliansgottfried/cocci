module generate

using Intervals # https://invenia.github.io/Intervals.jl/latest/#API-1
import StatsBase
import Distributions

# include("estimate.jl")
include("/scratch/users/jgottf/cocci/estimate.jl")

initialize = function(n, l1)
    intervalsn = [[] for _ in 1:l1]
    for i in 1:n push!(intervalsn[i], IntervalSet([0..1])) end

    leavesn = [[] for _ in 1:l1]
    for i in 1:n
        filler = falses(n)
        filler[i] = 1
        push!(leavesn[i], filler)
    end
  
    ancestor = falses(l1)
    ancestor[1:n] .= 1

    timen = zeros(Float64, l1)
    
    intervalsb = [[] for _ in 1:l1]

    leavesb = [[] for _ in 1:l1]
    
    timeb = zeros(Float64, l1)

    (intervalsn, leavesn, ancestor, timen, intervalsb, leavesb, timeb)
end

expand = function(intervalsn, leavesn, ancestor, timen, intervalsb, leavesb, timeb, l1)
    newintervalsn = [[] for _ in 1:2l1]
    newintervalsn[1:l1] = intervalsn

    newleavesn = [[] for _ in 1:2l1]
    newleavesn[1:l1] = leavesn

    newancestor = falses(2l1)
    newancestor[1:l1] = ancestor

    newtimen = zeros(Float64, 2l1)
    newtimen[1:l1] = timen

    newintervalsb = [[] for _ in 1:2l1]
    newintervalsb[1:l1] = intervalsb

    newleavesb = [[] for _ in 1:2l1]
    newleavesb[1:l1] = leavesb

    newtimeb = zeros(Float64, 2l1)
    newtimeb[1:l1] = timeb

    (newintervalsn, newleavesn, newancestor, newtimen, newintervalsb, newleavesb, newtimeb, 2l1)
end

simulator = function(n, l1, ρ)
    intervalsn, leavesn, ancestor, timen, intervalsb, leavesb, timeb = initialize(n, l1)
    
    numbern = n
    numberb = 0
    time = 0
    K = n
    while K > 1
        # println(K)

        λ = (K * (K - 1 + ρ)) / 2
        time -= log(1 - rand()) / λ
        
        if rand() > ρ / (K - 1 + ρ)
            chosen = StatsBase.sample((1:l1)[ancestor], 2, replace = false)
            ancestor[chosen] .= 0
            numbern += 1
            numberb += 2
            
            if numbern > l1 || numberb > l1
                intervalsn, leavesn, ancestor, timen, intervalsb, leavesb, timeb, l1 = expand(intervalsn, leavesn, ancestor, timen, intervalsb, leavesb, timeb, l1)
            end
            
            ancestor[numbern] = 1
            timen[numbern] = time
            
            blacklen = size(intervalsn[chosen[1]])[1]
            whitelen = size(intervalsn[chosen[2]])[1]

            tmpintervals = []
            tmpleaves = []

            for i in 1:blacklen, j in 1:whitelen
                blackset = intervalsn[chosen[1]][i]
                whiteset = intervalsn[chosen[2]][j]
                intersection = intersect(blackset, whiteset)
                if !isempty(intersection)
                    push!(tmpintervals, intersection)
                    newlist = leavesn[chosen[1]][i] .| leavesn[chosen[2]][j]
                    push!(tmpleaves, newlist)
                end
            end

            interlen = size(tmpintervals)[1]

            # println("# black intervals: ", blacklen)
            # println("# white intervals: ", whitelen)
            # println("# intersections: ", interlen)

            for i in 1:blacklen
                tmpblack = copy(intervalsn[chosen[1]][i])
                for j in 1:interlen
                    interset = tmpintervals[j]
                    tmpblack = setdiff(tmpblack, interset)
                end
                if !isempty(tmpblack)
                    push!(tmpintervals, tmpblack)
                    push!(tmpleaves, leavesn[chosen[1]][i])
                end
            end

            for i in 1:whitelen
                tmpwhite = copy(intervalsn[chosen[2]][i])
                for j in 1:interlen
                    interset = tmpintervals[j]
                    tmpwhite = setdiff(tmpwhite, interset)
                end
                if !isempty(tmpwhite)
                    push!(tmpintervals, tmpwhite)
                    push!(tmpleaves, leavesn[chosen[2]][i])
                end
            end

            # check for disjointedness
            tmplen = size(tmpintervals)[1]
            counter = 0
            for i in 2:tmplen, j in 1:(i - 1)
                counter += !isempty(intersect(tmpintervals[i], tmpintervals[j]))
            end
            # println("untoward intersections: ", counter)

            tmplen = size(tmpintervals)[1]
            for i in 1:tmplen
                newlen = size(intervalsn[numbern])[1]
                merged = false
                for j in 1:newlen
                    if tmpleaves[i] == leavesn[numbern][j]
                        intervalsn[numbern][j] = union(tmpintervals[i], intervalsn[numbern][j])
                        merged = true
                        break
                    end
                end
                if !merged
                    push!(intervalsn[numbern], tmpintervals[i]) 
                    push!(leavesn[numbern], tmpleaves[i]) 
                end
            end

            # println("# new intervals: ", size(intervalsn[numbern])[1])

            timeb[[numberb - 1; numberb]] = timen[numbern] .- timen[chosen]
            intervalsb[[numberb - 1; numberb]] = intervalsn[chosen]
            leavesb[[numberb - 1; numberb]] = leavesn[chosen]

            # check for discrepancy
            counter = 0
            for x in intervalsb
                nint = size(x)[1]
                for j in 2:nint, k in 1:(j - 1)
                    counter += x[j] == x[k]
                end
            end
            # println("untoward matches: ", counter)

        else
            chosen = StatsBase.sample((1:l1)[ancestor])
            brkpt = rand()
            
            len = size(intervalsn[chosen])[1]
            blacksets = []
            whitesets = []
            blackleaves = []
            whiteleaves = []
            blacktally = 0
            whitetally = 0
            
            for i in 1:len
                tmpblack = intersect(intervalsn[chosen][i], 0..brkpt)
                tmpwhite = intersect(intervalsn[chosen][i], brkpt..1)
                if !isempty(tmpblack)
                    push!(blacksets, tmpblack)
                    push!(blackleaves, leavesn[chosen][i])
                    blacktally += sum(leavesn[chosen][i]) != n
                end
                if !isempty(tmpwhite)
                    push!(whitesets, tmpwhite)
                    push!(whiteleaves, leavesn[chosen][i])
                    whitetally += sum(leavesn[chosen][i]) != n
                end
            end
            
            if blacktally == 0 || whitetally == 0
                continue
            end

            ancestor[chosen] = 0

            numbern += 2
            numberb += 1
            
            if numbern > l1 || numberb > l1
                intervalsn, leavesn, ancestor, timen, intervalsb, leavesb, timeb, l1 = expand(intervalsn, leavesn, ancestor, timen, intervalsb, leavesb, timeb, l1)
            end

            ancestor[[numbern - 1; numbern]] .= 1
            timen[[numbern - 1; numbern]] .= time
        
            intervalsn[[numbern - 1; numbern]] = [blacksets, whitesets]
            leavesn[[numbern - 1; numbern]] = [blackleaves, whiteleaves]

            timeb[numberb] = time - timen[chosen]

            intervalsb[numberb] = intervalsn[chosen]
            leavesb[numberb] = leavesn[chosen]
        end

        K = sum(ancestor)
    end

    # println(K)

    (intervalsb, leavesb, timeb, numberb)
end

generatemut = function(n, l1, ρ, θ)
    intervalsb, leavesb, timeb, numberb = simulator(n, l1, ρ)

    mutants = [[] for _ in 1:n]
    allmutations = []
    for i in 1:numberb
        mutations = Distributions.rand(Distributions.Poisson(θ * timeb[i]))
        locs = rand(mutations)
        for j in 1:mutations push!(allmutations, locs[j]) end
        locintervals = [IntervalSet([locs[j]..locs[j]]) for j in 1:mutations]
        len = size(intervalsb[i])[1]
        for j in 1:len
            for k in 1:mutations
                intersection = intersect(intervalsb[i][j], locintervals[k])
                if !isempty(intersection)
                    push!.(mutants[leavesb[i][j]], locs[k])
                end
            end
        end
    end

    nmut = size(allmutations)[1]
    contains = falses(nmut, n)
    for j in 1:n
        if length(mutants[j]) == 0 continue end
        for i in 1:nmut
            contains[i, j] = sum(allmutations[i] .== mutants[j])
        end
    end

    (contains, allmutations)
end

getconfig = function(input)
    config = sum(input, dims = 2)[:, 1]
    push!(config, sum(input[1, :] .& input[2, :]))
end

getconfigs = function(n, l1, ρ, θ)
    contains, allmutations = generatemut(n, l1, ρ, θ)

    loci = (1:size(allmutations)[1])[0 .< sum(contains, dims = 2) .< n]
    subset = contains[loci, :]

    nloci = size(loci)[1]
    configs = zeros(Int, Int(nloci * (nloci - 1) / 2), 3)
    dists = zeros(Float64, size(configs)[1])

    for i in 2:nloci
        for j in 1:(i - 1)
            c = Int((i * (i - 3)) / 2 + 1 + j)
            configs[c, :] = getconfig(subset[[i; j], :])
            dists[c] = abs(allmutations[loci[i]] - allmutations[loci[j]])
        end
    end

    (configs, dists)
end

repeated = function(ρs, collect, pseudo, n, l1, ρ, θ, J)
    ρhat = zeros(Float64, J)
    for j in 1:J
        println("sample $(j)")
        configs, dists = generate.getconfigs(n, l1, ρ, θ)
        while length(configs) == 0 
            configs, dists = generate.getconfigs(n, l1, ρ, θ) 
        end
        loglik = estimate.getl(ρs, n, collect, configs, dists, pseudo)
        ρhat[j] = ρs[argmax(loglik)]
    end
    
    ρhat
end

end