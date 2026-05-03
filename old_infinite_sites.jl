# ok. initialize n nodes with one interval each, [0,1]. Set the interval leaves to themselves.
# if coalescence, choose two ancestor nodes. take intersection of intervals. how do we preserve the leaf info efficiently?
# for each interval, you check with every other interval. each intersection creates a new interval, adding both together, and
# each intersection you remove from the interval. at the end, you have new intersection intervals and old exclusion intervals,
# which each have the old leaf lists
# if any interval has n leaves, remove it.

# keep track of the lengths of all intervals.
# recombination is scaled by total length. when recombination occurs, select an interval weighted by length,
# and then choose a point within interval. the node that contains that interval is then split into two, with one
# node containing all intervals to the left (and the left split of the chosen interval) and the second node containing all
# to the right (and the right split), with leaves preserved.



# maintain each interval separately
# when you have a coalescent event

# each node has an array of intervals, represented by two points
# sorted by left point
# for each interval in A, you compare with intervals in B
# for each interval in A:
# scan B intervals until you find one with B left interval less than A right interval
# and and B right interval greater than A left interval
# take the intersection: max(left intervals):min(right intervals). set that interval to have all leaves
# take setdiff by taking min(left intervals):max(left intervals) and also
# min(right intervals):max(right intervals) and checking A or B for both by checking
# left for the first and right for the second, and assigning accordingly



# and we assign mutations by keeping track of branch lengths, and
# placing a poisson number of mutations on branches, each of which is at a different
# random place between 0 and 1, inherited by all leaves



# so we need an array with every node, and a list of intervals for each node,
# and left and right points for each interval
# and we need the same as left and right points, but a leaf list instead

bounds = zeros(Float64, number_of_nodes, number_of_intervals, 2)
leaves = zeros(Float64, number_of_nodes, number_of_intervals, number_of_leaves)

# so in fact it's the exact same simulation except instead of having
# black and white lists we have one list of intervals

# not all that much more complicated
# just more things to initialize and keep track of

# when coalescence happens, intervals are joined.
# when recombination happens, intervals are split along a point
# the branch for coalescent inherits intervals
# ok so branches need all the same things

# when we do the intersection, we combine leaf lists - and we do so by the same method
# hm. we can't just delete an interval. we could set each to ancestral or non


n = 20 # number of leaves
l1 = 100 # number of pre-allocated nodes and branches
l2 = 10 # number of pre-allocated intervals
l3 = 100 # length of pre-allocated leaf lists
m = 100 # monte-carlo size
ρ = 1

initialize = function(n, l1, l2)
    # node interval bounds
    boundsn = zeros(Float64, l1, l2, 2)
    # leaves contain one interval (the unit interval)
    for i in 1:n boundsn[i, 1, :] = [0; 1] end

    # node leaf lists
    leavesn = zeros(Int, l1, l2, l3)
    # leaf nodes have themselves as leaves, in their unit interval
    leavesn[1:n, 1, 1] = 1:n

    # number of leaves in each node's leaf lists
    lengthn = zeros(Int, l1, l2)
    # each node contains only one leaf, in only one interval, to start
    lengthn[1:n, 1] .= 1

    # node times
    timen = zeros(Float64, l1)

    # whether each node is currently ancestral
    ancestor = falses(l1)
    # each node begin as an ancestor
    ancestor[1:n] .= 1

    # initial number of nodes
    numbern = n

    # branch interval bounds
    boundsb = zeros(Float64, l1, l2, 2)

    # branch leaf lists
    leavesb = zeros(Int, l1, l2, l3)

    # number of leaves in each branches' leaf lists
    lengthb = zeros(Int, l1, l2)

    # branch times
    timeb = zeros(Float64, l1)

    # initial number of branches
    numberb = 0

    (boundsn, leavesn, lengthn, timen, ancestor, numbern, boundsb, leavesb, lengthb, timeb, numberb)
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

simulator = function(n, l1, l2, ρ)
    boundsn, leavesn, lengthn, timen, ancestor, numbern, boundsb, leavesb, lengthb, timeb, numberb = initialize(n, l1, l2, l3)

    # begin time
    time = 0

    # begin ancestor count
    K = n

    # stop simulation when a single ancestor has been found
    while K > 1
        println(K)

        # interval widths of ancestors
        widths = boundsn[ancestor, :, 2] .- boundsn[ancestor, :, 1]

        # remove widths of intervals that have found an mrca
        widths[lengthn[ancestor, :] .== n] .= 0

        # scale ρ by the proportion of ancestral material
        ρ_eff = ρ * sum(widths) / K
      
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
            # if numbern > l1
            #     leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1 = expand_l(leavesn, lengthn, timen, ancestor, leavesb, lengthb, timeb, l1, l2)
            # end

            # new node is ancestral
            ancestor[numbern] = 1

            # set node's time
            timen[numbern] = time

            

            # now i have to compute new intervals
            lengthn[chosen, :]
            leavesn[chosen, :, :]
            boundsn[chosen, :, :]

            

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

    println(K)

    (leavesb, lengthb, timeb)
end