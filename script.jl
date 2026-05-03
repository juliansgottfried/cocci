#!/usr/bin/env julia

using Distributed, SlurmClusterManager
addprocs(SlurmManager())
@everywhere println("hello from $(myid()):$(gethostname())")
all_results = pmap(1:20) do index
	println("index: $(index), id: $(myid()), host: $(gethostname())")
end
