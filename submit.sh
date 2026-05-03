ml julia

sbatch \
    --nodes=5 \
    --ntasks-per-node=10 \
    --mem=8G \
    --time=60 \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    parallel.jl
