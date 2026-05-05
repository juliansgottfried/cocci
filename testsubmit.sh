ml julia

sbatch \
    --nodes=1 \
    --ntasks-per-node=2 \
    --mem=4G \
    --time=10 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    testparallel.jl
