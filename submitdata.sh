ml julia

sbatch \
    --nodes=3 \
    --ntasks-per-node=5 \
    --mem=8G \
    --time=60 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    paralleldata.jl
