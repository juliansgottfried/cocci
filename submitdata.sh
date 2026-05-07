ml julia

sbatch \
    --nodes=10 \
    --ntasks-per-node=10 \
    --mem=10G \
    --time=70 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    paralleldata.jl
