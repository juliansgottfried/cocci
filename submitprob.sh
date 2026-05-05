ml julia

sbatch \
    --nodes=5 \
    --ntasks-per-node=10 \
    --mem=4G \
    --time=120 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    parallelprob.jl
