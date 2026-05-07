ml julia

sbatch \
    --nodes=10 \
    --ntasks-per-node=10 \
    --mem=10G \
    --time=180 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    parallelprob.jl
