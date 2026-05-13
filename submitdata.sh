ml julia

sbatch \
    --nodes=5 \
    --ntasks-per-node=10 \
    --mem=30G \
    --time=01:30:00 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    paralleldata.jl
