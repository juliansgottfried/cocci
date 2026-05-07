ml julia

sbatch \
    --nodes=3 \
    --ntasks-per-node=5 \
    --mem=10G \
    --time=24:00:00 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    paralleldata.jl
