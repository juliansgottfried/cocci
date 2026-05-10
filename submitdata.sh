ml julia

sbatch \
    --nodes=5 \
    --ntasks-per-node=100 \
    --mem=5G \
    --time=12:00:00 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    paralleldata.jl
