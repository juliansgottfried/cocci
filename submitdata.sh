ml julia

sbatch \
    --nodes=8 \
    --ntasks-per-node=20 \
    --mem=15G \
    --time=03:30:00 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --error=/scratch/users/jgottf/cocci/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    paralleldata.jl
