ml julia

sbatch \
    --nodes=6 \
    --ntasks-per-node=13 \
    --mem=20G \
    --time=12:30:00 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --error=/scratch/users/jgottf/cocci/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    paralleldata.jl
