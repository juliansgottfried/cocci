ml julia

sbatch \
    --nodes=5 \
    --ntasks-per-node=15 \
    --mem=50G \
    --time=15:00:00 \
    --output=/scratch/users/jgottf/cocci/output/%j.out \
    --error=/scratch/users/jgottf/cocci/output/%j.out \
    --partition=normal,hns \
    --mail-type=ALL \
    --mail-user=juliansgottfried@gmail.com \
    paralleldata.jl
