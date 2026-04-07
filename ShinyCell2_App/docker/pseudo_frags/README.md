docker build -t h5mu2frag .
docker run --rm \
  -v "$(pwd)/data":/app/data \
  h5mu2frag:latest \
  python /app/h5mu_to_pseudofragments.py \
  --h5mu /app/data/HT.h5mu --out /app/data/HT_pseudo.tsv


#!/bin/bash

# Base directory where the data is mounted
BASE_DIR="/hive/hubmap/data/public/hubmap-data-products"

# List of input files
files=(
"./35de4ef6-7dfc-49cf-bb8c-3573a3d768a6/HT_processed.h5mu"
"./89dc5602-a511-464b-9c67-c874a8b20db7/LF_processed.h5mu"
"./80ebfff5-65a9-4991-8287-630b4b36891b/LO_processed.h5mu"
"./3d150b98-f129-4a89-89e1-eb9c576ed2e0/RF_processed.h5mu"
)

# Loop through each file and run the Docker command
for f in "${files[@]}"; do
    fname=$(basename "$f" .h5mu)
    docker run --rm \
        -v "${BASE_DIR}":/app/data \
        h5mu2frag:latest \
        python /app/h5mu_to_pseudofragments.py \
        --h5mu /app/data/"$f" \
        --out /app/data/"${fname}_pseudo.tsv"
done


# SLURM (interactive) example

# Start an interactive session on WORKSPACES-CPU, node l004
srun -p WORKSPACES-CPU --nodelist=l004 --time=06:00:00 --cpus-per-task=4 --mem=128G --pty bash

# Then run the converter inside the allocated node
docker run --rm \
    -v /hive/hubmap/data/public/hubmap-data-products:/data \
    h5mu_frag \
    bash -c "bash /data/pseudo_fragments/ShinyCell2_multiomic/pseudo_frags/run_parallel_psedu.sh"
