#!/usr/bin/env bash

set -euo pipefail

#echo "Listing *.h5mu files under /data (sorted by size ascending):"
#mapfile -d '' raw_paths < <(find /data -name "*HT.h5mu" -print0)

files=(
/data/06881207-692c-4ec9-8255-5ffbafc92b4d/HT.h5mu
/data/5fa94f5d-8ca0-4352-9a2f-21b42c75ece4/SI.h5mu
/data/f79ce9ea-0a44-4374-b1e0-1979bdbfe0b0/LI.h5mu
#/data/89dc5602-a511-464b-9c67-c874a8b20db7/LF_raw.h5mu  #pass
#/data/80ebfff5-65a9-4991-8287-630b4b36891b/LO_raw.h5mu  #pass
#/data/3d150b98-f129-4a89-89e1-eb9c576ed2e0/RF_raw.h5mu  #pass
#/data/855e6696-7c54-44ea-b176-ce28bdf41ad2/UT_raw.h5mu  #pass
#/data/bd092435-ce5f-4c54-80af-fb63c2d471e9/RO_raw.h5mu  #pass
#/data/3f1906f0-f0e2-40c4-bda8-cb5f89d1d9d7/SP_raw.h5mu  #Found modalities: ['3f1906f0-f0e2-40c4-bda8-cb5f89d1d9d7_raw']
#/data/3769013e-ccad-43be-9be1-577dbcc6d600/LI_raw.h5mu  # Found modalities: ['3769013e-ccad-43be-9be1-577dbcc6d600_raw']
#/data/e80e6a14-d43b-43da-9ae6-fb33e938934a/SI_raw.h5mu  #Found modalities: ['e80e6a14-d43b-43da-9ae6-fb33e938934a_raw']
#/data/cded847a-6064-4e68-a44b-d5f5578338cd/LY_raw.h5mu #Found modalities: ['cded847a-6064-4e68-a44b-d5f5578338cd_raw']
#/data/a0d5b879-18ff-4d1e-8061-52b01d63e659/TH_raw.h5mu  #Found modalities: ['a0d5b879-18ff-4d1e-8061-52b01d63e659_raw']
)

for f in "${files[@]}"; do
    echo "Running: $f"
    #python check_h5mu_bins.py --h5mu "$f"
done


if ((${#files[@]})); then
  printf '%s\0' "${files[@]}" | xargs -0 ls -lhSr

  mapfile -t sorted_paths < <(
    printf '%s\0' "${files[@]}" \
      | xargs -0 stat --format='%s %n' \
      | LC_ALL=C sort -n -k1,1 \
      | awk '{size=$1; $1=""; print substr($0,2)}'
  )

  for h5mu in "${sorted_paths[@]}"; do
    base=$(basename "${h5mu%.h5mu}")
    echo "Processing ${h5mu} → /data/pseudo_fragments/${base}_new_pseudo"
    python /data/pseudo_fragments/ShinyCell2_multiomic/pseudo_frags/h5mu_to_pseudofragments.py \
      --h5mu "$h5mu" \
      --out "/data/pseudo_fragments/${base}_new_pseudo"
  done
else
  echo "No *raw.h5mu files found."
fi
