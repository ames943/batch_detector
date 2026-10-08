#!/bin/bash
# Rebuild the colorectal cancer (Duvallet et al. 2017, microbiomeHD) input
# tables that scripts/revision/task1_* and task2_* read.
#
# Run from the paper/ folder:   bash data/get_duvallet_crc.sh
#
# The four study archives come from the microbiomeHD Zenodo record
# (10.5281/zenodo.840333, CC BY-NC 4.0). They are not stored in this repo.
# Output: results/batch_detector_external_validation/input_feature_matrix.tsv
#         results/batch_detector_external_validation/input_labels.tsv

set -e
mkdir -p microbiomeHD_crc
for study in zeller baxter zackular zhao; do
    curl -L -o "microbiomeHD_crc/${study}.tar.gz" \
        "https://zenodo.org/api/records/840333/files/crc_${study}_results.tar.gz/content"
    tar xzf "microbiomeHD_crc/${study}.tar.gz" -C microbiomeHD_crc
    rm "microbiomeHD_crc/${study}.tar.gz"
done

python scripts/prepare_duvallet_input.py
