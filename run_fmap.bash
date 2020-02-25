!/bin/env bash

module load tools
module load oshell/10.0.1.69

BASE_DIR="/home/projects/cu_10160/people/ricmic/"
FMAP_SCRIPT="${BASE_DIR}/gene_fusion/FusionMap.sh"
FMAP_CONFIG="${BASE_DIR}/gene_fusion/FusionMap_options.txt"
OUTPUT_DIR="${BASE_DIR}/output/"

function run_fmap() {
    ${FMAP_SCRIPT} -c ${FMAP_CONFIG} -o ${OUTPUT_DIR} -t 4
}

echo "Running FusionMap from ${FMAP_CONFIG} ..."
echo "Writing output to ${OUTPUT_DIR}"

run_fmap

exit 0
