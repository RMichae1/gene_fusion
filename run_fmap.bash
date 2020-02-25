!/bin/env bash

module load tools
module load oshell/10.0.1.69

BASE_DIR="/home/projects/cu_10160/people/ricmic/"
FMAP_SCRIPT="${BASE_DIR}/gene_fusion/FusionMap.sh"
FMAP_CONFIG="${BASE_DIR}/gene_fusion/FusionMap_options.txt"
OUTPUT_DIR="${BASE_DIR}/output/"

LEFT_FILES="${DATA_DIR}/synth/*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_1.fq*"
RIGHT_FILES="${DATA_DIR}/synth/*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_2.fq*"

function run_fmap() {
    ${FMAP_SCRIPT} -c ${FMAP_CONFIG} -o ${SAMPLE_OUTPUT_DIR} -t 4 ${LEFT_FASTQ} ${RIGHT_FASTQ}
}

for left in ${LEFT_FILES}; do
   for right in ${RIGHT_FILES}; do
      if [[ ${left%_*} == ${right%_*} ]]; then
        BASENAME=$(basename ${left})
        SAMPLE_OUTPUT_DIR=${OUTPUT_DIR}/${BASENAME%_*}/fmap
        
        if [[ ! -d ${SAMPLE_OUTPUT_DIR} ]]; then
            echo "creating dir ${SAMPLE_OUTPUT_DIR}"
            mkdir ${SAMPLE_OUTPUT_DIR}
        fi
        
        LEFT_FASTQ=${left}
        RIGHT_FASTQ=${right}
        
        echo "Running FusionMap on ${BASENAME}..."
        echo "Mapping ${LEFT_FASTQ} with ${RIGHT_FASTQ}..."
        echo "Writing output to ${SAMPLE_OUTPUT_DIR}"
        run_fmap
        
      fi
    done
done

exit 0
