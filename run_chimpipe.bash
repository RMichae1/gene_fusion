!#/bin/env bash

# SCRIPT TO RUN chimpipe

module load tools
moudle load perl/5.24.0
module load bedtools/2.28.0
module load samtools/1.9
module load ncbi-blast/2.8.1+
module load chimpipe/0.9.5

BASE_DIR="/home/projects/cu_10160/people/ricmic/"
DATA_DIR="${BASE_DIR}/data/"
LEFT_FILES="${DATA_DIR}/synth/sim50/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_1.fq*"
GENOME_LIB="${DATA_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/"

# this tool requires a .gem index of the reference genome used - use gem-tools for that (index)

function run_chimpipe() {
    ChimPipe.sh --bam ${BAM_STAR_FILE} --genome-index ${DATA_DIR}/ref/ref_genome.gem --annotation ${GENOME_LIB}/ref_annot.gtf --sample-id ${BASENAME} -o ${SAMPLE_OUTPUT_DIR} -t 8
}

for left in ${LEFT_FILES}; do
       BASENAME=$(basename ${left})
       BAM_STAR_FILE=${BASE_DIR}/output/${BASENAME%_*}/Aligned.out.bam
       SAMPLE_OUTPUT_DIR=${OUTPUT_DIR}/${BASENAME%_*}/chimpipe
        # check if output directory exists and create
        if [[ ! -d ${SAMPLE_OUTPUT_DIR} ]]; then
            echo "creating dir ${SAMPLE_OUTPUT_DIR}"
            mkdir ${SAMPLE_OUTPUT_DIR}
        fi
        run_chimpipe
done

exit 0
