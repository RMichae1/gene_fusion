# SCRIPT TO RUN STARSEQR IN LOCAL ENVIRONMENT

module load tools
module load anaconda3/4.4.0

source activate starseqr

BASE_DIR="/home/projects/cu_10160/people/ricmic/"
DATA_DIR="${BASE_DIR}/data/"
GENOME_LIB="${DATA_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/"
LEFT_FILES="${DATA_DIR}/synth/sim*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_1.fq*"
RIGHT_FILES="${DATA_DIR}/synth/sim*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_2.fq*"
OUTPUT_DIR="${BASE_DIR}/output/"

# build index with STAR first

function star_seqr(){
starseqr.py -1 ${LEFT_FASTQ} -2 ${RIGHT_FASTQ} -m 1 -p ${BASENAME}_SEQR -t 8 -i ${GENOME_LIB}/ref_genome.fa.star.idx/ -g ${GENOME_LIB}/ref_annot.gtf -r ${GENOME_LIB}/ref_genome.fa -vv
}


for left in ${LEFT_FILES}; do
   for right in ${RIGHT_FILES}; do
      if [[ ${left%_*} == ${right%_*} ]]; then
         BASENAME=$(basename ${left})

         LEFT_FASTQ=${left}
         RIGHT_FASTQ=${right}
         echo "Run SEQR on ${BASENAME} -> output: ${SAMPLE_OUTPUT_DIR}"
         echo "align... ${LEFT_FASTQ} and ${RIGHT_FASTQ} ..."
         star_seqr
      fi
   done
done

exit 0
