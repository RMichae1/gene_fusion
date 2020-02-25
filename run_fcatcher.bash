 
module load tools
module load anaconda2/4.4.0
module load jdk/13.0.1
module load openjdk/13.0.1
module load java/1.8.0
module load jre/1.8.0
module load fusioncatcher/1.20

DATA_DIR="/home/projects/cu_10160/people/ricmic/data"
OUTPUT_DIR="/home/projects/cu_10160/people/ricmic/output/"
LEFT_FILES="${DATA_DIR}/synth/*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_1.fq*"
RIGHT_FILES="${DATA_DIR}/synth/*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_2.fq*"
REF_DIR="${DATA_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/"

# FUSIONCATCHER RELIES ON PYTHON2

function run_catcher() {
    fusioncatcher -d ${REF_DIR} -i ${FASTQ_FILE} -o ${SAMPLE_OUTPUT_DIR}
}

for left in ${LEFT_FILES}; do
   for right in ${RIGHT_FILES}; do
      if [[ ${left%_*} == ${right%_*} ]]; then
        BASENAME=$(basename ${left})
        SAMPLE_OUTPUT_DIR=${OUTPUT_DIR}/${BASENAME%_*}
        # check if output directory exists and create
         if [[ ! -d ${SAMPLE_OUTPUT_DIR} ]]; then
                echo "creating dir ${SAMPLE_OUTPUT_DIR}"
                mkdir ${SAMPLE_OUTPUT_DIR}
         fi
         # Fusioncatcher also takes single FASTQ files to read
         FASTQ_FILES="${left},${right}"
         
         echo "Running FusionCatcher on ${BASENMAE}..."
         echo "With ${FASTQ_FILES}"
         
         fusioncatcher
      fi
   done
done

exit 0



