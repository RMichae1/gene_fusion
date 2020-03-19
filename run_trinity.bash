# SCRIPT TO RUN TRINTIY FUSION

module load tools
module load gmap/20180325
module load bowtie2/2.3.4.1
module load samtools/1.9
module load jellyfish/2.2.7
module load salmon/1.1.0
module load java/1.8.0
module load jre/1.8.0
module load openjdk/13.0.1
module load jdk/13.0.1
module load trinityrnaseq/2.9.1
module load trinityfusion/0.3.4

# pysam needs to be loaded via anaconda

BASE_DIR="/home/projects/cu_10160/people/ricmic/"
DATA_DIR="${BASE_DIR}/data/"
export CTAT_GENOME_LIB="${DATA_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/"
LEFT_FILES="${DATA_DIR}/synth/*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_1.fq*"
RIGHT_FILES="${DATA_DIR}/synth/*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_2.fq*"
OUTPUT_DIR="${BASE_DIR}/output/"

# reference CTAT genome index with GMAP
# build index with STAR -> generate fusion junctions and alignment bam w/ fusions enabled

function run_trinity(){
TrinityFusion --left_fq ${LEFT_FASTQ} --right_fq ${RIGHT_FASTQ} --CPU 4 --aligned_bam ${BAM_STAR_FILE} --chimeric_junctions ${JUNCT_STAR_FILE} --genome_lib_dir ${GENOME_LIB} --output_dir ${SAMPLE_OUTPUT_DIR}
}

for left in ${LEFT_FILES}; do
   for right in ${RIGHT_FILES}; do
      if [[ ${left%_*} == ${right%_*} ]]; then
       BASENAME=$(basename ${left})
       SAMPLE_OUTPUT_DIR=${OUTPUT_DIR}/${BASENAME%_*}/trinity
       JUNCT_STAR_FILE=${OUTPUT_DIR}/${BASENAME%_*}/Chimeric.out.junction
       BAM_STAR_FILE=${OUTPUT_DIR}/${BASENAME%_*}/Aligned.out.bam

        # check if output directory exists and create
         if [[ ! -d ${SAMPLE_OUTPUT_DIR} ]]; then
                echo "creating dir ${SAMPLE_OUTPUT_DIR}"
                mkdir ${SAMPLE_OUTPUT_DIR}
         fi

         LEFT_FASTQ=${left}
         RIGHT_FASTQ=${right}
         echo "Run TRINITY on ${BASENAME} -> output: ${SAMPLE_OUTPUT_DIR}"
         echo "align... ${LEFT_FASTQ} and ${RIGHT_FASTQ} ..."
         echo "using BAM: ${BAM_STAR_FILE} and Junctions: ${JUNCT_STAR_FILE}..."
         run_trinity
      fi
   done
done

exit 0
