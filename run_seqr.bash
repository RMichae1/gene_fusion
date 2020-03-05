# SCRIPT TO RUN STARSEQR IN LOCAL ENVIRONMENT

module load tools
module load anaconda3/4.4.0

source activate starseqr

BASE_DIR="/home/projects/cu_10160/people/ricmic/"
DATA_DIR="${BASE_DIR}/data/"
GENOME_LIB="${DATA_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/"
LEFT_FILES="${DATA_DIR}/synth/sim50/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_1.fq*"
RIGHT_FILES="${DATA_DIR}/synth/sim50/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_2.fq*"
OUTPUT_DIR="${BASE_DIR}/output/"

# build index with STAR first

function star_seqr(){
starseqr.py -1 ${LEFT_FASTQ} -2 ${RIGHT_FASTQ} -m 1 -p RNA_${BASENAME%_*} -t 4 -sb ${BAM_STAR_FILE} -sj ${JUNCT_STAR_FILE} -g ${GENOME_LIB}/ref_annot.gtf -r ${GENOME_LIB}/ref_genome.fa -vv
}

for left in ${LEFT_FILES}; do
   for right in ${RIGHT_FILES}; do
      if [[ ${left%_*} == ${right%_*} ]]; then
       BASENAME=$(basename ${left})
       SAMPLE_OUTPUT_DIR=${OUTPUT_DIR}/${BASENAME%_*}/star_seqr
       JUNCT_STAR_FILE=${OUTPUT_DIR}/${BASENAME%_*}/Chimeric.out.junction
       BAM_STAR_FILE=${OUTPUT_DIR}/${BASENAME%_*}/Aligned.out.bam

        # check if output directory exists and create
         if [[ ! -d ${SAMPLE_OUTPUT_DIR} ]]; then
                echo "creating dir ${SAMPLE_OUTPUT_DIR}"
                mkdir ${SAMPLE_OUTPUT_DIR}
         fi

	# STAR-SEQR cannot handle header or commenting tail lines => clean files
	HEADER_FIRST_ELEM="chr_donorA"
	FOOTER_ELEM="#"
	FILE_FIRST_ELEM=$(head -n 1 ${JUNCT_STAR_FILE} | awk -F '\t' '{print $1}')
	FILE_FOOTER_ELEM=$(tail -n 1 ${JUNCT_STAR_FILE} | awk -F '\t' '{print substr($1,1,1)}')
	# We are going to change the file around.. lets make a backup first
	if [[ ${HEADER_FIRST_ELEM} == ${FILE_FIRST_ELEM} ]] || [[ ${FOOTER_ELEM} == ${FILE_FOOTER_ELEM} ]]; then
		cp ${JUNCT_STAR_FILE} ${JUNCT_STAR_FILE}.bak
	fi
	# use first value to test if header exists
	 if [[ ${HEADER_FIRST_ELEM} == ${FILE_FIRST_ELEM} ]]; then
		# clean header
		sed 1,1d ${JUNCT_STAR_FILE} > ${JUNCT_STAR_FILE}	
	 fi
	 if [[ ${FOOTER_ELEM} == ${FILE_FOOTER_ELEM} ]]; then
		# clean tail commenting lines
		head -n -1 ${JUNCT_STAR_FILE} > ${JUNCT_STAR_FILE}
	 fi
         LEFT_FASTQ=${left}
         RIGHT_FASTQ=${right}
         echo "Run SEQR on ${BASENAME} -> output: ${SAMPLE_OUTPUT_DIR}"
         echo "align... ${LEFT_FASTQ} and ${RIGHT_FASTQ} ..."
         echo "using BAM: ${BAM_STAR_FILE} and Junctions: ${JUNCT_STAR_FILE}..."
         star_seqr
      fi
   done
done

exit 0
