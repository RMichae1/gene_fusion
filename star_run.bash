# STAR ALIGNER

# load dependencies
module load tools
module load perl
module load samtools/1.9
module load star/2.7.2b
module load star-fusion/1.8.1

DATA_DIR="/home/projects/cu_10160/people/ricmic/data"
OUTPUT_DIR="/home/projects/cu_10160/people/ricmic/output/"
CTAT_GENOME_LIB="/home/projects/cu_10160/people/ricmic/data/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx"
LEFT_FILES="${DATA_DIR}/synth/*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_1.fq*"
RIGHT_FILES="${DATA_DIR}/synth/*/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/*/reads/*_2.fq*"


# First run STAR aligner - this is input for other tools (e.g. Arriba)
run_star() {
STAR --genomeDir ${CTAT_GENOME_LIB} \
--readFilesIn ${LEFT_FASTQ} ${RIGHT_FASTQ} \
--runThreadN 1 \
--outReadsUnmapped None \
--chimSegmentMin 12 \ # minimum overlap fusion detection
--chimJunctionOverhangMin 12 \ # minimum required overhang fusion detection
--chimOutJunctionFormat 1 \ # **essential** includes required metadata in Chimeric.junction.out file.
--alignSJDBoverhangMin 10 \ # minimum overhang for annotated spliced alignments
--alignMatesGapMax 100000 \   # avoid readthru fusions within 100k
--alignIntronMax 100000 \
--alignSJstitchMismatchNmax 5 -1 5 5 \   # settings improved certain chimera detections
--twoPassMode Basic \
--readFilesCommand "gunzip -c" \
--outSAMstrandField intronMotif \
--outSAMunmapped Within \
--outSAMtype BAM unsorted \
--outSAMattrRGline ID:GRPundef \
--chimMultimapScoreRange 10 \ # (recomm. 3) chimeric alignments multimap score **TEST** 
--chimScoreJunctionNonGTAG -4 \ # Test if this is good for junction
--chimMultimapNmax 10 \ # (recomm. 20) max mappings
--chimNonchimScoreDropMin 10 \
--peOverlapNbasesMin 12 \ # TEST publication has NbasesM 0.1 value...
--peOverlapMMp 0.1
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
         LEFT_FASTQ=${left}
         RIGHT_FASTQ=${right}
         echo "Run Aligner on ${BASENAME} -> output: ${SAMPLE_OUTPUT_DIR}"
         echo "align... ${LEFT_FASTQ} and ${RIGHT_FASTQ}"
         run_star
      fi
   done
done
          
