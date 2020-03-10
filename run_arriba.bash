# CLUSTER SCRIPT TO RUN ARRIBA FUSION CALLING TOOL

module load tools
module load gcc
module load star/2.7.2b
module load arriba/1.1.0

BASE_DIR="/home/projects/cu_10160/people/ricmic/data/"
SYNTH_DATA="${BASE_DIR}/synth/sim50/data.broadinstitute.org/Trinity/CTAT_FUSIONTRANS_BENCHMARKING/on_simulated_data/sim_50/reads/"
STAR_INDEX_PATH="${BASE_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/"
FASTQ1_BRAIN="${SYNTH_DATA}/sim_brain_1.fq.gz"
FASTQ2_BRAIN="${SYNTH_DATA}/sim_brain_2.fq.gz"
ASSEMBLY_FA="${BASE_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"
ANNOT_GTF="${BASE_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
BLACKLIST="${BASE_DIR}/blacklist_hg38_GRCh38_2018-11-04.tsv.gz"

function arriba(){
# This is the default configuration as specified in the docs
STAR \
    --runThreadN 4 \
    --genomeDir ${STAR_INDEX_PATH} --genomeLoad NoSharedMemory \
    --readFilesIn ${FASTQ1_BRAIN} ${FASTQ2_BRAIN} --readFilesCommand zcat \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
    --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
    --chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimOutJunctionFormat 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 |
arriba \
    -x /dev/stdin \
    -o fusions.tsv -O fusions.discarded.tsv \
    -a ${ASSEMBLY_FA} -g ${ANNOT_GTF} -b ${BLACKLIST} \
    -T -P
}

echo "Changing into output directory..."
cd /home/projects/cu_10160/people/ricmic/output
echo "Running Arriba..."

arriba
