#!/bin/env bash

# LOCAL SCRIPT TO RUN ARRIBA FUSION CALLING TOOL

#module load tools
#module load gcc
#module load star/2.7.2b
#module load arriba/1.1.0

STAR_INDEX_PATH="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/"
FASTQ1_BRAIN="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/data/fastq/sim_brain_1.fq.gz"
FASTQ2_BRAIN="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/data/fastq/sim_brain_2.fq.gz"
ASSEMBLY_FA="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"
ANNOT_GTF="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf"
BLACKLIST="/home/rimichael/System/arriba_v1.2.0/database/blacklist_hg38_GRCh38_2018-11-04.tsv.gz"

function arriba(){
# This is the default configuration as specified in the docs
STAR \
    --runThreadN 4 \
    --genomeDir ${STAR_INDEX_PATH} --genomeLoad NoSharedMemory \
    --readFilesIn ${FASTQ1_BRAIN} ${FASTQ2_BRAIN} --readFilesCommand zcat \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
    --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
    --chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 |
arriba \
    -x /dev/stdin \
    -o fusions.tsv -O fusions.discarded.tsv \
    -a ${ASSEMBLY_FA} -g ${ANNOT_GTF} -b ${BLACKLIST} \
    -T -P
}

echo "Running Arriba..."
arriba
