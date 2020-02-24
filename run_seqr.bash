#!/bin/env/ bash

# SCRIPT TO RUN STARSEQR IN LOCAL ENVIRONMENT

pyenv activate strubi

BASE_DIR="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/"
DATA_DIR="${BASE_DIR}/data/"
GENOME_LIB="${BASE_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/"
LEFT_FILES="${DATA_DIR}/fastq/*_1.fq*"
RIGHT_FILES="${DATA_DIR}/fastq/*_2.fq*"

# build index with STAR first
# RUN ON CLUSTER - 31GB not enough to write index to file
function make_index(){
STAR --runMode genomeGenerate --genomeFastaFiles /home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa --genomeDir /home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ --sjdbGTFfile /home/rimichael/Uni/KU_BioInf/projects/gene_fusion/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf --runThreadN 6 --sjdbOverhang 150 --genomeSAsparseD 1
}

function star_seqr(){
starseqr.py -1 ${left} -2 ${right} -m 1 -p RNA_test -t 12 -i ${sample}/STAR_INDEX -g ${GENOME_LIB}/gencode.gtf -r ${GENOME_LIB}/hg19.fa -vv
}

for left in ${LEFT_FILES}; do
   for right in ${RIGHT_FILES}; do
      if [[ ${left%_*} == ${right%_*} ]]; then
       SAMPLE_OUTPUT_DIR=${OUTPUT_DIR}/${BASENAME%_*}/star_seqr
        # check if output directory exists and create
         if [[ ! -d ${SAMPLE_OUTPUT_DIR} ]]; then
                echo "creating dir ${SAMPLE_OUTPUT_DIR}"
                mkdir ${SAMPLE_OUTPUT_DIR}
         fi
         LEFT_FASTQ=${left}
         RIGHT_FASTQ=${right}
         echo "Run SEQR on ${BASENAME} -> output: ${SAMPLE_OUTPUT_DIR}"
         echo "align... ${LEFT_FASTQ} and ${RIGHT_FASTQ}"
         star_seqr
      fi
   done
done
