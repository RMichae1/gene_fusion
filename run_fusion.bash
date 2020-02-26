# SCRIPT TO RUN STAR FUSION ON THE LOCAL PROJECT STRUCTURE

BASE_DIR="/home/rimichael/Uni/KU_BioInf/projects/gene_fusion/"
GENOME_LIB="${BASE_DIR}/ref/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/"
STAR_DIRS=$(ls ${BASE_DIR}/data/)

function star_fusion(){
STAR-Fusion --genome_lib_dir ${GENOME_LIB} -J ${BASE_DIR}/data/${dir}/Chimeric.out.junction --output_dir ${OUTPUT_DIR}
}

for dir in ${STAR_DIRS}; do
   echo "Running STAR-Fusion on ${dir}"
   OUTPUT_DIR=${BASE_DIR}/output/${dir}/star_fusion
   if [[ ! -d ${OUTPUT_DIR} ]]; then
        mkdir ${OUTPUT_DIR}
   fi
   star_fusion
done
