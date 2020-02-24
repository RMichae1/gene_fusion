#!/bin/env bash

# load required modules
module load tools
module load SRA-Toolkit

# Global variables - default values for script
OUTPUT_PATH="/home/projects/cu_10160/people/ricmic/data/real"
OUTPUT_DIR="fastq"

print_help_exit () {
    echo "Options not set correctly"
    echo "---"
    echo "Option -h for this help menu."
    echo "Option -i for input txt file of sample names."
    echo "Option -o for output path for dump. DEFAULT: ${OUTPUT_PATH}"
    exit 1
}

dump_id () {
    # check if id in right format
fastq-dump --outdir ${OUTPUT_PATH}/${OUTPUT_DIR} --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-files --clip ${SAMPLEITEM} &
  printf 'Dumping... %s\n' "${SAMPLEITEM}"
}

while getopts "i:o:h:" opt; do
    case $opt in
	i) INPUT_FILE=$OPTARG  ;;
	o) OUTPUT_PATH=$OPTARG ;;
	h) print_help_exit     ;;
	*) print_help_exit >&2
    esac
done

# check for required flags
if [ ! -f ${INPUT_FILE} ]; then
    echo "Error - Input file does not exist!"
    print_help_exit
fi

echo "TEST ${OUTPUT_PATH}/${OUTPUT_DIR}"
if [ ! -d ${OUTPUT_PATH}/${OUTPUT_DIR} ]; then
   echo "Output path not found..."
   echo "Creating ${OUTPUT_PATH}/${OUTPUT_DIR}"
   mkdir ${OUTPUT_PATH}/${OUTPUT_DIR}
fi


# parse input file and use lines as arguments for dump function
while IFS="" read -r input || [ -n "${input}" ]
do
    SAMPLEITEM=${input}
    dump_id 
done < ${INPUT_FILE}
