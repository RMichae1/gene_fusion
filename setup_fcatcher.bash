# Setting up fusion catcher reference directory for later processing

module load tools
module load anaconda2/4.4.0
module load jdk/13.0.1
module load openjdk/13.0.1
module load java/1.8.0
module load jre/1.8.0
module load fusioncatcher/1.20

REF_DIR="${DATA_DIR}/ref/catcher_data/" 

fusioncatcher-build --output ${REF_DIR} --organism='homo_sapiens' --threads=4
