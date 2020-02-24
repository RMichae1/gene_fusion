# gene_fusion
Gene Fusion Detection Tooling Benchmarking Project.

## File Structure
run_*.bash - Run files for tools in the benchmarking process

sample_ids.txt - List of samples used for benchmarking (publicly available through GEO browser or SRA)

## STAR
Used for alignment with CTAT genome provided by the Broad Institute (GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play)

## STAR-Fusion
Fusion caller which relies on STAR output.

Version used 1.8.1

## Arriba
Fusion Caller tool which relies on STAR output.

Version used 1.1.0

## STAR-SEQR

Using conda for starseqr greatly simplifies the installation process and required dependencies.

starseqr_requirements.txt - contains requirements for the conda environment
