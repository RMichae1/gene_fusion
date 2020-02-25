# Gene Fusion Benchmarking
Gene Fusion Detection Tooling Benchmarking Project.

![](https://encrypted-tbn0.gstatic.com/images?q=tbn%3AANd9GcS22lqiVjXS_ZcXQ25RVTVgC6YVgdR-W4Np3fJBtl277mEFauVp)

## File Structure
+ run_*.bash - Run files for tools in the benchmarking process
+ sample_ids.txt - List of samples used for benchmarking (publicly available through GEO browser or SRA)
+ ${BASE_DIR}/data/real contains real data samples from SRA as specified in sample_ids.txt
+ ${BASE_DIR}/data/synth contains synthetic benchmark data
+ ${BASE_DIR}/data/ref contains reference file (the compiled GRCh38_gencode_v31_CTAT_lib)
+ ${BASE_DIR}/data/ref/ref_genome.gem reference genome in .gem format created with gemtools (used by ChimPipe)

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

## FusionMap

FusionMap is a fusion caller tool and part of the oshell suite. This is a Windows.exe tool, which can be run in the mono environment on Linux.
This is the tool which the others are compared to.

FusionMap_options.txt - contains the specified options for the run of fusionmap

FusionMap.sh - Boilerplate Code that creates a config and command to run FusionMap from specified parameters


## ChimPipe

Requires a genome index in .gem format (use gemtools for that).
Runs with either fastq files or bam files.

## Fusion Catcher

Fusion Calling tool for somatic fusion genes, translocations and chimeras in RNAseq data
Takes fastq files and has capability to generate reference directory with `fusioncatcher-build`
