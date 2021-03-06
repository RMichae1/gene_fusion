# Gene Fusion Benchmarking
Gene Fusion Detection Tooling Benchmarking Project.

![](https://encrypted-tbn0.gstatic.com/images?q=tbn%3AANd9GcS22lqiVjXS_ZcXQ25RVTVgC6YVgdR-W4Np3fJBtl277mEFauVp)

## File Structure
+ run_*.bash - Run files for tools in the benchmarking process
+ sample_ids.txt - List of samples used for benchmarking (publicly available through GEO browser or SRA)
+ ./truth contains the truth data to synthetic fusions
+ ${BASE_DIR}/data/real contains real data samples from SRA as specified in sample_ids.txt
+ ${BASE_DIR}/data/synth contains synthetic benchmark data
+ ${BASE_DIR}/data/ref contains reference file (the compiled GRCh38_gencode_v31_CTAT_lib)
+ ${BASE_DIR}/data/ref/ref_genome.gem reference genome in .gem format created with gemtools (used by ChimPipe)

## STAR
Used for alignment with CTAT genome provided by the Broad Institute (GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play)

## STAR-Fusion
Fusion caller which relies on STAR output.

Version used 1.8.1

*Status*: Ran to completion on local setup with correct Perl dependencies.

## Arriba
Fusion Caller tool which relies on STAR output.

Version used 1.1.0

*Status*: STAR aligner output of version used on the cluster compute node in conflict with Arriba
required input. Errors in formatting.

## STAR-SEQR

Using conda for starseqr greatly simplifies the installation process and required dependencies.

starseqr_requirements.txt - contains requirements for the conda environment

## FusionMap

FusionMap is a fusion caller tool and part of the oshell suite. This is a Windows.exe tool, which can be run in the mono environment on Linux.
This is the tool which the others are compared to.

FusionMap_options.txt - contains the specified options for the run of fusionmap

FusionMap.sh - Boilerplate Code that creates a config and command to run FusionMap from specified parameters

Make sure to run with correct reference file version (hg38) with right annotation.

*Status*: Ran to completion on cluster-compute node.

## ChimPipe

Requires a genome index in .gem format (use gemtools for that).
Runs with either fastq files or bam files.

*Status*: gemtools generated index faulty and key-file empty. Error could not be identified.
No pipeline run possible.

## Fusion Catcher

Fusion Calling tool for somatic fusion genes, translocations and chimeras in RNAseq data
Takes fastq files and has capability to generate reference directory with `fusioncatcher-build`.
Due to requirement conflicts (Python 2.7.X) a pre-built data directory was used.
Setup of reference directory from:
```
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.aa
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.ab
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.ac
wget http://sourceforge.net/projects/fusioncatcher/files/data/human_v98.tar.gz.ad
```

*Status*: Ran to completion on cluster-compute node

## Trinity Fusion

Trinity Fusion Pipeline is a multi-stage pipeline 
requirements:
 - python environment with pysam installed
 - perl >v.5 with Set::IntervalTree installed
 - gmap index of your reference required
 
 *Status*: Ran to completion on local set-up after using anaconda environment for first 2 stages of pipeline 
 and exiting the environment for the last stage and required perl dependencies
 
 ## TopHat Fusion
 
 Requirements:
 - samtools
 - bowtie2
 
 *Status*: Conflict between TopHat2 and bowtie. Error when reading from index after ~4-6h runtime.
