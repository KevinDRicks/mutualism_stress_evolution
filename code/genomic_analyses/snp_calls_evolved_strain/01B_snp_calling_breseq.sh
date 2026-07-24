#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --output=output_breseq.txt
#SBATCH --error=error_breseq.txt
#SBATCH --mem=96G

#this script is one of multiple methods used to detect introduced mutations in our evolved strains
#in this script, raw sequences from evolved strains are passed through the breseq pipeline 
#and compared against the high quality ancestor assembly that have previously been built
#raw long reads for can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1294393/
#ancestral assemblies can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA310138/

#this script must by run for each evolved strain (could be modified to run in a loop)



#load libraries
module load StdEnv/2023
module load breseq/0.38.2
module load bowtie2/2.5.2

#path to the short reads that have been passed through quality control
#will include both forward and reverse reads
#additional path to assembly of ancestor, to map back to
target_evolve1="/path/to/long_reads/[STRAIN_ID]_R1_longread_001.fastq.gz"
target_evolve2="/path/to/long_reads/[STRAIN_ID]_R2_longread_001.fastq.gz"
target_ancestor="/path/to/ancestor/assembly/anc_assembly.fasta"

#pass path to raw sequences and ancestor assembly into breseq 
breseq -j 8 -p -r "$target_ancestor" "$target_evolve1" "$target_evolve2" > log.txt