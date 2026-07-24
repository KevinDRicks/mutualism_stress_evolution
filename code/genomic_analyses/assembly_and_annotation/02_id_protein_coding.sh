#!/bin/bash
#code passes the reference ancestral rhizobia genomes through prokka for annotation
#primary utility here is for an easy passage through prodigal for identifying protein-coding genes 
#note: this script was initially run on Compute Canada server
#note: these were on the ancestor genomes-we did not run this on the descendant strains

#ancestral assemblies can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA310138/

#this was looped to run for all 28 of our ancestor strains


#load libraries
module load StdEnv/2020
module load gcc/9.3.0
module load prokka/1.14.5

prokka /path/to/ancestor/rhizobiu_ancestor_assembly.fasta --outdir /prokka/output/predicted/proteincoding/genes
#while there is a lot of output, the salient pieces are: prokka_output.gff & prokka_output.tsv & prokka_output.faa
#these are used in downstream analyses
