#!/bin/bash
#SBATCH --time=02:00:00
####SBATCH --output=output_spine.txt  
#SBATCH --error=error_spine.txt
#SBATCH --mem=96G

#this script is one of multiple methods used to detect introduced mutations in our evolved strains
#in this script, assemblies of our evolved strains are aligned and compared against previously 
#constructed high quality ancestor assemblies 
#this follows a custom pipeline previous developed in our working group, described here: 
#https://github.com/Alan-Collins/Spine-Nucmer-SNPs
#custom perl script, spine.pl, was downloaded and added as a path

#assemblies for our evolved strains can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1294393/
#ancestral assemblies can be found at https://www.ncbi.nlm.nih.gov/bioproject/PRJNA310138/

#this script must by run for each evolved strain (modified to run in a loop)

module load StdEnv/2023
module load perl/5.36.1
module load mummer/4.0.0rc1
  

mkdir ASSEMBLIES
mkdir NUCMER
mkdir SPINE

target_evolved_assembly="/path/to/evolved/strain/assembly/[STRAIN_ID]_evolved_assembly.fasta"
target_ancestor_assembly="/path/to/ancestral/assembly/[ANCESTOR_STRAIN_ID]_ancestor_assembly.fasta"


cp $target_ancestor ASSEMBLIES/target_evolved_assembly
cp $target_evolve ASSEMBLIES/target_evolved_assembly
  
cd ASSEMBLIES
ls | awk 'BEGIN { FS="\t"; OFS="\t" } { print "../ASSEMBLIES/"$1, $1, "fasta" }' > ../SPINE/config.txt
  
cd ../SPINE
perl ~/modules/added_path/spine.pl  -f config.txt -t 16

cd ../NUCMER
ls ../SPINE/*.core.fasta | while read i; do acc=${i%.core*}; acc=${acc#../SPINE/output.}; nucmer --prefix=${acc}_core ../SPINE/output.backbone.fasta $i; delta-filter -r -q ${acc}_core.delta > ${acc}_core.filter; show-snps -Clr ${acc}_core.filter > ${acc}_core.snps; done

