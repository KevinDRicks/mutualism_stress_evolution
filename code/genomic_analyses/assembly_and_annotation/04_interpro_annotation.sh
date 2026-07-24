#!/bin/bash
#codes passes the predicted protein sequences from the ancestor rhizobia genomes through interproscan for annotation
#these protein sequences were extracted from the ancestor genome in the following script: 
#02_id_protein_coding.sh & 03_split_protein_coding_by_replicon.R
#note: this script was initially run on Compute Canada server

#!/bin/bash
#code passes the protein sequences from the ancestral rhizobia genomes through interproscan for annotation
#these protein sequences were extracted from the ancestor genome in the following script: 
#id_protein_coding.sh & split_protein_coding_by_replicon.R
#note: this script was initially run on Compute Canada server
#note: these were on the ancestor genomes-we did not run this on the descendant strains


#this was looped to run for all 28 of our ancestor strains
#we can annotate these ancestors to understand where introduced mutations are occuring on our samples


#load libraries
module load interproscan/5.73-104.0
module load interproscan_data/5.73-104.0
module load StdEnv/2023
module load java/21.0.1


#note, the path will need to be altered path in order to call the interproscan.sh script
#th enumber of plasmids will also depend on plasmid number

interproscan.sh -i chromosome_AA.faa -b chromosome_interpro -goterms
interproscan.sh -i plasmidA_AA.faa -b plasmidA_interpro -goterms
interproscan.sh -i plasmidB_AA.faa -b plasmidB_interpro -goterms
interproscan.sh -i plasmidC_AA.faa -b plasmidC_interpro -goterms
interproscan.sh -i plasmidD_AA.faa -b plasmidD_interpro -goterms