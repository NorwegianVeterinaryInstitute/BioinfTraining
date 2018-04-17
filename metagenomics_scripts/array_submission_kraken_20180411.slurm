#!/bin/bash
## data submission script
## created by T.H.A. Haverkamp
## creation date: 20180411

## a slurm script to classify metagenomics sequences using Kraken.
## This script will use an array to classify shotgun reads from different samples
## each dataset is run in one slurm job and submitted with "arrayrun"
## this submission script needs a long time and little memory.
## final step is generation of kronaplots from all samples

##specify job settings
#SBATCH --account=cees
#SBATCH --time=140:00:00
#SBATCH --mem-per-cpu=200M
#SBATCH --job-name=submit
#SBATCH --mail-user=thhaverk@ibv.uio.no
#SBATCH --mail-type=ALL

## creating output file for debugging and error analysis.
#SBATCH --output=/work/users/thhaverk/slurm_out/kraken_array_submit-%j-%N.out

## Set up job environment
source /cluster/bin/jobsetup

## This script will use an arrayrun to process multiple files for the same analysis
## each file is run in one slurm job

#Where is the script started
echo checking location
pwd

# changing to the data directory
cd /work/users/thhaverk/kraken_test/fastq

pwd
echo

# collecting all the dataset files into an array called FOLDERS
FOLDERS=($(ls -1))

# Than we count the number of folders (samples) to determine the number of tasks
NUM_TASKS=${#FOLDERS[*]}

# determine the maximum number of tasks
MAX_ID=$((NUM_TASKS - 1))

# running the arrayrun command with tasks id's ranging from 0-MAX_ID
arrayrun 0-$MAX_ID /work/users/thhaverk/slurm_scripts/workscript_kraken_20180411.slurm

# move the slurm.output files in the data to output_slurms
mv slurm-*.out /work/users/thhaverk/slurm_out

## change to results directory
cd /work/users/thhaverk/kraken_test/results/

#loading module kronatools 
module load kronatools/2.7

# generating krona html document containing all samples
ktImportTaxonomy */*.krona -o Kraken_classifications.html

echo finished



