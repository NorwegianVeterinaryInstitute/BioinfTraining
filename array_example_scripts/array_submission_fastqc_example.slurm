#!/bin/bash
## data submission script

##specify job settings
#SBATCH --account=nn9305k
#SBATCH --time=5:00:00

##memory specs
#SBATCH --mem-per-cpu=3800M
#SBATCH --cpus-per-task=4
#SBATCH --job-name=fastqc
#SBATCH --job-name=submit
#SBATCH --mail-user=thhaverk@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=/work/users/thhaverk/slurm_out/submit.%j.out

## Set up job environment
source /cluster/bin/jobsetup

## This script will use an arrayrun to process multiple files for the same analysis
## each file is run in one slurm job

echo jobstarting...

echo checking location
pwd -L

#changing to directory with input data
cd ../fastq_files

# collecting all the PAIRED dataset files into an array called FOLDERS
FOLDERS=($(ls -1))

# Than we count the number of files to determine the number of tasks
NUM_TASKS=${#FOLDERS[*]}

# determine the maximum number of tasks
MAX_ID=$((NUM_TASKS - 1))

# running the arrayrun command with tasks id's ranging from 0-MAX_ID
arrayrun 0-$MAX_ID ../slurm_scripts/workscript_fastqc_example.slurm

echo finished

