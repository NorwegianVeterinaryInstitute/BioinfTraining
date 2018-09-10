#!/bin/bash
#SBATCH --account=nn9305k
#SBATCH --time=5:00:00

## Memory specs
#SBATCH --mem-per-cpu=3800M
#SBATCH --cpus-per-task=4
#SBATCH --job-name=fastqc
#SBATCH --output=/work/users/thhaverk/slurm_out/fastqc.%j.out

## Set up job environment
source /cluster/bin/jobsetup

## Make sure you know where u are when u submit the script
JOBHOME=$SUBMITDIR
JOBNAME="fastqc"

## A slurm script to run fastqc on illumina data, written by Thomas Haverkamp,
## version 1.0.

echo checking location
pwd -L

#changing to directory with input data
cd ../fastq_files

# The variable FILES contains all the files for analysis
FILES=($(ls -1))

#create a variable using the filename and combine it with the TASK_ID from arrayrun
CURRENT_FILE=${FILES[$TASK_ID]}

echo current dataset is $CURRENT_FILE

# loading modules
module load fastqc/0.11.2 
module list

#running fastqc and save in folder: fastqc_results
fastqc $CURRENT_FILE -o ../fastqc_results -t 4

echo finished with fastqc