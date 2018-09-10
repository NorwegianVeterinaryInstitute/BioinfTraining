#!/bin/bash
## data submission script
#SBATCH --account=nn9305k
#SBATCH --time=1:00:00

##memory specs
#SBATCH --array=0-29  # array with 100 jobs
#SBATCH --mem-per-cpu=3800M
#SBATCH --cpus-per-task=4
#SBATCH --job-name=fastqc
#SBATCH --mail-user=thomas.haverkamp@vetinst.no
#SBATCH --mail-type=ALL
#SBATCH --output=/work/users/thhaverk/slurm_out/fastqc.%j.out

## Set up job environment
source /cluster/bin/jobsetup

## a slurm script to run fastqc on unprocessed illumina sequences

## created by Thomas H.A. Haverkamp
## version 1.0 (20180907)

echo jobstarting...

echo checking location
pwd -L

#changing to directory with input data
cd ../fastq_files

# loading modules
module load fastqc/0.11.2

#Generate array of unique datasets
FILES=($(ls -1))

# counting number of datasets
echo "Number of samples:" ${#FILES[*]}

#Generate name of current sample
#[$SLURM_ARRAY_TTASK_ID] is the element index. Increases with one per job, i.e. one new sample for every job)
CURRENT_FILE=${FILES[$SLURM_ARRAY_TASK_ID]} 

echo "Current input-file name is" $CURRENT_FILE

#running fastqc and save in folder: fastqc_results
fastqc $CURRENT_FILE -o ../fastqc_results --threads 4

echo finished with fastqc