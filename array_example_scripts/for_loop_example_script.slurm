#!/bin/bash
#SBATCH --account=nn9305k
#SBATCH --time=5:00:00

##memory specs
#SBATCH --mem-per-cpu=3800M
#SBATCH --cpus-per-task=4
#SBATCH --job-name=fastqc
#SBATCH --output=/work/users/thhaverk/slurm_out/fastqc.%j.out

# run script by calling a for loop, where each folder to process is called.
# example: 
# for file in *; do sbatch fastq.slurm $file; done

# listing input data
echo processing dataset ${1}   

# load module
module load fastqc/0.11.2

#running fastqc
fastqc ${1} -o ../../fastqc_results -t 32;

echo finished with folder ${1}
