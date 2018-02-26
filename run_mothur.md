# a slurm script to run mothur

The script below can be used on the computing cluster abel, or any other cluster using the slurm queing system, to run mothur on a dataset with a textfile that contains the commands to analysis the dataset of interest.

```
#!/bin/bash

## a script to run mothur on abel with an input file having commands

## Version 20180226

#SBATCH --account=nn9305k
#SBATCH --time=0:30:00
#SBATCH --mem-per-cpu=3800M
#SBATCH --cpus-per-task=16
#SBATCH --job-name=mothur

## adding email notification
## replace mail address with your own.
#SBATCH --mail-user=YOURMAIL@DOMAIN
#SBATCH --mail-type=ALL

## adding logfile creation
## replace username with your own username
#SBATCH --output=/work/projects/nn9305k/home/thhaverk/amplicons/slurm_out/mothur_run_%j.txt


## check location
pwd

## changing to MiSeq_SOP directory

cd ../MiSeq_SOP/

## check input file
ls

##loading modules
module load mothur/1.38.1

mothur commands_edit.txt

echo finished 
```


Make sure you have the following datastructure in your home area on abel

`amplicons`

In that folder create three directories:

```
MiSeq_SOP
slurm_out
slurm_scripts
```

save the above script in the folder: `slurm_scripts`

run the script from the folder with the command:

```
sbatch run_mothur.slurm
```
