#!/bin/bash

## created by T.H.A. Haverkamp
## creation date: 20180410

## a slurm script to classify metagenomics sequences using Kraken.
## This script will use an array to classify shotgun reads from different samples
## each dataset is run in one slurm job
## the job requires hugemem due to the size of the classification database.

##specify job settings
#SBATCH --array=0-1  ## this asks for 2 jobs to process 2 samples.
#SBATCH --account=nn9305k
#SBATCH --time=1:00:00

##memory specs
#SBATCH --mem-per-cpu=32000M
#SBATCH --cpus-per-task=16
#SBATCH --partition=hugemem
#SBATCH --job-name=kraken
#SBATCH --mail-user=thhaverk@ibv.uio.no
#SBATCH --mail-type=ALL

## creating output file for debugging and error analysis.
#SBATCH --output=/work/users/thhaverk/slurm_out/krak_classification-%j-%N.out

## Set up job environment
source /cluster/bin/jobsetup

## make sure you know where u are when u submit the script
JOBHOME=$SUBMITDIR

echo jobstarting...

## change to work directory
cd /work/users/thhaverk/kraken_test/fastq

echo checking location
pwd -L

## Creating variables for the input data, database and output
FOLDERS=($(ls -1))
MYFOLDER=${FOLDERS[$SLURM_ARRAY_TASK_ID]}
KRAKEN_DB=(/work/users/thhaverk/kraken_test/kraken_db_20180410/)
OUT_RESULT=(/work/users/thhaverk/kraken_test/results/$MYFOLDER)  ## output folder

## starting analysis
echo current datasets is: $MYFOLDER

## loading modules
echo Loading required modules
module purge
module load kraken/1.0
module load KronaTools/2.7

echo loaded modules
module list

## Running kraken on the paired fastq.gz files in the folder
cd $MYFOLDER
pwd -L

echo making an array of the fastq.gz files
RAWFILE=($(ls -1 ${prefix}*.fastq.gz))

echo running kraken on $SCRATCH
echo fastq input pair1: ${RAWFILE[0]}, pair2: ${RAWFILE[1]}
 
kraken -db $KRAKEN_DB \
--threads 16 \
--fastq-input \
--gzip-compressed \
--classified-out $MYFOLDER.classified.fastq \
--unclassified-out $MYFOLDER.unclassified.fastq \
--output $MYFOLDER.Kr_out.txt \
--preload \
--paired \
${RAWFILE[0]} \
${RAWFILE[1]}

echo finished with classification, generating reports

## Parsing the kraken output to generate reports of the classification
echo the taxonomic labels 
kraken-translate --db $KRAKEN_DB $MYFOLDER.Kr_out.txt > $MYFOLDER.Kr_out.labels.txt

echo  generating confidence reports - level 0.00
kraken-filter --db $KRAKEN_DB \
 --threshold 0 $MYFOLDER.Kr_out.txt > $MYFOLDER.Kr_out.conf_000.txt

echo classification report
kraken-report --db $KRAKEN_DB $MYFOLDER.Kr_out.conf_000.txt > $MYFOLDER.Kr_classification_c000.txt

echo mpa classification report
kraken-mpa-report --db $KRAKEN_DB $MYFOLDER.Kr_out.conf_000.txt > $MYFOLDER.Kr_classification.mpa.c000.txt

## Generating a Krona plot to vizualize the results

## generating a report with a higher precision : 99.13
echo  generating confidence reports - level 0.20
kraken-filter --db $KRAKEN_DB \
 --threshold 0.20 $MYFOLDER.Kr_out.txt > $MYFOLDER.Kr_out.conf_0.20.txt

echo classification report
kraken-report --db $KRAKEN_DB $MYFOLDER.Kr_out.conf_0.20.txt > $MYFOLDER.Kr_classification_c0.20.txt

echo mpa classification report
kraken-mpa-report --db $KRAKEN_DB $MYFOLDER.Kr_out.conf_0.20.txt > $MYFOLDER.Kr_classification.mpa.c0.20.txt

##Generating a Krona plot to vizualize the results
cat $MYFOLDER.Kr_out.conf_000.txt | cut -f 2,3 > $MYFOLDER.Kr_out.conf_000.krona
cat $MYFOLDER.Kr_out.conf_0.20.txt | cut -f 2,3 > $MYFOLDER.Kr_out.conf_0.20.krona

ktImportTaxonomy $MYFOLDER.Kr_out.conf_000.krona $MYFOLDER.Kr_out.conf_0.20.krona -o krona.$MYFOLDER.html


## cleaning up the data
echo gzipping the fastq output files
gzip -v *.fastq
echo 

echo moving output data to kraken_results
mkdir $OUT_RESULT
rsync -r --remove-source-files *.txt $OUT_RESULT
rsync -r --remove-source-files *.fastq.gz $OUT_RESULT
rsync -r --remove-source-files *.htm* $OUT_RESULT
rsync -r --remove-source-files *.krona $OUT_RESULT

## checking if the destination directory contains all the data
echo checking contents of $OUT_RESULT
cd $OUT_RESULT
pwd
ls

echo finished.




