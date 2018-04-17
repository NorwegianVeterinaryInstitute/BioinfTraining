#!/bin/bash

## created by T.H.A. Haverkamp
## creation date: 20180411

## a slurm script to classify metagenomics sequences using Kraken.
## the job requires hugemem due to the size of the classification database.

#SBATCH --account=nn9305k
#SBATCH --time=1:00:00

##memory specs
#SBATCH --mem-per-cpu=32000M
#SBATCH --cpus-per-task=16
#SBATCH --partition=hugemem
#SBATCH --job-name=kraken

## creating output file for debugging and error analysis.
#SBATCH --output=/work/users/thhaverk/slurm_out/kraken_classification-%j-%N.out

## Set up job environment
source /cluster/bin/jobsetup

##make sure you know where u are when u submit the script
JOBHOME=$SUBMITDIR
JOBNAME="kraken"

echo jobstarting...

echo checking location
pwd -L

# Creating variables to use
FOLDERS=($(ls -1))   ##  FILES contains all the files for analysis
MYFOLDER=${FOLDERS[$TASK_ID]}  ## identify which dataset is analyzed
KRAKEN_DB=(/work/users/thhaverk/kraken_test/kraken_db_20180410/)  ## the kraken database location
OUT_RESULT=(/work/users/thhaverk/kraken_test/results/$MYFOLDER)  ## output folder

## starting analysis
echo current dataset is: $MYFOLDER

echo Loading required modules
module purge
module load kraken/1.0
module load kronatools/2.7

echo loaded modules
module list

## running kraken on the paired fastq.gz files in the folder
cd $MYFOLDER
pwd -L

echo making an array of the fastq.gz files
RAWFILE=($(ls -1 ${prefix}*.fastq.gz))

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

## commands to generate kraken & krona reports at the two confidence levels
echo the taxonomic labels 
kraken-translate --db $KRAKEN_DB $MYFOLDER.Kr_out.txt > $MYFOLDER.Kr_out.labels.txt

#running reports with different confidence levels.
for number in 0 0.20 ; do
	echo filtering kraken output;
	kraken-filter --db $KRAKEN_DB \
	--threshold $number \
	$MYFOLDER.Kr_out.txt > $MYFOLDER.Kr_out.C.$number.txt;
	
	echo generating classification report;
	kraken-report --db $KRAKEN_DB \
	$MYFOLDER.Kr_out.C.$number.txt > $MYFOLDER.Kr_classification.C.$number.txt;
	
	echo generating mpa classification report;
	kraken-mpa-report --db $KRAKEN_DB \
	$MYFOLDER.Kr_out.C.$number.txt > $MYFOLDER.Kr_classification.mpa.C.$number.txt;ls
	
	echo generating krona plot;
	cat $MYFOLDER.Kr_out.C.$number.txt | cut -f 2,3 > $MYFOLDER.Kr_out.C.$number.krona
	ktImportTaxonomy $MYFOLDER.Kr_out.C.$number.krona -o krona.$MYFOLDER.html
done

## Cleaning up the folders.
echo gzipping the fastq output files
gzip -v $MYFOLDER.*.fastq
echo 

echo moving output data to kraken_results
mkdir $OUT_RESULT
rsync -r --remove-source-files *.txt $OUT_RESULT
rsync -r --remove-source-files $MYFOLDER.*.fastq.gz $OUT_RESULT
rsync -r --remove-source-files *.htm* $OUT_RESULT
rsync -r --remove-source-files *.krona $OUT_RESULT

#checking if the destination directory contains all the data
echo checking contents of $OUT_RESULT
cd $OUT_RESULT
pwd
ls

echo finished.




