#!/bin/bash

## created by T.H.A. Haverkamp
## creation date: 20180416

## a slurm script to classify metagenomics sequences using Kraken.
## We use the minikraken database with 4 cpus for classification

## requesting job
#SBATCH --account=nn9305k
#SBATCH --time=1:00:00

##memory specs
#SBATCH --mem-per-cpu=3800M
#SBATCH --cpus-per-task=4
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

## Creating variables to use

MYFOLDER=(/projects/nn9305k/home/thhaverk/kraken_test)  ## identify which dataset is analyzed
SAMPLE=sample_K1  ## sample identification
KRAKEN_DB=(/work/projects/nn9305k/tmp/minikraken_20171013_4GB/)  ## the (mini-) kraken database location
OUT_RESULT=(/projects/nn9305k/home/thhaverk/kraken_test_results)  ## output folder

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
--threads 4 \
--fastq-input \
--gzip-compressed \
--classified-out $SAMPLE.classified.fastq \
--unclassified-out $SAMPLE.unclassified.fastq \
--output $SAMPLE.Kr_out.txt \
--preload \
--paired \
${RAWFILE[0]} \
${RAWFILE[1]}

echo finished with classification, generating reports

## commands to generate kraken & krona reports at the two confidence levels
echo the taxonomic labels 
kraken-translate --db $KRAKEN_DB $MYFOLDER.Kr_out.txt > $SAMPLE.Kr_out.labels.txt

## running reports with different confidence levels.
for number in 0 0.20 ; do
	echo filtering kraken output;
	kraken-filter --db $KRAKEN_DB \
	--threshold $number \
	$SAMPLE.Kr_out.txt > $SAMPLE.Kr_out.C.$number.txt;
	
	echo generating classification report;
	kraken-report --db $KRAKEN_DB \
	$SAMPLE.Kr_out.C.$number.txt > $SAMPLE.Kr_classification.C.$number.txt;
	
	echo generating mpa classification report;
	kraken-mpa-report --db $KRAKEN_DB \
	$SAMPLE.Kr_out.C.$number.txt > $SAMPLE.Kr_classification.mpa.C.$number.txt;ls
	
	echo generating krona plot;
	cat $SAMPLE.Kr_out.C.$number.txt | cut -f 2,3 > $SAMPLE.Kr_out.C.$number.krona
	ktImportTaxonomy $SAMPLE.Kr_out.C.$number.krona -o krona.$SAMPLE.html
done

## Cleaning up the folders.
echo gzipping the fastq output files
gzip -v $SAMPLE.*.fastq
echo 

echo moving output data to kraken_results
mkdir $OUT_RESULT
rsync -r --remove-source-files *.txt $OUT_RESULT
rsync -r --remove-source-files $SAMPLE.*.fastq.gz $OUT_RESULT
rsync -r --remove-source-files *.htm* $OUT_RESULT
rsync -r --remove-source-files *.krona $OUT_RESULT

## checking if the destination directory contains all the data
echo checking contents of $OUT_RESULT
cd $OUT_RESULT
pwd
ls

echo finished.




