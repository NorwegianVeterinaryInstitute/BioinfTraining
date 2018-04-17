#!/bin/bash
## This script generate the Kraken database.
## In addition to installing the standard Kraken libraries (bacteria, archaea, plasmids and viruses), 
## it will download PROTOZOAN and FUNGAL sequences from RefSeq, and add them to the libraries.
## Finally, it will build the actual database.
## If needed, hugemem nodes can be used.

##Specify job settings
#SBATCH --account=nn9305k
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=30000M
#SBATCH --cpus-per-task=8
#SBATCH --partition=hugemem
#SBATCH --job-name=build_db
#SBATCH --mail-user=thhaverk@ibv.uio.no   ##CHANGE
#SBATCH --mail-type=ALL
#SBATCH --output=/work/users/thhaverk/slurm_out/krak_build-%j-%N.out

## Set up job environment
source /cluster/bin/jobsetup

# creating variables for folders
KRAKEN_DB=(/work/users/thhaverk/kraken_test/kraken_db)

## Go to database folder
cd $KRAKEN_DB

## Load necessary modules
module load kraken/1.0

## Download taxonomy
kraken-build --download-taxonomy --db kraken_db

## Download and install standard genomic libraries
echo Downloading standard Kraken libraries

kraken-build --download-library plasmid --db kraken_db
kraken-build --download-library viral --db kraken_db
kraken-build --download-library archaea --db kraken_db
kraken-build --download-library bacteria --db kraken_db

echo Done with downloading and installing standard Kraken libraries

## Download protozoan and fungal genomes (no need to call perl, done in module)
#echo "Downloading custom libraries (fungi and protozoa)"

#$DOWNLOAD_SCRIPTS/download_fungi.pl
#$DOWNLOAD_SCRIPTS/download_protozoa.pl

#echo Done with downloading fungal and protozoan genomes

## Download the taxonomy.
## kraken-build --download-taxonomy --db kraken_db

## Add protozoa and fungi to library
#for dir in fungi protozoa; do
#        for fna in `ls $dir/*.fna`; do
#                kraken-build --add-to-library $fna --db kraken_db
#        done
#done

## Build the actual database (from libraries)
kraken-build --build --db kraken_db --threads 8

## Remove all intermediate files
kraken-build --clean --db kraken_db --threads 8

echo Checking files

ls

echo Finished

