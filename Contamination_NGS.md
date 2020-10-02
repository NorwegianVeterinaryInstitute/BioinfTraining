# Removing contamination

## Tutorial datasets
The data for this tutorial can be found on Saga in the folder: 

```
/cluster/projects/nn9305k/tutorial/20201005_contamination/
```
The data consists of two raw Illumina datasets from the Airsample project that contain the reads for two *Campylobacter* isolates:

```
17-Cjejuni-CCUG11284T: 17-Cjejuni-CCUG11284T_Subsampled_L008_R1/2_0087.fastq.gz
18-Cjejuni-927: 18-Cjejuni-927_Subsampled_L008_R1/2_0089.fastq.gz
```

Both datasets were sequenced vry deeply, and here we will use a subsampled dataset for a tutorial on how to clean such datasets.

Besides these two datasets we also need datasets that are from common contaminants.
* the phage PhiX - Phix.fasta
* The human genome - hg19_main_mask_ribo_animal_allplant_allfungus.fasta 

And we need a database that we can use to classify our reads or our contigs.
We will use a minikraken database retrieved from here: (http://ccb.jhu.edu/software/kraken2/downloads.shtml)
We use the database with the human genome in there. You can find it on Saga at this location:
` /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312`


#### Question: Why can we expect both PhiX and Human DNA to be present in the sequence data of these two *Campylobacter* isolates?

Both genomes are found in the location: `/cluster/projects/nn9305k/db_flatfiles/contamination_hosts`.

## The Tutorial
To start the tutorial move to your work directory, and create a directory called : `contamination`
The commands:

```
cd $USERWORK
mkdir contamination
cd contamination
```
Now we retrieve the fastq files for both datasets so each of us can do the tutorial

```
 rsync -rauWP /cluster/projects/nn9305k/tutorial/20201005_contamination/*.gz ./
 ls
```
Now we have the two files ready for this excercise.

Normally we would start with a proper clean-up of the datasets, using trimmomatic or similar tool, to remove bases and reads of bad quality, too short reads ect. But for the purpose of this excercise we skip that for now.

# Kraken 2 classification.
We will start with the classification of the raw reads that came straight from the sequencing machine. We do this to get an idea of the problems. We will use kraken2 for this and classify the reads against the minikraken database

The commands

```
screen
```
set up a interactive job, using the development mode (--qos=devel)
```
srun --account=nn9305k --qos=devel --mem-per-cpu=4800M --cpus-per-task=4 --time=0:30:00 --pty bash -i
```
Now activate the kraken environment.
```
conda activate kraken2
```
Since the location of the minikraken database is a long line, we want to shorten that a bit. We can do that by putting the directory PATH into a a variable. The variable is kept only for this session.

KR2-DB=(/cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312



kraken2 --db /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312 \
  --threads 4 \
  --output 18-Cjejuni-927.kraken2.out \
  --report 18-Cjejuni-927.kraken2_report.txt \
  --minimum-base-quality 20 \
  --paired \
  --gzip-compressed \
  18-Cjejuni-927_Subsampled_L008_R*_0089.fastq.gz 
```
and
```
kraken2 --db /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312 --threads 4 --output 17-Cjejuni-CCUG11284T.kraken2.out --report 17-Cjejuni-CCUG11284T.kraken2_report.txt --minimum-base-quality 20 --paired --gzip-compressed 17-Cjejuni-CCUG11284T_Subsampled_L008_R1_0087.fastq.gz 17-Cjejuni-CCUG11284T_Subsampled_L008_R2_0087.fastq.gz
```



## Removing PhiX


### Input files
IF1="/cluster/projects/nn9305k/tutorial/20201005_contamination/18-Cjejuni-927_Subsampled_L008_R1_0089.fastq.gz"
IF2="/cluster/projects/nn9305k/tutorial/20201005_contamination/18-Cjejuni-927_Subsampled_L008_R2_0089.fastq.gz"

### Output files
OF1="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step1_PhiX_Cleaned_Reads/Step1_18-Cjejuni-927_Subsampled_L008_R1_0089.fastq"
OF2="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step1_PhiX_Cleaned_Reads/Step1_18-Cjejuni-927_Subsampled_L008_R2_0089.fastq"

conda activate BBTools

bbduk.sh threads=5 ref=phix_location,adapter_location in1=INF1 in2=INF2 out=OF1 out2=OF2 k=31 ktrim=r mink=11 hdist=1 tbo=f tpe=f qtrim=r trimq=15 maq=15 minlen=36 forcetrimright=149 stats=stats.txt > log_file

## Removing Human DNA 

### Input files
IF1="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step1_PhiX_Cleaned_Reads/Step1_18-Cjejuni-927_Subsampled_L008_R1_0089.fastq"
IF2="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step1_PhiX_Cleaned_Reads/Step1_18-Cjejuni-927_Subsampled_L008_R2_0089.fastq"

### Output files
OF1="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step2_Dehumaned_Reads/Step2_18-Cjejuni-927_Subsampled_L008_R1_0089.fastq"
OF2="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step2_Dehumaned_Reads/Step2_18-Cjejuni-927_Subsampled_L008_R1_0089.fastq"

bbduk.sh threads=5 ref=$human_genome in1=$IF1 in2=$IF2 out=$OF1 out2=$OF2 k=31 ktrim=r mink=11 hdist=1 tbo=f tpe=f qtrim=r trimq=15 maq=15 minlen=36 forcetrimright=149 stats=stats.txt > log_file

conda deactivate

## Other Contamination

screen

srun --account=nn9305k --qos=devel --mem-per-cpu=4800M --cpus-per-task=4 --time=0:30:00 --pty bash -i

### Input files
IF1="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step2_Dehumaned_Reads/Step2_18-Cjejuni-927_Subsampled_L008_R1_0089.fastq"
IF2="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step2_Dehumaned_Reads/Step2_18-Cjejuni-927_Subsampled_L008_R1_0089.fastq"

conda activate kraken2

kraken2 --db /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312 --threads 4 --output 18-Cjejuni-927.kraken2.out --report 18-Cjejuni-927.kraken2_report.txt --minimum-base-quality 20 --paired $IF1 $IF2


kraken2 --db /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312 --threads 4 --output 17-Cjejuni-CCUG11284T.kraken2.out --report 17-Cjejuni-CCUG11284T.kraken2_report.txt --minimum-base-quality 20 --paired --gzip-compressed 17-Cjejuni-CCUG11284T_Subsampled_L008_R1_0087.fastq.gz 17-Cjejuni-CCUG11284T_Subsampled_L008_R2_0087.fastq.gz

conda deactivate
