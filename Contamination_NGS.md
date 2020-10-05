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
Their location on Saga is here: `/cluster/projects/nn9305k/db_flatfiles/contamination_hosts`.

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

```
MY_DATABASE=/cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312
echo $MY_DATABASE
```

Now we can use the variable to indicate where to find the database for our kraken analysis. To use a variable we create we have to add a `$` to the variable name. In our case that would be: $MY_DATABASE. So the command we would write without our variable looks like this:

```
kraken2 --db /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312 \
  --threads 4 \
  --output 18-Cjejuni-927.kraken2.out \
  --report 18-Cjejuni-927.kraken2_report.txt \
  --minimum-base-quality 20 \
  --paired \
  --gzip-compressed \
  18-Cjejuni-927_Subsampled_L008_R*_0089.fastq.gz 
```
Or with the variable we just created:

```
kraken2 --db $MY_DATABASE \
  --threads 4 \
  --output 18-Cjejuni-927.kraken2.out \
  --report 18-Cjejuni-927.kraken2_report.txt \
  --minimum-base-quality 20 \
  --paired \
  --gzip-compressed \
  18-Cjejuni-927_Subsampled_L008_R*_0089.fastq.gz 
```

Now, Set-up the same analysis for the other dataset: `17-Cjejuni-CCUG11284T_Subsampled_L008_R1_0087.fastq.gz`.

#### Question: What is the effect of the option `--minimum-base-quality 20`?

The kraken run, created two files: 

* `18-Cjejuni-927.kraken2.out` - the classification for each single read.
* `18-Cjejuni-927.kraken2_report.out`. - a summary of the classifications per taxon.

Open the later one with the `less` command to explore the classifications of the reads belonging to the *Campylobacter* isolate

```
less 18-Cjejuni-927.kraken2_report.out
```

#### Question: Can you spot species that should not be in there? Why do think these are contaminating sequences?.

#### Question: What is the case for the unclassified reads? Should those be removed?

The Illumina sequences of the *Campylobacter jejuni* strain 18-Cjejuni-927, have contamination reads. The most obvious ones (Phix and *Homo sapiens*) can be removed using the reference genomes for both species.

**Note!**: Under normal circumstances we would always do a quality trimming and removal of the obvious contaminations, before we use Kraken2 or something similar to see if there is contaminating data in our datasets. In addition, most people find out that their data is contaminated with something, when something does not work in the way they expected it.

## Removing PhiX

Phix is a small phage originally isolated from *`E.coli`* (https://en.wikipedia.org/wiki/Phi_X_174). Nowadays it is know as the piece of DNA that is added to Illumina sequencing libraries. Most of the PhiX reads that are sequenced together with your samples are removed from the data. Nonetheless, there will always be reads that are from Phix that manage to end up in your dataset. One way of removing them from your data, is to do the genome assembly, and then look among the small contigs for a contig of about 5386 bases and remove that. There are examples of submitted microbial genomes where phiX was included in the submission as well altough it was never part of the original genome.
Here we try to prevent that by removing reads that match the PhiX genome before assembly. This is how.

We first make a variable that contains the full path the the PhiX genome, just to make the command for cleaning easier

```
PHIX=/cluster/projects/nn9305k/db_flatfiles/contamination_hosts/PhiX.fasta
echo $PHIX
````
We will use a tool called BBduk with is part of the BBTools packages (https://jgi.doe.gov/data-and-tools/bbtools/)

We activate the tool like this:
```
conda activate BBTools
bbduk.sh -h 
```
Now we use our input datasets to clean our samples.

```
bbduk.sh threads=4 \
 ref=$PHIX \
 in1=18-Cjejuni-927_Subsampled_L008_R1_0089.fastq.gz in2=18-Cjejuni-927_Subsampled_L008_R2_0089.fastq.gz \
 out=18.phix_R1.fastq out2=18.phix_R2.fastq \
 stats=18.Phix_stats.txt
```
The following is optional, it does add quality trimming and control to the analysis. To add this to the previous command you can add the following:
```
k=31 ktrim=r mink=11 hdist=1 tbo=f tpe=f qtrim=r trimq=15 maq=15 minlen=36 forcetrimright=149 
```

#### Question: How many reads were identify as PhiX in both genomes for this analysis?

#### Question: If you use quality trimming and control, what does that do with identifying PhiX? 

## Removing Human DNA 

To remove human contamination, we use the human reference genome `Hg19`. The maker from BBTools, has taken this genome and then masked all the regions in the genome that match ribosomal sequences, animals genomes, all plant genomes, fungal genomes. That gives a human genome with regions that are only found in the human genome. So the number of false positive matches is reduced.
The removal of reads matching the Human genome is done with bbduk.

With the current interactive job we have a problem with the memory when we do this with the human genome. We need more memory. Thus we need to cancel our current interactive job and recreate a new one.

```
srun --account=nn9305k --qos=devel --mem-per-cpu=4800M --cpus-per-task=10 --time=0:30:00 --pty bash -i
```


We first make a variable that contains the full path the the PhiX genome, just to make the command for cleaning easier

```
HUMAN=/cluster/projects/nn9305k/db_flatfiles/contamination_hosts/hg19_main_mask_ribo_animal_allplant_allfungus.fasta
echo $HUMAN
```
The input for this is the output from the PhiX removal step.
```

bbduk.sh threads=10 \
 ref=$HUMAN \
 in1=18.phix_R1.fastq in2=18.phix_R2.fastq \
 out=18.phix_human_R1.fastq out2=18.phix_human_R2.fastq \
 stats=18.human_stats.txt \
 k=31 ktrim=r mink=11 hdist=1 tbo=f tpe=f qtrim=r trimq=15 maq=15 minlen=36 forcetrimright=149 
```
When your job went through than you can check the results file.

## Detecting other possible contaminating sequences.

screen

srun --account=nn9305k --qos=devel --mem-per-cpu=4800M --cpus-per-task=4 --time=0:30:00 --pty bash -i

### Input files
IF1="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step2_Dehumaned_Reads/Step2_18-Cjejuni-927_Subsampled_L008_R1_0089.fastq"
IF2="/cluster/projects/nn9305k/tutorial/20201005_contamination/Step2_Dehumaned_Reads/Step2_18-Cjejuni-927_Subsampled_L008_R1_0089.fastq"

conda activate kraken2

kraken2 --db /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312 --threads 4 --output 18-Cjejuni-927.kraken2.out --report 18-Cjejuni-927.kraken2_report.txt --minimum-base-quality 20 --paired $IF1 $IF2


kraken2 --db /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312 --threads 4 --output 17-Cjejuni-CCUG11284T.kraken2.out --report 17-Cjejuni-CCUG11284T.kraken2_report.txt --minimum-base-quality 20 --paired --gzip-compressed 17-Cjejuni-CCUG11284T_Subsampled_L008_R1_0087.fastq.gz 17-Cjejuni-CCUG11284T_Subsampled_L008_R2_0087.fastq.gz

conda deactivate
