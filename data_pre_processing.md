# Data pre-processing

## Datasets used for this session

Dataset used during this session can be found in the following location within abel:

`/work/projects/nn9305k/tmp/Files_for_Dec14/`

## NB: Replace <your_user_name> with your abel username

Create a new folder called _Data_pre_processing_Dec14_ in your home area _/work/projects/nn9305k/home/<your_user_name>/_ and move there.

`cd /work/projects/nn9305k/home/<your_user_name>/Data_pre_processing_Dec14`

Create three folders here.

`mkdir data`
`mkdir raw_fastqc`
`mkdir trim`

Move to the _data_ folder.

`cd /work/projects/nn9305k/home/<your_user_name>/Data_pre_processing_Dec14/data`

Type the following command to link the files (not copy):

`ln -s /work/projects/nn9305k/tmp/Files_for_Dec14/*fq.gz .`


## Fastq quality check

We will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to check the quality of raw sequenced data 

--------

**Task**
1. Run _fastqc_ on the file _Ha_R1.fq.gz_ (in the login node).
2. type in the following (don't include the `$`).

```
$ fastqc Ha_R1.fq.gz
```

3. Try to find help for `fastqc` and discuss what flags one can use to process multiple samples
4. Use `SLURM` to process the other four files


```
#!/bin/bash
#
# Job name:
  
#SBATCH --job-name=raw_fastq
  
#
  
# Project:
  
#SBATCH --account=nn9305k
  
#
  
# Wall clock limit:
  
#SBATCH --time=01:00:00
  
#
  
#SBATCH --ntasks=4
  
#
  
# Max memory usage:
  
## A good suggestion here is 4GB
  
#SBATCH --mem-per-cpu=4Gb
  
## Set up job environment
  
  source /cluster/bin/jobsetup
  
module load fastqc
fastqc -t 4 Br_R* Ed_R*
```


5. Move the `html` and `zip` files to `raw_fastqc`
`cd ../raw_fastqc`
`mv ../data/*html .`
`mv ../data/*.zip .`

6. Copy the raw_fastqc folder to Biolinux
In Biolinux 
`scp -r <your_user_name>@abel.uio.no:/work/projects/nn9305k/home/<your_user_name>/Data_pre_processing_Dec14/raw_fastqc .`

7. Go through the html files and discuss.

--------

## Trimmomatic - adapter trimming and removing

We wll use [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) to trim and remove adapter and low quality reads.
This tool is NOT available via `module load` in abel but available at `/work/projects/nn9305k/bin/`. Make sure you know where the adapter sequences are available.

--------

**Task**
1. Move to `trim`folder.
`cd ../trim`

2. Run `trimmomatic-0.36.jar` on Ha_R1.fq.gz file.


<details>
 <summary>Click here for SLURM script (trim.slurm) to trim Ha_R1.fq.gz</summary>
  
  \#!/bin/bash
  
  \#
  
  \# Job name:
  
  \#SBATCH --job-name=trim
  
  \#
  
  \# Project:
  
  \#SBATCH --account=nn9305k
  
  \#
  
  \# Wall clock limit:
  
  \#SBATCH --time=01:00:00
  
  \#
  
  \#SBATCH --ntasks=12
  
  \#
  
  \# Max memory usage:
  
  \## A good suggestion here is 4GB
  
  \#SBATCH --mem-per-cpu=4Gb
  
  \## Set up job environment
  
  source /cluster/bin/jobsetup
  
  java -jar /work/projects/nn9305k/bin/trimmomatic-0.36.jar SE -threads 12 -phred33 ../data/Ha_R1.fq.gz Ha_trim_R1.fq.gz ILLUMINACLIP:/work/projects/nn9305k/db_flatfiles/trimmomatic_adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 CROP:75*
</details>


3. Run `fastqc`on the output fastq files and copy the html and zip to BioLinux and view them in the browser
--------

## FastX-Toolkit

One quick example with [FastX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html)


# Homework

Use trimmomatic to trim/remove adapters and low quality reads in Br_R1.fq.gz and Br_R2.fq.gz (or Ed_R1.fq.gz and Ed_R2.fq.gz)
1. Remember that you are working with paired end data.
2. Change 'SE' to 'PE'
3. There are two input files and four output files
4. Use 'TruSeq3-PE-2.fa' instead of 'TruSeq3-SE.fa'
4. Change MINLEN parameter to 36
5. Use appropriate value for CROP (check the fastqc output and use the correct value)

6. Also, remember to change `#SBATCH --mem-per-cpu=4Gb`to `#SBATCH --mem-per-cpu=12Gb`. This is a bigger job and needs more memory (12Gb instead of 4Gb)

<details>
 <summary>Click here for SLURM script for homework. Please do check this AFTER performing the task yourself</summary>
  
  \#!/bin/bash
  
  \#
  
  \# Job name:
  
  \#SBATCH --job-name=trim
  
  \#
  
  \# Project:
  
  \#SBATCH --account=nn9305k
  
  \#
  
  \# Wall clock limit:
  
  \#SBATCH --time=01:00:00
  
  \#
  
  \#SBATCH --ntasks=12
  
  \#
  
  \# Max memory usage:
  
  \## A good suggestion here is 4GB
  
  \#SBATCH --mem-per-cpu=12Gb
  
  \## Set up job environment
  
  source /cluster/bin/jobsetup
  
  java -jar /work/projects/nn9305k/bin/trimmomatic-0.36.jar PE -threads 12 -phred33 ../data/Br_R1.fq.gz ../data/Br_R2.fq.gz Br_trim_R1.fq.gz Br_trim_R1_UNPAIRED.fq.gz Br_trim_R2.fq.gz Br_trim_R2_UNPAIRED.fq.gz ILLUMINACLIP:/work/projects/nn9305k/db_flatfiles/trimmomatic_adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:150
  
  java -jar /work/projects/nn9305k/bin/trimmomatic-0.36.jar PE -threads 12 -phred33 ../data/Ed_R1.fq.gz ../data/Ed_R2.fq.gz Ed_trim_R1.fq.gz Ed_trim_R1_UNPAIRED.fq.gz Ed_trim_R2.fq.gz Ed_trim_R2_UNPAIRED.fq.gz ILLUMINACLIP:/work/projects/nn9305k/db_flatfiles/trimmomatic_adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 CROP:150
  </details>
