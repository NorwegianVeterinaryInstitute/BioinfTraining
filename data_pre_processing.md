# Data pre-processing

## Datasets used for this session

Dataset used during this session can be found in the following location within abel:

`/work/projects/nn9305k/tmp/Files_for_Dec14/`

## NB: Replace <your_user_name> with your abel username

Create a new folder called `Data_pre_processing_Dec14`at your home area `/work/projects/nn9305k/home/<your_user_name>/` and move there.

`cd /work/projects/nn9305k/home/<your_user_name>/Data_pre_processing_Dec14`

Create three folder here.

`mkdir data`
`mkdir raw_fastqc`
`mkdir trim`

Move to the `data` folder.

`cd /work/projects/nn9305k/home/<your_user_name>/Data_pre_processing_Dec14/data`

Type the following command to link the files (not copy):

`ln -s /work/projects/nn9305k/tmp/Files_for_Dec14/*fq.gz .`


## Fastq quality check

We will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to check the quality of raw sequenced data 

--------

**Task**
1. Run `fastqc` on the file `Ha_R1.fq.gz` (in the login node).
2. type in the following (don't include the `$`).

```
$ fastqc Ha_R1.fq.gz
```

3. Try to find help for `fastqc` and discuss what flags one can use to process multiple samples
4. Use `SLURM` to process the other four files


<details>
 <summary>Click here for SLURM script (fastq.slurm) to process four fastq files</summary>
  
  \#!/bin/bash
  
  \#
  
  \# Job name:
  
  \#SBATCH --job-name=raw_fastq
  
  \#
  
  \# Project:
  
  \#SBATCH --account=nn9305k
  
  \#
  
  \# Wall clock limit:
  
  \#SBATCH --time=01:00:00
  
  \#
  
  \#SBATCH --ntasks=4
  
  \#
  
  \# Max memory usage:
  
  \## A good suggestion here is 4GB
  
  \#SBATCH --mem-per-cpu=4Gb
  
  \## Set up job environment
  
  source /cluster/bin/jobsetup
  
  module load fastqc
  
  fastqc -t 4 Br_R* Ed_R*
 </details>


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
  
  java -jar /work/projects/nn9305k/bin/trimmomatic-0.36.jar SE -threads 12 -phred33 ../data/Ha_R1.fq.gz Ha_trim_R1.fq.gz ILLUMINACLIP:/work/projects/nn9305k/db_flatfiles/trimmomatic_adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 CROP:75
  </details>



3. Run `fastqc`on the output fastq files and copy the html and zip to BioLinux and view them in the browser
--------

## FastX-Toolkit

One quick example with [FastX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html)

