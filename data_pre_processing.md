# Data pre-processing

## Datasets used for this session

Dataset used during this session can be found in the following location within abel:

`/work/projects/nn9305k/tmp/Files_for_Dec14/`


Create a new folder called `Data_pre_processing_Dec14`at your home area `/work/projects/nn9305k/home/<your_user_name>/`

Create a folder `data` within the above new folder and navigate there.

Type the following command to link the files (not copy):

`ln -s /work/projects/nn9305k/tmp/Files_for_Dec14/*fq.gz .`

## Fastq quality check

We will use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to check the quality of raw sequenced data 

--------

**Task**
1. Run `fastqc` on the file `Ha_R1.fq.gz`.
2. type in the following (don't include the `$`).

```
$ fastqc Ha_R1.fq.gz
```

3. Try to find help for `fastqc` and discuss what flags one can use to process multiple samples
4. Use `SLURM` to process the other four files
5. Move the `html` and `zip` files to BioLinux and dicuss the results
--------

## Trimmomatic - adapter trimming and removing

We wll use [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) to trim and remove adapter and low quality reads.
This tool is NOT available via `module load` in abel but available at `/work/projects/nn9305k/bin/`. Make sure you know where the adapter sequences are available.

--------

**Task**
1. Run `trimmomatic-0.36.jar` on the three sets of fastq files.
2. Run `fastqc`on the output fastq files
--------

## FastX-Toolkit

One quick example with [FastX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/index.html)
--------

**Task**
1. Run `trimmomatic-0.36.jar` on the three sets of fastq files.
2. Run `fastqc`on the output fastq files
--------
