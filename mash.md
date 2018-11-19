# Screening read files for contaminant reads

This is a short tutorial for the use of Mash screen to identify contaminant reads in read files. Mash screen utilizes a .msh library of refseq genome sketches to identify which organisms the reads in the file originates from.
An explanation of the program can be found here: https://genomeinformatics.github.io/mash-screen/.
A tutorial of the program can be found here: https://mash.readthedocs.io/en/latest/tutorials.html.

## Database
First, you need a database that mash screen can use. You can use the refseq genomes found here: https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh, but note that this is all the genomes in the refseq database, both incomplete and complete genomes from prokaryotes, eukaryotes and viruses. Depending on what you want to use the program for, you might want to create for own database by downloading complete genomes and running sketches on all to create the .msh file for mash screen.
To download the necessary files from refseq, look here: https://github.com/NorwegianVeterinaryInstitute/BioinfTraining/blob/master/genome_downloads.md. 
(Note that this download will take long time and will take a lot of space!). The download will create folders for each genome file.
When the wanted genomes are downloaded, use ```mash sketch``` to create the .msh file:

```
source activate mash

mash sketch -o name_of_file */*.fna.gz
```

Where name_of_file will become the .msh file.

## Mash screen
Unless you expect a difference in read content between R1 and R2, it is probably enough to use the R1 file for the screen.
To use mash screen with default parameters, run the following command:

```
mash screen refseq.genomes.msh fastq_file.fq.gz > mash.out
```

The results can be found in the mash.out file. If you want to run mash screen on multiple read files, it is recommended to create an arrayrun slurm on a hugemem node. See here for more information: https://github.com/NorwegianVeterinaryInstitute/BioinfTraining/blob/master/array_jobs_abel.md.

To hasten the analysis and clean the output, the following settings can be used:

```
mash screen -w -p 4 refseq.genomes.msh fastq_file.fq.gz > mash.out
```
The -p specifies amount of threads to use, while the -w setting reduces redundancy in the resulting file by a "winner takes it all" strategy (see website for more details).

## Results
The mash.out file can be secure copied over to your local computer and imported into excel to view the results. The results will look like this:

identity | shared-hashes | median-multiplicity | p-value | query-ID | query-comment
--------|----------------|---------------------|---------|-----------|--------------
0.99957 | 991/1000 | 26 | 0 | GCF_000841985.1_ViralProj14228_genomic.fna.gz | NC_004313.1 Salmonella phage ST64B, complete genome
0.99957 | 991/1000 | 24 | 0 | GCF_002054545.1_ASM205454v1_genomic.fna.gz | [57 seqs] NZ_MYON01000010.1 Salmonella enterica strain BCW_4905 NODE_10_length_152932_cov_1.77994, whole genome shotgun sequence [...]
0.999522 | 990/1000 | 102 | 0 | GCF_900086185.1_12082_4_85_genomic.fna.gz | [51 seqs] NZ_FLIP01000001.1 Klebsiella pneumoniae strain k1037, whole genome shotgun sequence [...]

- **Identity**: fraction of bases that are shared between the genome and your sequencing reads. Sequencing errors and gaps in coverage will reduced the identity estimate.

- **Shared hashes**: 

- **Median multiplicity**: 

- **P-value**:

- **Query-ID**: accession number of the reference

- **Query-comment**: Name of the reference


## Further analysis
Script for analyzing multiple mash screen reports and identifying files with significant contamination have already been created and can be found here: https://github.com/hkaspersen/misc-scripts.
The script can be found under vi_src/misc_scripts in the project directory.
