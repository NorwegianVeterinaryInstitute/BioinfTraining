# Bioinformatics training: transcriptomics

## Protocol
We will follow the protocol described in [Tophat2 bioinformatic protocol](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3334321/) published in [Nature protocol, 2012](https://www.nature.com/articles/nprot.2012.016) for this course. 

## Day 1
* Introduction of transcriptomic analyses 
* Discusssion regarding reference-based and de novo approaches
* How to use tophat v2 - based on the Nature protocol (2012) paper
* Download genome reference and create an index using bowtie2

* in /work/projects/nn9305k/home/<username>/transcriptome/ref
```
$ ln ../../../../bioinf_course/transcriptomics/ref/Dm.BDGP6.dna.toplevel.fa .
$ ls ../../../../bioinf_course/transcriptomics/ref/Dm.BDGP6.91.gtf .
$ bowtie2-build Dm.BDGP6.dna.toplevel.fa Dm_BDGP6_genome
```
  
* USAGE
```
$ bowtie2-build <Reference_fasta_file_name> <bowtie2-build_ref_index_name_that_you_will_use_later>
```

# To do
* Align the given reads to the Drosophila genome using tophat v2

## Day 2
* Discuss results from tophat v2 alignemnt
