# Variant calling with SNIPPY

## Manuals

  * [SNIPPY manual](https://github.com/tseemann/snippy)
  * [BWA manual](http://bio-bwa.sourceforge.net/bwa.shtml)
  * [samtools](http://www.htslib.org/doc/samtools.html)
  * [FreeBayes](https://github.com/ekg/freebayes)
  * [bcftools](http://www.htslib.org/doc/bcftools.html) 

## Getting reference data

  * [NCBI Ecoli ST38](https://www.ncbi.nlm.nih.gov/search/all/?term=escherichia%20coli%20ST38)

## Practicals - single genome:

  1. Create a directory to work in, and go into it
  2. Copy data from reference site: 
  ``` bash
  cp -r /projects/nn9305k/bioinf_course/variantcalling/snippy .
  ```
  3. Look to see what contents are in the various directories (analysis is pre-cooked)
  4. Run a simple analysis, with one isolate against the reference
      * Create a screen
      * Set up a qlogin
      * activate snippy 4, with conda (conda activate snippy4)
      * Run a single genome
  5. While it is running, try to identify what commands that is being run. We will 
     dissect them one by one. 
  6. Look at output files, and see what they are - compare to SNIPPY webpage.
  7. How to look at bam files:
      * samtools view snps.bam | less
          * find name of contig, second column
      * samtools tview snsps.bam reference/ref.fa
      * use `q` to get a search field, input contig name
  
## Practicals - run with multiple genomes

SNIPPY can also run multiple genomes against the same reference. This is done with the
`snippy-multi` option. As input to that, we create a list of genomes as seen in `list.txt`.

We can use this as input to a snippy-multi run. A sample result from such a run can be found
in the `analysis` directory. 




