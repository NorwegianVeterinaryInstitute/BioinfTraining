# Visualize your assembly with [IGV](http://software.broadinstitute.org/software/igv/home)

 **IGV**: Integrative Genomics Viewer.

 We will use this softare to evaluate our assemblies, among others to view how the reads map to the assembly.

 We will also add tracks with genome annotation and gaps.

We have assebled reads and annotated the assemblies using Bifrost pipeline. The data for today's practical are to be found in `/project/... `

## 1. Mapping the reads against a reference/assenbly

### 1.1 What is read-mapping

It is simply finding a/the matching locus/area of a read on a sequence. Usually few mismatches are allowed (think about the consequences).

Reads can be mapped as paired or single. If paired is used, then the matching regions are defined by the insert size and the length of each read

NB: note that there is a difference between mapping and aligning what about that? <https://link.springer.com/chapter/10.1007/978-3-319-54064-1_6>

- [ ] discuss (some say aligning reads...) https://en.wikibooks.org/wiki/Next_Generation_Sequencing_(NGS)/Alignment


<html>
<head>
<style>
* {
  box-sizing: border-box;
}

.column {
  float: left;
  width: 50%;
  padding: 5px;
}

/* Clearfix (clear floats) */
.row::after {
  content: "";
  clear: both;
  display: table;
}

</style>

<div class="row">
  <div class="column">
    <a href="https://commons.wikimedia.org/wiki/File:Mapping_Reads.png"><img src="https://upload.wikimedia.org/wikipedia/commons/2/2e/Mapping_Reads.png" width="300">
  </div>
   
  <div class="column">
    <a href="https://www.biostars.org/p/95803/"><img src=" http://www.frontiersin.org/files/Articles/77572/fgene-05-00005-HTML/image_m/fgene-05-00005-g001.jpg" width=300>
  </div>
</div>
<br>


### 1.2 why mapping reads**

1) To evaluate how good an assembly is, statistics are necessary but might not be sufficient to evaluate the quality of an assembly (or to compare different methods used to assemble your reads). statistics might not be sufficient to identify problematic areas.

We want to look at:
- the coverage regularity (ex: some repeated regions might have increased coverage)
- the coverage at the beginning and end of scaffolds - is it good enough?
- are they positions where pilep reads that show polymorphism?
- ...

2) to identify SNPs: some methods to identify SNPs use read mapping of your reads against a reference to detect SNPs

3) Assembly polishing such as `pilon` and `reapr` use mapped reads to identify potential improvement in assemblies.
  > information that can be of use: variation in coverage, wrong read pairs orientation, discrepancy between expected insert size and actual insert optained after read pairs mapping (longer than expected)

 - [ ] other?

### 1.3 Re-map the reads to the assembly optained after pilon.

We use [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) option from [bwa tools](http://bio-bwa.sourceforge.net/) software.


```bash
# If you want to map also unpared reads: we need to concatenate the unpaired reads (otherwise bwa will consider them as PE
cat  <U1.fasq.gz> <U2.fastq.gz>  > <U_reads.fasq>

# 1) We need to index the reference
bwa index <reference.fasta>

# 2) We map the reads (attribute the position of the reads according to the idenx) as PE
bwa mem -t 4 <reference.fasta> <read1.fq.gz> <read2.fq.gz> \
| samtools sort -o <outfile:PE_mapped_sorted.bam> -

# NB: We `sort` directly the mapping by index position
#`-` means that the output of the pipe is used as input

#For unpaired reads
bwa mem -t 4 <references.fasta> <U_reads.fastq.gz> \
| samtools sort -o <outfile:U_mapped_sorted.bam> -

#we need to merge those two files as one
samtools merge <out:all_merged.bam> <in1:_.bam> <in2:.bam> |
samtools sort -o <out:final.bam> <all_merged.bam>
#NB I like to resort to be sure its still sorted

#NB:some software like pilon can bug if the index is not updated - so do:
samtools index <final.bam>
```
- [ ] test sorting mapping if can merge those steps




### 1.4 The [sam/bam file format](https://dash.harvard.edu/bitstream/handle/1/10246875/2723002.pdf?sequence=1)

[wiki-summray](https://en.wikipedia.org/wiki/SAM)

Maybe put a sam file, wont see much with a bam file
- a bam file is a compressed sam file (binary format - computer but not human frendly)

- link to the explanations
- explain difference between formats

-look at the reads
convert bam to sam



### 1.5 looking at the reads with samtools pileup

# 2.Assembly visualisation with IGV

## 2.1 Install IGV in conda

`conda create -n <envname> -c bioconda igv`


## 2.2 creating gap tracks

We use a little python script found in [sequencetools repository](https://github.com/lexnederbragt/sequencetools) to insert gap locations into a file and load it as a track in IGV. This will allow to easily locate the different scaffolds and potentially problematic regions.

- [ ] need to check which version we can use
If you use python 3 : you need to run this version here:
`conda activate <env_with_python3>`
`
`path/scaffoldgap2bed.py -i <assembly.fasta> > <name_gaps>.bed`

## 2.3 Fetch annotation files `.gff`
In Bifrost we annotated the assembly with `Prokka` (using annotations derived from a reference genome)

We will use the `.gff` file created by Prokka: the file containing the annotated `tags` for the IGV visualisation

NB: here need to look that the tags still correspond (corrected VS uncorrected reads) > - [ ] check if this correspond
otherwise will not work

## 2.4 Loading files in [IGV](https://software.broadinstitute.org/software/igv/)

1. Create a [genome File](Creating a .genome File  (Advanced)) as this allows associating tracks to the assembly : `Genomes > create.genome file` - select your assembly file
and fill the different tabs

![picture]()

- [ ] do

2. Load your mapped reads, gap file using: `file > load from file`

[Options and interpretation](http://software.broadinstitute.org/software/igv/PopupMenus#AlignmentTrack), particularly usefull for color alignments


3. To be able to easily reoppen (without re-importing everything you can do: `file > save session` and choose a session name

Now you are ready to navigate and explore your assembly.

![navigation]()
- [ ] do: insert image

You can look at your assembly general or look at specific scaffolds (little meny with all -> here you can choose specifi scaffolds)

You can look how the assmbly looks at the gaps positions

You can look at the coverage - ie. a strange scaffold: repeated gene:

![tRNA]()
- [ ] my scaffold nr 10?/12

Have a look at:
- the coverage
- the gaps positions
- some strange scaffolds?
- PE orientations: in detail how the reads map to your assembly (you will need to zoom a lot)
  - are some PE reads missoriented? reported as having abnormal insert sizes?


---------------------------------------------------------------------
# 3 Additional: assembly evaluation
#### visualizing the insert sizes:
- [ ] - need to make the script functionning if we do that

Additional (outside IGV): we will view the distribution of insert size (estimated from how your reads map to your assembly)
for a subset of the alignments. We modified the script from the Jupyter notebook mentioned in the [Uio course].
/share/inf-biox121/data/assembly/Plot_insertsizes.ipynb : Link to the [original jupyter notebook](https://github.com/lexnederbragt/INF-BIOx121_fall2014_de_novo_assembly/blob/master/practicals/Plot_insertsizes.ipynb)




__________________________________________________________________

# Going further

NB: you can look here at this [Uio course] for more details or if you want to do things slighly differently.
We will reuse several of their [scripts](https://inf-biox121.readthedocs.io/en/2017/Assembly/practicals/Sources.html) so do not be surprized.





??? other things to look at?


For interpretation help: [PE orientations](http://software.broadinstitute.org/software/igv/interpreting_pair_orientations)



[Uio course]:https://inf-biox121.readthedocs.io/en/2017/Assembly/practicals/03_Mapping_reads_to_an_assembly.html
