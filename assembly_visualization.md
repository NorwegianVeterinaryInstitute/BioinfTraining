# Visualize your assembly with [IGV](http://software.broadinstitute.org/software/igv/home)

 **IGV**: Integrative Genomics Viewer.

 We will use this software to evaluate our assemblies, among others to view how the reads map to the assembly.

 We will also add tracks with genome annotation and gaps.

We have assembled reads and annotated using Bifrost pipeline. The data for today's practical are found in `/projects/nn9305k/Bioinfo_mapping_training`

## 1. Mapping the reads against a reference/assembly

### 1.1 What is read-mapping

It is simply finding a/the matching locus/area of a read on a sequence. Usually few mismatches are allowed (think about the consequences).

Reads can be mapped as paired or single. If paired is used, then the matching regions are defined by the insert size and the length of each read

<center> <a href="https://commons.wikimedia.org/wiki/File:Mapping_Reads.png"> <img src="https://upload.wikimedia.org/wikipedia/commons/2/2e/Mapping_Reads.png" width="400">
<br>


### 1.2 why mapping reads

1) Statistics are necessary but might not be sufficient to evaluate the quality of an assembly (or to compare different methods used to assemble your reads). Read mapping can help identifying problematic areas.

We want to look at:
- the coverage regularity (ex: some repeated regions might have increased coverage)
- the coverage at the beginning and end of scaffolds - is it good enough?
- are they positions where pileup of reads show polymorphism?
- ...

2) to identify SNPs: some methods use reads mapping against a reference to identify and type SNPs

3) Assembly polishing software such as `pilon` and `reapr` use mapped reads to identify potential improvement in assemblies. It is good to visualize what information they are using ...
  > eg. variation in coverage, wrong read pairs orientation, discrepancy between expected insert size and actual insert size obtained from read mapping (ie. longer than expected)

 - [ ] other?

### 1.3 Re-map the reads to the assembly

We use [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) option from [bwa tools](http://bio-bwa.sourceforge.net/) software.

You will use the trimmed-filtered reads in: `/projects/nn9305k/Bioinfo_mapping_training/bifrost_run/bbduk_trimmed` and the assembly: `/projects/nn9305k/Bioinfo_mapping_training/bifrost_run/prokka/MiSeq_Ecoli_MG1655_50x.fna`

We could use the assembly output from `Pilon` **BUT** then we would have to create aliases between labels used by `Prokka` and by `Pilon` to be able to visualize the assembly with IGV. `Prokka` assembly is identical to `Pilon` assembly, only labels have been changed.

In your home make a directory for today's work
and a folder called `mapping` where you will copy the input files

```bash
# softlink the required input files in the mapping folder:
ln -s <path>/bifrost_run/bbduk_trimmed/*.fq.gz .
ln -s <path>/bifrost_run/pilon/MiSeq_Ecoli_MG1655_50x_pilon_spades.fasta .

# we will use Bifrost with the command line
# As bwa and samtools are used to map reads for Pilon ...

qlogin --account=nn9305k --time=00:30:00 --ntasks=4 --mem-per-cpu=4G
source activate bifrost

# 1) We need to index the reference
bwa index <reference.fasta>

# 2) We map the reads (attribute the position of the reads according to the index) as PE
bwa mem -t 4 <reference.fasta> <in1:read1.fq.gz> <in2:read2.fq.gz> \
| samtools sort -o <out:PE_mapped_sorted.bam> -

# NB: We `sort` directly the mapping by index position
#`-` means that the output of the pipe is used as input in samtools

#For unpaired reads
bwa mem -t 4 <references.fasta> <in:U_reads.fastq.gz> \
| samtools sort -o <out:U_mapped_sorted.bam> -

# We need to merge those two files as one
samtools merge <out:all_merged.bam> <in1:_.bam> <in2:.bam>

# To be sure reads are still sorted: we resort
samtools sort -o <out:final_all_merged.bam> <in:all_merged.bam>

#NB:Some software like Pilon need that we update the index (you can do:). We will need index after
samtools index <final_all_merged.bam>
```
### 1.4 The [sam/bam file format](https://samtools.github.io/hts-specs/SAMv1.pdf)

You can also have a look at [Samtools article](https://academic.oup.com/bioinformatics/article/25/16/2078/204688) and at [samtools manual](http://www.htslib.org/doc/samtools.html)

`.bam` files (position indexed mapped-reads) are in a compressed binary format. We need to transform the `.bam` (to a `.sam` file) to be able to see how indexed mapped-reads are represented in the file.  

To decompress: chose f.eks. `PE_mapped_sorted.bam` that we did in the first step:

`samtools view -h -o <out.sam> <in.bam>`

Look at your `.sam` file with:

 `less <filename.sam>`

### 1.5 looking at how the reads maps against the reference with [samtools](http://www.htslib.org/doc/samtools.html) pileup module

```bash
# To see the list of nodes:
grep "^>" <reference.fasta>
# look at the pileup in one node
samtools mpileup <in:final_all_merged.bam> --max-depth 0 --fasta-ref <reference.fasta> -r <ex:NODE_76_length_637_cov_0.071429_pilon> | less
```
- [ ] here if you want to do something different, that's up to you ...
I just tried

# 2.Viewing assembly and mapping with IGV

## 2.1 Install IGV in conda
At your pc:

```bash
conda create -n <envname> -c bioconda igv
conda activate <IGV:envname>
conda install -c anaconda biopython
```

Transfer the file: containing the assembly, the .bam file and some useful scripts to your pc, as well as files generated to visualize with IGV. How those files have been created is explained below.

```bash
scp <username>@abel.uio.no:/projects/nn9305k/Bioinfo_mapping_training/for_igv.tar.gz <destination_path_pc>
# Then:
tar -xzvf for_igv.tar.gz
```
- [ ] this contains all files necessary to run IGV (renamed because I had to redo everyting because of labels) - then can also fetch the files, but that might be longer...

## 2.2 creating gap tracks

We use a little python script from [sequencetools repository](https://github.com/lexnederbragt/sequencetools) to insert gap locations into a file and load it as a track in IGV. This will allow to easily locate the different scaffolds and potentially problematic regions. This is how we generated the file:

```bash
conda activate <IGV:envname>
python path/scaffoldgap2bed_py3.py -i <assembly.fasta> > <name_gaps>.bed
```

## 2.3 Fetch annotation files `.gff`

In `Bifrost` we annotated the assembly with `Prokka` (using annotations derived from a reference genome) `.gff` file. It contains the annotated gene locus `tags` that we use for visualization with IGV

The file is here: `/projects/nn9305k/Bioinfo_mapping_training/bifrost_run/prokka`

## 2.4 Loading files in [IGV](https://software.broadinstitute.org/software/igv/)

![IGV](./figures/IGV.png)

1. Create a [genome File](Creating a .genome File  (Advanced) as this allows associating tracks to the assembly : `Genomes > create.genome file`
 > use the menu to select your assembly file `.fasta`and the gene file: `.gff`

2. Load your `mapped reads` and the `gap file` using: `file > load from file`

3. To be able to easily re-open (without re-importing everything you can do: `file > save session`

Now you are ready to navigate and explore your assembly.

> Find a gap.
> - NB: To zoom while staying centered on the gap: click with the mouse at the gap position

You can look here for [Options and interpretation](http://software.broadinstitute.org/software/igv/PopupMenus#AlignmentTrack),
and here: [PE orientations](http://software.broadinstitute.org/software/igv/interpreting_pair_orientations).

Have a look at:
- coverage
- gaps positions
- some strange scaffolds?
- PE orientations: in detail how the reads map to your assembly (you will need to zoom a lot)
- are some PE reads miss-oriented? reported as having abnormal insert sizes?

---------------------------------------------------------------------
# 3 Additional: assembly evaluation

- [ ] do
#### visualizing the insert sizes:
- [ ] - need to make the script functioning if we do that

Additional (outside IGV): we will view the distribution of insert size (estimated from how your reads map to your assembly)
for a subset of the alignments. We modified the script from the Jupyter notebook mentioned in the [Uio course].
/share/inf-biox121/data/assembly/Plot_insertsizes.ipynb : Link to the [original jupyter notebook](https://github.com/lexnederbragt/INF-BIOx121_fall2014_de_novo_assembly/blob/master/practicals/Plot_insertsizes.ipynb)

otherlink <https://github.com/lexnederbragt/INF-BIOx121/blob/2017/Assembly/practicals/Plot_insertsizes.ipynb>

_________________________________________________________________

# Going further

You can look here at the [Uio course] for more details or if you want to do things slightly differently. We reuse some of their [scripts](https://inf-biox121.readthedocs.io/en/2017/Assembly/practicals/Sources.html).

[Uio course]:https://inf-biox121.readthedocs.io/en/2017/Assembly/practicals/03_Mapping_reads_to_an_assembly.html
