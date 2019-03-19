# Visualize your assembly with IGV: Integrative Genomics Viewer 

For [IGV installation](./assembly_visualization.md#igv-install): look down the page

NB: you can look here at this [Uio course] for more details or if you want to do things slighly differently. 
We actually reused several of their scripts so do not be surprized. 

> We assume that you generated you assemblies using Bifrost
... and you want to have a look at your assembly and how the reads are mapped to your assembly, look at your coverage.... and more generaly evaluate your assembly visually. 

** To do so, we need to do several things:**

1) re-map the reads to the assembly optained after pilon. 
 > We will use PE mapping and `bwa` option `mem` to map your reads to the scaffolds producted by pilon 

2) we will use a little python script to insert gap locations into a file. This will allows use to add gap locations
asa track in IGV and therefore will allows you to easily locate the different scaffolds and potentially problematic regions.

If you use python 3 : you need to run this version here: (here just need to modify the print)
(just copy the script and save it in a text file that you name scaffoldgap2bed.py"
scaffoldgap2bed.py -i ASSEMBLY.FASTA >gaps.bed
> python version? 

3) In Bifrost you generated annotations of your assembly with `Prokka`
(by using annotations derived from a reference genome you chosed) 
We will use the `.gff` file created by Prokka: the file contening the annotated `tags` for the IGV visualisation

NB: here need to look that the tags still correspond (corrected VS uncorrected reads) > - [ ] check if this correspond
otherwise will not work


4) Load everything into IGV 


5) Additional (outside IGV): we will view the distribution of insert size (estimated from how your reads map to your assembly)
for a subset of the alignments. We modified the script from the Jupyter notebook mentioned in the [Uio course].
/share/inf-biox121/data/assembly/Plot_insertsizes.ipynb
 
[Options and interpretation in IGV](ttp://software.broadinstitute.org/software/igv/PopupMenus#AlignmentTrack)


# The code part

PS: look at the comments to understand what we do at eatch step

```bash
#we map our reads to the optained assembly
## first creating the index for mapping
bwa index <my_scaffolds.fasta>

## second we mapped PE and U reads separately (NB: 2 Unpaired reads files after spades)
## and sort the mapped reads by position (index) using samtools

bwa mem -t 4 <my_scaffolds.fasta> <R1_reads.fasq.gz> <R1_reads.fasq.gz> \
| samtools sort -o <PE_mapped_sorted.bam> -

bwa mem -t 4 <my_scaffolds.fasta> <U1_reads.fasq.gz> <U2_reads.fasq.gz> \
| samtools sort -o <U_mapped_sorted.bam> -

## We only want one final file (both with PE and U mapped reads) - so we merged the .bam files
## Maybe can make it in one step have to test
samtools merge <out:PE&U_all_merged.bam> <in1:_.bam> <in2:.bam> 

## to be sure everything is still sorted (it should be but I do not take risks) 
samtools sort -o <out:final.bam> <PE&U_all_merged.bam>

## refresh the index of the bam file
samtools index <final.bam>


#creating gap track with same names
conda activate python3
python scaffoldgap2bed.py -i <my_scaffolds.fasta> > <mygaps.bed>
conda deactivate

#visualizing the insert sizes:
- do [ ]

```
# Results interpreation
## Loading your files into IGV

Lauching IGV 
```bash
conda activate IGV
igv
```

1. Load your assembly as a reference genome: `Genomes > create.genome file` - select your assembly file 
and fill the different files

2. Load your mapped reads, gap file using: `file > load from file`

3. To be able to easily reoppen (without re-importing everything you can do: `file > save session` 
and choose a session name

Now you are ready to navigate and explore your assembly.

You can look at your assembly general or look at specific scaffolds (little meny with all -> here you can choose specifi scaffolds) 

![navigation]() 
- [ ] do: insert image

You can look how the assmbly looks at the gaps positions

You can look at the coverage - ie. a strange scaffold: repeated gene: 

![tRNA]() 
- [ ] my scaffold nr 10?/12

You can have a look at the PE orientations: in detail how the reads map to your assembly (you will need to zoom a lot) 

??? other things to look at? 


For interpretation help: [PE orientations](http://software.broadinstitute.org/software/igv/interpreting_pair_orientations)

# [Installation of IGV](#igv-install)

> we do that in conda: `conda create -n myenv scipy`

`conda create -n IGV igv` 

[Uio course]:https://inf-biox121.readthedocs.io/en/2017/Assembly/practicals/03_Mapping_reads_to_an_assembly.html
