# Comparative genomics

Teacher: Karin Lagesen


## Format for this course

The course takes place on Mondays, at 9am, in Resepten. Remember to bring your
laptops. We will be running analyses on abel, and using tools within the 
virtual machine for evaluation of the results.

Each week we will do the following:
  * Introduction to today's theme
  * Practical exercise on a small dataset
  * Homework, where you will perform the same types of exercise on a different 
  given dataset
  
Most weeks, two students will be asked to review some more material regarding
how something works, and then prepare a short (10-15 minutes) presentation on
the subject.


## Practicals

For the first part of this course, we will follow the assembly exercises
that are part of the INF-BIOx121 course at the University of Oslo.

[Link to the relevant pages](https://github.com/karinlag/INF-BIOx121/tree/2017/Assembly/practicals) 

You will also find a cheat list for commands at the bottom of this page.


### 2018-02-05 and  2018-02-12 

#### Preparatory work

Make sure that access to abel works, and that the jupyter notebook (anaconda)
is installed on the virtual machine.

#### Today's practical

Program for today:
  * Short lecture on how assembly works
  * Experiment with de Bruijn graphs, using the DeBruijnGraph.ipynb notebook
  * Go through most of the Velvet exercise from the INF-BIO site
  
The dataset that we will use for the in-class exercises can be found here:

`/work/projects/nn9305k/bioinf_course/compgenomics/rawdata/inclass`

Note, for the INF-BIO course, the students were working on a non-abel computer.
To make these exercises work, we have to do some things to make things work.

1. Log onto abel
2. Check which login node you are on (`hostname`)
3. Start the `screen` program
4. Go to your project home directory
5. Create a course directory, call it `2018_compgenome`, and go into it
6. Ask for a qlogin session:  
`qlogin --account=nn9305k --time=03:00:00 --ntasks=2 --mem-per-cpu=4G`
7. To run the velvet program, we have to use the `module` system to load it:  
`module load velvet`

We can now work through the exercises more or less as described. Note,
we only have paired end data, and will only be working with that data.

We will be using a google spreadsheet to record our assembly data:

[Velvet k-mer test](https://docs.google.com/spreadsheets/d/1mvIV0jenKBWGxIVyHTMe2Stb2OkNuLRv0iqpliC3URY/edit?usp=sharing)


When it comes to downloading and running notebooks:
1. On github, click on the filename so that the notebook is displayed.
2. Click on `raw` in the upper corner
3. Copy the URL in the address field
4. Go to a terminal on your local computer, and type in `wget URL`
5. To run the notebook: `jupyter notebook filename.ipynb`


#### Homework
 
There are six genomes in the directory 

`/work/projects/nn9305k/bioinf_course/compgenomics/rawdata/homework`

You will work with a partner. Each team will pick three k-mer sizes and use 
velvet to assemble the six genomes. Each person in the team is responsible
for three genomes each.

[Record your results here](https://docs.google.com/spreadsheets/d/124Eb6IQ44coSKMH0kRLU18AJ5FZ7-ijwsxqf1NsC9Ys/edit?usp=sharing)


### 2018-02-17 

#### Today's practical

Today's practical consists of different sections:
  * first, we will assemble the test genome with SPAdes, using two different
  options
  * second, we will run the program QUAST to evaluate our assemblies
  * third, we will spend some time understanding the QUAST output

We are again working on abel, in the `2018_compgenome` directory we created
last time. We want to be working inside of a `qlogin`, since some of the
things we will do might take a bit of time.

Note: to use SPAdes and QUAST, we have to use `module load`. There
are some extra files we need for QUAST, these are:
   * A genome sequence: 
   `/work/projects/nn9305k/genome_references/genomes/ecoli/GCF_000005845.2_ASM584v2_genomic.fna `
   * A genome annotation file, in `GFF` format:   
    `/work/projects/nn9305k/genome_references/genomes/ecoli/GCF_000005845.2_ASM584v2_genomic.gff`

We will examine the QUAST results using a browser. To do that, we need to 
transfer the QUAST result directory to our local computer. We can do that using 
this command line on our local vm:

`scp -r username@abel.uio.no:/full/path/to/directory .`

1. In your webbrowser, go to the course pages listed above. Find the SPAdes module. 
   * Create a new directory _next_ to the velvet directory, called spades, and 
     go into it
   * Figure out what the SPAdes command line looks like. 
   * Create two assemblies from the inclass data, one using the `--careful` 
     option, and one that does not use this option. Name their output
     directories `spades_wo` and `spades_careful`
   * Google for the SPAdes manual and try to figure out two things:
     * Which k-mers are used
     * What does --careful do?
   * Go into the results directory, and and look at the results
   * Run the assemblathon script on the two assemblies, and do a quick comparison between the
     * velvet assembly
     * the SPAdes assembly
     * the SPAdes assembly with `--careful`
     
2. In your webbrowser, go to the course pages listed above. Find the module
   where we are comparing to a reference. 
   * Have a look at the gff file that contains the genome annotation using 
   `less`. What does it contain? What's the format like?
   * Using the webpage from the course, figure out what the command line for
   QUAST should look like. In your comparison, include the velvet assembly,
   the SPAdes assembly, and the SPAdes careful assembly. IMPORTANT! You
   should include the `--scaffold` option.
   * Run the QUAST comparison.
   
3. While QUAST is running, have a look at the 
   [QUAST paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624806/) and 
   the [QUAST manual](http://quast.bioinf.spbau.ru/manual.html).
   Figure out the following things:
   * What basic statistics does QUAST produce? 
   * What is a missassembly according to QUAST?
   * What is the Duplication ratio, number of mismatches per 100k, and number of genes?
   * What are the Nx, NAx, NGx and so on values?

4. Transfer the resulting QUAST output directory to your local computer, and
   open the `report.html` file using firefox.
   * Have a look at the basic statistics
     * Which assembly has the highest N50?
     * Which assembly has the most found genes?
     * Which genome has the fewest misassemblies?
   * Open the Icarus browser
     * Select the first 300 000 bases of the reference genome
     by typing in the numbers `1` and `300000` in the boxes above
     * Click on the first contig in the velvet assembly. Have a
       look at the box on the right. It seems like this contig is
       split up, but it kind of isn't. What's going on here?
       Use the arrows that is shown in the section below to understand
       what is going on.
     * The velvet assembly is first, and it has a red contig in the 
       unbroken assembly. Figure out what the difference between
       the two assemblies here, by examining the boxes on the right.
     * The two SPAdes assemblies look very similar. Can you find some
       regions where there is a difference between the non-careful and
       the careful assembly?  
   
   
### Homework

Remember to use `qlogin` when doing the exercises.

1. Redo the velvet assemblies using the command in the cheat list. We did get 
   to the place where I should have introduced the `shortPaired` option last
   week, which is why you didn't use it for your assemblies. This option helps
   a lot, which is why we should use it. 

2. Do SPAdes assemblies for all the genomes that we gave you. Remember to use 
   --careful
   
3. Use QUAST and compare each of your SPAdes assemblies to the velvet assembly
   of the same genome. Use the same reference genome and annotation as
   in the exercise. Which assembly is `best` and why?
   

   
## Command line cheat list

Note: anything in CAPITAL LETTERS should be replaced with something, usually
a file name, an output directory name, a fasta file name or something similar.
You can figure out what the various options do by googling for the manual,
by using the INFBIO course pages linked to above, or (quite often) typing
in the name of the program, without any options, and pressing enter.

Note, in the command below, long commands will be broken up with a `\`. If you
write in the commands in one long line, you do not include this slash.

#### qlogin

XX should be replaced with how many hours you want, CPUS with the number
of cpus you want to use. Remember, if you ask for more than one cpu, you
should actually use those CPUs when running commands. I.e. for SPAdes
for instance, you should specify the `-t` option. You are also likely to want
to use `screen` before requesting a `qlogin` session.

```
qlogin --account=nn9305k --time=XX:00:00 --ntasks=CPUS --mem-per-cpu=4G
```
#### Velvet

```
velveth ASM_NAME VALUE_OF_K \  
-shortPaired -fastq -separate \  
PATH/TO/READ_1.FASTQ \  
PATH/TO/READ_2.FASTQ \  

velvetg ASM_NAME -exp_cov auto -cov_cutoff auto  

```

#### SPADdes

```
spades.py -t CPUS --careful --pe1-1 PATH/TO/READ_1.FASTQ -pe1-2 PATH/TO/READ_2.FASTQ \  
-o OUTPUTDIRECTORY > LOGFILENAME.LOG 2>&1 

```
   
#### QUAST

Note, you can compare as many assemblies as you want, from 1 to hundreds. 
You need to add them sequentially to the command line, and remember to 
name them prperly inside the `-l` option.
 
```
quast.py -t CPUS -o OUTPUTDIRECTORY \ 
-R /work/projects/nn9305k/genome_references/genomes/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \ 
-G /work/projects/nn9305k/genome_references/genomes/ecoli/GCF_000005845.2_ASM584v2_genomic.gff \ 
--scaffolds PATH/TO/ASSEMBLY_1.FSA PATH/TO/ASSEMBLY_2.FSA PATH/TO/ASSEMBLY_3.FSA \ 
-l "ASSEMBLY_1, ASSEMBLY_3, ASSEMBLY_3" > QUAST.LOG 2>&1
```
