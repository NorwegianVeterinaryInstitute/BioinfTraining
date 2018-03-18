# Genome Finding and Annotaton

In this session, we will look at genome finding and annotation and how that works.

After assembly, we have a file containing scaffolds. Genome annotation
is the process of figuring out the location of genes in the scaffolds,
and what these genes are. Thus, we can say that genome annotation consists
of two different processes:

1. Gene finding - finding the location of genes in scaffolds
2. Gene annotation - figuring out what the genes are.

We will use a program that is called PROKKA to do genome annotation.
This program merges the results from several other programs to find and
annotate the genes. Some of these programs do both finding and annotation,
others only do one of the two things. The programs used to find rRNAs and
tRNAs, for instance, do both in one. The annotation of proteins is however
a two step process, where we first us a program called `prodigal` to find
the genes, and then `blast` and `hmmer` against several different databases
to figure out what the genes are.

## Links to resources

  * [Presentation made by fellow students](https://docs.google.com/presentation/d/1vKbxuXWcrvvcj5Pi_mE8MzoHJdzVf_g_Fz_QcR5dhoM/edit?usp=sharing)
  * [Link to the PROKKA software](https://github.com/tseemann/prokka/)
  * [How to understand a Genbank record](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html)

## Today's practical

### Filtering a genome

As is evident from our results, it is not unusual to have some very short
contigs. These are usually filtered out. We will today work on the genome
we've worked on in class so far, an ecoli genome. The first thing we will do
is to filter away contigs that are below 500 bp.

The basic spades assembly lives here:

`/work/projects/nn9305k/bioinf_course/compgenomics/assemblies/inclass/spades_careful/scaffolds.fasta`

1. Log in to abel, go to the course directory in your project home area.
2. Create a directory named `annotation` next to the `assembly` directory.
3. Start a new `screen`, and ask for a `qlogin` session for 3 hours and 4 CPUs.
4. The script we're using to do this is a python version 2 script, so we have
   to actiate the python2 module.
5. The name of the script is `filter.contigs.py`. Try running this script from
   the command line with only the option `-h`. Can you figure out which option to use. 
6. Filter the genome using these options: `-c 0`, and input filename and an output 
   filename. Store the genome in this directory.

### Basic PROKKA use

First, we are going to do a first PROKKA run. We will annotate the SPAdes
assembly that we have now filtered.

1. Go read the website mentioned above, to learn about the options to PROKKA.
2. Enable the two following modules: prokka and hmmer/3.1b2 .
3. Type in `prokka -h` and look at the options, specifically, 
   find out what the two mentioned below mean. 
7. Start prokka with the following options: --compliant, --cpus and the assembly above.

It will take a bit of time for the results to come in. Try following along 
on the log file that comes on the screen, and try to see if you can figure out
what is happening.

### How to understand a Genbank file.

Let's have a look at the description of a Genbank record as linked to above. 
You will see that it consists of three different sections:  
* the header region
* the header section
* the Sequence section
   
All Genbank record end with two "//". We can have more than one Genbank record
inside one file.

Click on at least these elements in the sample Genbank file to see what the
components mean: 
* LOCUS
* ACCESSION
* VERSION
* FEATURES
* CDS
* <1..206 (which is in the same line as the first CDS)
* gene

### Looking at the results from PROKKA

Once the results are in, we will look at the results.

1. Have a look at the `.txt` file. What kinds of genes have been found, 
   and how many?
2. Have a look at the `.gbk` produced by PROKKA using `less`.
3. If you have some favorite genes, see if you can find them. If not,
   see if you can find the `gyrA` protein. To search for a text string
   in the program less`less`, type in `/` followed what you're searching 
   for, and press enter.
4. Look at the `.log` file. Can you see how many pseudogenes were detected?
   Note, these are not recorded anywhere else.
5. Have a look at the `.gff` file. Try figure out how the format works.
   Can you find any hypothetical proteins?

### How genes are found and annotated 

Protein coding genes are found using the program `prodigal`. Open
the log file and see if you can find the command line for this program.

If you type in `prodigal -h` on the command line, you can find out 
what each of the options that PROKKA uses for this program.
Try figuring out the following things:
* Will we in our results have partial protein coding genes?
* Will we have protein coding genes that span regions with Ns?
    
### The use of databases to annotate genomes. 

Read the manual for PROKKA again. Try to figure out what databases
it uses to annotate genomes with, and in which order it uses the
various databases.

1. Type in `prokka -h` on the command line. Notice that we have several
   sections of options. Please note, the Organism details section contains
   options that affect the Organism fields in the output files, it does
   not affect what databases are used. We need to look at the Annotations
   details to do that.
2. Use the --listdb option to see which genus databases we have available.
3. Our genome is an E. coli genome. We will now do a new PROKKA run in such
   a way that we get the right name and also use a genus specifc database
   for our search. Use the following options in addition to what we did above:
   --outdir, --locustag, --genus, and --usegenus 
4. Once the run is done, open the Genbank file for the two runs, each in a 
   separate terminal. Find your favorite gene (or `gyrA`) and see if you 
   can see any difference in how they are annotated. 

## Homework

The assemblies for today's homework can be found here:

`/work/projects/nn9305k/bioinf_course/compgenomics/assemblies/homework/spades_covcutoff`

Here, you will find six directories, and inside of those you will find the assemblies,
that is, the scaffolds.fasta files.

Today's homework consist of two parts:

1. Filter the genomes that we have been using so far such that only contigs
   longer than 500 bp remain. Use the assemblathon program to compare the 
   N50, the number of scaffolds, the percentage above 1K and the total length
   to the other genomes.
2. Run prokka on the resulting assemblies. These are all E.coli genomes.
