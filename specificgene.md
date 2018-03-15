# Specific gene annotation

Goal ist to show:

- specific gene finding
  - two methods
    - assembly, then blast
    - mapping, then local assembly then comparison
  - using asm first may be risky, may transport asm errors into prediction
  - why not just base ourselves on prokka - specific tools give info at greater specificity
  - how ariba actually works
  - importance of databases  
- MLST sequences
  - what is an MLST sequence
  - MLST schemes
  - MLSt prediction from wgs data
  - ariba output
  - phyloviz?
- AMR prediction
  - sort out AMR genes vs mutations
  - ariba output of AMR genes


Practical

- mlst prediction of our test genome
- amr prediction of our test genome

//////////////////////////

In today's session, we will work on finding specific genes in our genomes.
In this case, we have a set of genes that we're interested in, and we want
to see if these are present in our genomes, and if so, we might be interested
in whether there are mutations in them. Example situations are to find the
genes that are used for Multi Locus Sequence Typing (MLST), antibiotic 
resistance and virulence.

There are two types of tools that to this kind of thing:
1. assemble genome, search for reference sequences in genomes
2. align reads to reference sequences, and then extract the 
   genes and variants from the results

As is evident, both approaches require a set of reference sequences that 
we compare to. This also means that these tools can not find sequences
that are not included in the set of sequences that we compare to. This
means that having updated sequences is important.

We will in this session use the program that is called ARIBA to do
both MLST and AMR detection. Virulence finding is also possible with this
program, but we will leave that as an exercise to the participants due to
time constraits. ARIBA uses the second of the two approaches above. ARIBA
is far from the only program that does these things. However, we have
landed on this program for two reasons:
1. It gives very detailed results
2. It enables us to use several databases, and to download the at any
   time current version of the database requested.

## Links to resources

* [ARIBA paper]()
* [ARIBA manual]()
* [PUBMLST, MLST database](https://pubmlst.org/)
* [CARD AMR database]()
* [Resfinder AMR database]()


## Activating ARIBA on abel

ARIBA uses three different programs to do its work, these are `cd-hit`, 
`mummer` and `bowtie2`. To run ARIBA, we to activat them by using `module load`.
ARIBA is a `python` program. To ensure that it runs correctly on abel, we have
had to package it in a specific `python` environment to avoid interactions with
other things on abel.

1. Log in on abel, and go to your project home area and in to your course area
2. Make a directory called `specific_genes` next to the annotation and assembly 
   directories, and go into it.
3. Start the `screen` program, and ask for a qlogin session for 3 hours and 1 cpu.
4. Module load these modules in this specific order: `cd-hit`, `mummer`, `bowtie2`.
5. Type in `source activate ariba211`.

## How to work in this session

We will use two of the homework genomes for this session. We will work two and
two, where each person in the pair will take one genome, and we will later on
in the session merge the results for the two.

* Person A: genome X
* Person B: genome Y


## Multi Locus Sequence Typing

We are first going to do sequence typing for each of our genomes.




Acthman is #1, use this