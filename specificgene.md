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

### Learning about the PUBMLST schemes and sequences

We will first explore the website and the MLST schemes a bit.

1. Go to the PUBMLST website mentioned above. 
2. Click on `Download MLST definitions`
3. Scroll down, and find the two Escherichia coli schemes. 
4. Question: Can you figure out how many genes are shared between the 
   two schemes?
5. Click on `profiles` under one of the schemes. This is what a MLST
   scheme looks like.
6. Now go back to the previous page, and click on one of the genes in
   the #1 scheme. You will now get a long list of fasta sequences on
   your screen. As you can see, each of them have a number behind them.
   These numbers correspond to the numbering in the MLST scheme.
   
7. OPTIONAL: download this file by right-clicking on the window and
   selecting `Save as...`.
8. From within the VM, go to the menu, and click `Applications - Bioinformatics - Jalview`.
9. Use the `File` menu and load in the file you just downloaded via `Input alignment - From file`. 
   Remember, yousaved it as text, so you have to ask it to show `All files`.
10. Go to the `Color` menu and select `Percentage identity`. You will now see
   the differences between the sequences.
11. Close `Jalview`.

### Running MLST ARIBA


1. Go read the ARIBA manual and have a look at the presentation made by your
   fellow students. Try to figure out the order the following commands have
   to be run in:
   * ariba pubmlstget "YOU NEED TO WRITE SOMETHING HERE" mlstdb
   * ariba pubmlstspecies
   * ariba run mlstdb/ref\_db reads\_1.fq reads\_2.fq ariba\_out
2. Figure out where the reads for your isolate is.
3. Run the commands in the order they should be in.
4. When you are done, you should have two files for your isolate.
   Open each of them with `less` and have a look at them.
5. Question: are your isolates from the same or different sequence
   types? Note, you can find explanations of the files on the 
   ARIBA wiki website.

!!!! something here about merging results, and visualization with phyloviz


## Antibiotic resistance finding

### AMR databases

!!!! something with showing them the website for CARD and Resfinder, and
showing them that one has one thing, and the other has the ther things.

### Running AMR ARIBA

1. Go read the ARIBA manual and have a look at the presentation made by your
   fellow students. Try to figure out the order the following commands have
   to be run in:
   * ariba getref card out.card
   * ariba run out.card.prepareref reads1.fastq reads2.fastq out.run
   * ariba prepareref -f out.card.fa -m out.card.tsv out.card.prepareref
2. Figure out where the reads for your isolate is.
3. Run the commands in the order they should be in.
4. Go into the result directory. Have a look at the `report.tsv` file. Answer
   the following questions. Note. the information found for the `run` and
   `prepareref` commands can be useful for finding answers.
   * Which AMR genes were found?
   * Are these presence/absence AMR genes or mutation AMR genes?
   * If there are any mutation AMR genes, which mutations do you have?
5. !!!! copy neighbhorfile, do summary, do fandagno.
   