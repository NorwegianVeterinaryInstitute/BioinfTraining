# Specific gene annotation

In today's session, we will work on finding specific genes in our genomes.
In this case, we have a set of genes that we're interested in, and we want
to see if these are present in our genomes. If so, we might be interested
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
time constraits. 

ARIBA uses the second of the two approaches above. ARIBA
is far from the only program that does these things. However, we have
landed on this program for two reasons:
1. It gives very detailed results
2. It enables us to use several databases, and to download the at any
   time current version of the database requested.

## Links to resources

* [Student presentation](https://docs.google.com/presentation/d/1FfiKvkPrwsNn9Dpqz-TucZH4VH6morUMHig2JgqzAq8/edit?usp=sharing)
* [ARIBA paper](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000131)
* [ARIBA manual](https://github.com/sanger-pathogens/ariba/wiki)
* [PUBMLST - MLST database](https://pubmlst.org/)
* [CARD AMR database](https://card.mcmaster.ca/home)
* [ResFinder AMR database](https://cge.cbs.dtu.dk//services/ResFinder/)
* [Phandango site](http://jameshadfield.github.io/phandango/#/) 


## Activating ARIBA on abel

ARIBA uses three different programs to do its work, these are `cd-hit`, 
`mummer` and `bowtie2`. To run ARIBA, we to activate them by using `module load`.
ARIBA is a `python` program. To ensure that it runs correctly on abel, we have
had to package it in a specific `python` environment to avoid interactions with
other things on abel.

1. Log in on abel, and go to your project home area and in to your course area. 
2. Make a directory called `specific_genes` next to the annotation and assembly 
   directories, and go into it.
3. Start the `screen` program, and ask for a qlogin session for 3 hours and 4 cpus.
4. Module load these modules in this specific order: `cd-hit`, `mummer`, `bowtie2`.
5. Type in `source activate ariba211`. Note, this step might take a minute or two!
6. Explore the options to ariba
   * ariba -h
   * ariba getref -h
   * ariba pubmlstget -h
   * ariba run -h
   * ariba summary -h

## How to work in this session

We will use two of the homework genomes for this session. We will work two and
two, where each person in the pair will take one genome. 
* Person A: pick F27
* Person B: pick S19

## Multi Locus Sequence Typing

We are first going to do sequence typing for each of our genomes.

### Learning about the PUBMLST schemes and sequences

We will first explore the website and the MLST schemes a bit.

1. Go to the PUBMLST website mentioned above. 
2. Click on `Download MLST definitions`.
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

1. Go read the ARIBA manual and also have a look at the presentation made by your
   fellow students. Try to figure out the order the following commands have
   to be run in, and what you have to fill in. Note. Our reads are E. colis,
   and you should use the first schema.
   * ariba pubmlstget "YOU NEED TO WRITE SOMETHING HERE" mlstdb
   * ariba pubmlstspecies
   * ariba run --threads 4 mlstdb/ref\_db READ1.FQ READ2.FQ OUTPUTDIRNAME
2. Figure out where the reads for your isolate is.
3. Run the commands in the order they should be in.
4. When you are done, you have a new directory with results in it.
   The main files are `report.tsv`, `mlst_report.tsv` and `mlst_report.details.tsv`.
   Open each of them with `less` and have a look at them. Note, you can find 
   explanations of the files on the ARIBA wiki website.
5. Question: are your isolates from the same or different sequence
   types? 
6. Note, the results found in `mlst_report.tsv` can, when aggregated,
   be used as input for tools such as `Phyloviz`, which is a java
   vizualisation program that we won't go into here.

## Antibiotic resistance finding

### AMR databases

It is important to know that different databases contain different things.
Considering how specific gene finding works, we can only find things if
those things are in the database that we're using for our finding. Thus,
it is imperative for you to know that the things you are interested in are
in the database that you are using.

ARIBA allows us to use several databases for AMR finding. Go have a look 
at the ARIBA wiki page to see which databases are available.

#### Exploring ResFinder

1. Go to the ResFinder webpage
2. Scroll down, and click on `The ResFinder database download site`.
3. Find the ResFinder database in the list, and click on `notes.txt`.
   This will download a text file to your computer. Click on it.
3. This file should now open in a text editor window. Have a look at it.
   If you choose to search with the ResFinder database, these are
   the genes that you will be searching for.
4. Have a look at the Quinolone resistance section. Can you find the
   gyrA gene there?

#### Exploring CARD

1. Go to the CARD database webpage
2. Click on `Browse`. Then, look under `Antibiotic Resistance Ontology`, 
   click on `ARO Index`.
3. Set the `Show` to show 50 entries, and write `fluoroquinolone` in the 
   search box. 
4. As you see, gyrA is among the search results.

Thus, be aware of the contents in the data base that you are using. If
what you're looking for isn't in there, you cannot find it!


#### ARIBA and databases

As was evident from looking at the results earlier, ARIBA can work with several 
databases. It can also download the current contents of the database. This means
we can get up-to-date results. 

However, we should be aware that ARIBA requires that the data it searches with 
is in a specific format, this is something that all available tools do. This means
that after downloading data, it will format the contents. Throughout this process,
it might then eliminate some genes that come wrongly formatted from the database.
Thus, this is a second step where we need to check what we're searching with. 

### Running AMR ARIBA

1. Go read the ARIBA manual and have a look at the presentation made by your
   fellow students. Try to figure out the order the following commands have
   to be run in:
   * ariba getref card carddb
   * ariba run --threads 4 cardprepref READ1.FQ READ2.FQ OUTPUTDIRNAME
   * ariba prepareref -f carddb.fa -m carddb.tsv cardprepref
2. Figure out where the reads for your isolate is.
3. Run the first two steps, where we download and prepare a database.
   It might complain that it can't find spades. This is a bug/feature with this
   version, which will be fixed.
   Have a look at the log file for the prepareref step. This will tell you
   about genes and variants that were removed. Also have a look in the 
   `02.cdhit.clusters.tsv` file, this one will let you know about the clusters
   that we're actually including in our search.
3. Run the remaining `run` step.
4. Go into the result directory. Have a look at the `report.tsv` file. Answer
   the following questions. Note. the information found for the `run` and
   `prepareref` commands can be useful for finding answers.
   * Which AMR genes were found?
   * Are these presence/absence AMR genes or mutation AMR genes?
   * If there are any mutation AMR genes, which mutations do you have?
5. Copy your partner's result directory into your own directory. Note, 
   `copy -r` will let you copy a directory with its contents.
   You should now have the two directories side by side.
6. Run `ariba summary comparison TSV_REPORT_F27 TSV_REPORT_S19`.
7. Copy the three resulting files to your computer, and load them
   with the Phandango web browser.

   
## Homework

Run ARIBA mlst and amr detection on the remaining four read sets.

Put up a table that describes your results. This should contain
the MLST type for each genome, and their resistance profile.