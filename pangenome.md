# Pangenomes

There are two main ways of clustering genomes on content.

* SNP based comparisons. In this case we compare the entire genomes,
  find shared genomic regions, and examine these for small mutations.
  This is the analyses we will focus on next week.
* Gene based comparisons. In this case we are looking at the gene
  content, and finding genes that are shared between genomes. To do
  this, we have to have a method for evaluating two sequences and
  saying that these are "the same gene". This is what we will be
  exploring today.
  
In this session, we will be using a program called ROARY to do the
pangenome analysis. ROARY takes the gff files produced by PROKKA as 
input files. The main result from ROARY is a file showing the presence
and absence of genes within the set of genomes that we have given it.
Based on this file, we can deduce several things: 
* How many genes are shared between the genomes
* How many accessory genes are there
* How do the genomes cluster, according to their shared gene content

Note, pangenomes may seem similar to MLST, however they are not. 
With pangenomes, we are effectively taking sequences, and using certain rules
we are saying that these are "the same" (or not). Thus, we are effectively
finding out whether a genome has, for instance, an 'alkB' gene (or not). 
With MLST, we are looking at a gene that _should_ be shared between all
our isolates (as for instance the _E. coli_ MLST 'fumC' gene should be), 
and we are looking at what variant of that gene that we have. Thus, the 
main basis for MLST and its extensions core genome and whole genome MLST, 
are looking at what variants of core genes that we have, while with pangenomes 
we are looking at which genes that are shared between all (core genes), 
and which that are not (accessory genomes).
    
## Links to resources

* [Student presentation](https://docs.google.com/presentation/d/1dvyPqD7nkw-h4X58GNYHTsvdeyvY14ETmCAT2Wcxmyw/edit?usp=sharing)
* [Roary paper](https://academic.oup.com/bioinformatics/article/31/22/3691/240757)
* [Roary supplementary material](Roary_supplementary_material.pdf)
* [Roary website and manual](https://sanger-pathogens.github.io/Roary/)
* [Roary tutorial](https://github.com/microgenomics/tutorials/blob/master/pangenome.md)

## Today's practical

First, read through the Roary website, and the Roary tutorial linked to above.

Question: can you find the criteria that are employed when we have two sequences,
and we are trying to figure out whether these are "the same" or not?

### Running ROARY 

1. Log into abel, and go to your compgenome area
2. Create a directory called 'pangenome', next to the assembly and annotation areas
3. Copy in the gff files of the six _E. coli_ genomes we have been working with.
4. Start a screen, and request a qlogin for 3 hours, and 4 cpus
5. We need to activate roary, and we do that by running 'module load roary'
6. Check out which options are available to roary. Pay special attention 
   to the options '-e' and '-n', and also '-i' and '-cd'
7. Do a first roary run using this command: 
   'roary -p 4 –f basic \*.gff'
   You will see some messages about Citation of Parallell, this is ok. You 
   might also see something about an uninitialized value, this is also ok. 

### Exploring output files
1. Look at the 'summary\_statistics.txt' file. How many genes are shared,
   and how many do we have in total?
2. Copy the 'gene_presence_absence.csv', 'gene_presence_absence.Rtab'
   and the 'accessory_binary_genes.fa.newick' files onto your computer 
   using the 'scp' command.
3. Open the Rtab file using Libreoffice Calc. This file forms the basis
   for what Roary really is. 
   * Select the headings, and go to 'Data > Filter > Autofilter'. 
   * Then sort all in the second column by 'Ascending'.
   * Select the six genome columns, and go to 'Format > Conditional
     Formatting > Condition'. Set the options so that all cells with
     a content of 0 have a red background.
   * Scroll down, and find the first row that is marked with red.
     How does the row number correspond to what you found in 
     the statistics file? Do you understand what the 0s and 1s
     mean?
4. Open the csv file with the same program. This file is in many ways 
   a summary of what's in the Rtab file. 
5. Use the 'Dendroscope' program that you can find on the Bioinformatics
   menu to look at the Newick file. Find an unrooted way of looking at
   the tree. Can you see which isolates share the most genes?

### Visualizing the results

We will now use some python code to create some graphs to illustrate
the results.

1. Make sure you are in the 'basic' directory.
2. Run 'roary\_plots.py -h' and look at the options. Can you figure out 
   which two input files we need, and in which order?
3. Run 'roary\_plots.py' with the required two input files and the
   labels option. It will give some warnings regarding fonts, this is ok.
4. Scp the three resulting png files to your virtual machine, and look
   at them. 
   * The frequency plot will tell you how many genes are in how many
     genomes
   * The pie chart will show you how many genes are in the different
     categories, that is, core, soft-core, cloud, and shell.
   * The matrix file will show you a tree, based on the clustering of
     the presence/absence files. In addition, it will show you the
     "sideways" view of the Rtab files, with 1s shown with blue, and
     0s shown as white.

### Lowering the threshold on "the same"

The default similarity level for ROARY is a minimum of 95%. We will now 
explore what happens if we allow genes to be considered "the same" 
if we allow a lower similarity threshold. 

We will divide into groups, and each group will pick a number between
70 and 90, and use this instead of 'NUMBER' in your command below.

1. Re-run the roary command above, but do the following changes:
   * include '-i NUMBER'
   * change the output directory name to 'lowsim'
2. Go into the 'lowsim' directory and rerun the visualization step.
3. Copy the resulting files, including the newick file, to your laptop,
   while being careful to not overwrite your existing results! 
4. Look at your results, and evaluate with your neighbors. 

### Creating a core genome alignment

So far, all of our results have been focused on what is shared, and
what is not. ROARY works by comparing all of the genes within each 
genome to all the others, this means that it has access to the 
sequences of all the genes. ROARY can also output a core genome
alignment.

1. Read the supplementary document, and try to figure out how
   ROARY creates its alignment.
2. Run roary with the same options as in the first run, i.e. 
   with a similarity threshold of 95, and add the options
   '-e' and '-n', and name the new output directory 'with\_aln'.
   If you forget the '-n' option, you will still get a core
   alignment, but it will be much slower, and potentially 
   more accurate. 
3. Go into the new directory, and see which new files that
   have been produced.
4. Use either 'Jalview', or if that one is slow'



