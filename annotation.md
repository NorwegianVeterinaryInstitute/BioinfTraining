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

### Basic PROKKA use

First, we are going to do a first PROKKA run. We will annotate the SPAdes
assembly that we have been creating in class. You can find this assembly here:

`where spades asm is on abel`.

1. Go read the website above, to learn about the options to PROKKA.
2. Log in to abel, go to the course directory in your project home area.
3. Create a directory named `annotation` next to the `assembly` directory.
4. Start a new `screen`, and ask for a `qlogin` session for 3 hours and 4 CPUs.
5. Enable the two following modules: prokka and hmmer/3.1b2 .
6. Type in `prokka -h` and look at the options, specifically, 
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

1. Have a look at the `txt` file. What kinds of genes have been found, 
   and how many?
2. Have a look at the `gbk` produced by PROKKA using `less`.
3. If you have some favorite genes, see if you can find them. If not,
   see if you can find the `gyrA` protein. To search for a text string
   in `less`, type in `/` followed what you're searching for, and press
   enter.

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
   --outdir,--locustag, --genus, and --usegenus 
4. Once the run is done, open the Genbank file for the two runs, each in a 
   separate terminal. Find your favorite gene (or `gyrA`) and see if you 
   can see any difference in how they are annotated. 








1. What are genes
  * show anatomy of a gene
2. What kinds of genes do we have
  * Protein coding
  * rRNA genes
  * tRNA genes
  * ncRNA
  * mobile elements
  * CDS
3. Show a genbank file, describe how to read it

prokka -h - look at options

Read prokka -docs

Figure out output files

run first without anything except compliant and cpus
prokka scaffolds.fasta --compliant --cpus 4

look at txt file
look at gbk file

[karinlag@abel 2018-03-05-coursetest]$ prokka scaffolds.fasta --compliant --cpus 4 --outdir genus --genus Escherichia --locustag myeco


show difference with/without compliant

show difference with/without genus


prodigal -i PROKKA_03022018/PROKKA_03022018.fna -c -m -g 11 -p single -f sco -q

 -c:  Closed ends.  Do not allow genes to run off edges.
 -m:  Treat runs of n's as masked sequence and do not build genes across them
 -g:  Specify a translation table to use (default 11).
 -p:  Select procedure (single or meta).  Default is single.
 -f:  Select output format (gbk, gff, or sco).  Default is gbk.
 -q:  Run quietly (suppress normal stderr output).

Annotation of proteins

An optional user-provided set of annotated proteins. These are expected to be trustworthy curated datasets and will be used as the primary source of annotation. They are searched using BLAST+ blastp ( Camacho et al. , 2009 ).

All bacterial proteins in UniProt ( Apweiler et al. , 2004 ) that have real protein or transcript evidence and are not a fragment. This is ∼16 000 proteins, and typically covers >50% of the core genes in most genomes. BLAST+ is used for the search.

All proteins from finished bacterial genomes in RefSeq for a specified genus. This captures domain-specific naming, and the databases vary in size and quality, depending on the popularity of the genus. BLAST+ is used for this and is optional.

A series of hidden Markov model profile databases, including Pfam ( Punta et al. , 2012 ) and TIGRFAMs ( Haft et al. , 2013 ). This is performed using hmmscan from the HMMER 3.1 package ( Eddy, 2011 ).

If no matches can be found, label as ‘hypothetical protein’.

  
Have to find a way to make sure that they know that these are preedictions and not truth



read all the files
Prokka is 'dumb'. It searches different databases in the order of (--proteins, --genus, in built SPROT, then the HAMAP HMM etc). If it matches your one exacrtly in proteins it should annotate it as such.
ASSUMING that Prodigal actually finds the gene that is.

@rattei i have had some confusion here, yes, I apologise. The contig name length is only checked when you use --compliant. If you don't want compliant files, you get back what you put in. The GFF3 and FASTA files are all perfectly valid with any contig ID length. It's only GBK which is the problem.

THat's the crux of the issue. If you put in normal contig names, ALL output files are fine. if you put in long ones and don't use --compliant, you get ALL output files being good EXCEPT for GBK. I'm stuck between a rock and a hard place here.

Annotations:
  --kingdom [X]     Annotation mode: Archaea|Bacteria|Viruses (default 'Bacteria')
  --gcode [N]       Genetic code / Translation table (set if --kingdom is set) (default '0')
  --gram [X]        Gram: -/neg +/pos (default '')
  --usegenus        Use genus-specific BLAST databases (needs --genus) (default OFF)
  --proteins [X]    Fasta file of trusted proteins to first annotate from (default '')
  --metagenome      Improve gene predictions for highly fragmented genomes (default OFF)


[karinlag@abel 2018-03-05-coursetest]$ grep -A 10 "complement(3072..3635" PROKKA_03022018/PROKKA_03022018.gbk 
     gene            complement(3072..3635)
                     /gene="gntK_1"
                     /locus_tag="PROKKA_00005"
     CDS             complement(3072..3635)
                     /gene="gntK_1"
                     /locus_tag="PROKKA_00005"
                     /EC_number="2.7.1.12"
                     /inference="ab initio prediction:Prodigal:2.60"
                     /inference="similar to AA sequence:UniProtKB:P46859"
                     /codon_start=1
                     /transl_table=11
                     /product="Thermoresistant gluconokinase"
                     /protein_id="Prokka:PROKKA_00005"
                     /translation="MAGESFILMGVSGSGKTLIGSKVAALLSAKFIDGDDLHPAKNID
[karinlag@abel 2018-03-05-coursetest]$ grep -A 10 "complement(3072..3635" usegenus/myeco_03022018.gbk 
     gene            complement(3072..3635)
                     /gene="idnK"
                     /locus_tag="myeco_00005"
     CDS             complement(3072..3635)
                     /gene="idnK"
                     /locus_tag="myeco_00005"
                     /EC_number="2.7.1.12"
                     /inference="ab initio prediction:Prodigal:2.60"
                     /inference="similar to AA sequence:RefSeq:YP_002405683.1"
                     /codon_start=1
                     /transl_table=11
                     /product="D-gluconate kinase"
                     /protein_id="Prokka:myeco_00005"
                     /translation="MAGESFILMGVSGSGKTLIGSKVAALLSAKFIDGDDLHPAKNID
[karinlag@abel 2018-03-05-coursetest]$ 
