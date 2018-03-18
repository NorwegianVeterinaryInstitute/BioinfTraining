# Teaching notes

## 2018-03-19

### MLST databases

Genes shared between MLST schemes:

Escherichia coli#1	7413 profiles	2018-03-15
adk; fumC; gyrB; icd; mdh; purA; recA

Escherichia coli#2	851 profiles	2018-03-15
dinB; icdA; pabB; polB; putP; trpA; trpB; uidA

Maybe one, icd

Basic thought: we are finding what we are looking for, i.e. 
we have to check what we're looking for. 

Camilla can explain why we use #1.

READS:

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


### MLST ariba finding






## 2018-03-12

### Homework

* Transfer the files to your computer, and open them in firefox. Most 
  of the data is untrimmed, and 150 bp or shorter. Only the ESBL is longer.
    
  Re fastqc results: you have to know your data. This time, we went into 
  it blind. We didn't check what we had. We could have had very crappy
  data. As it was, most of it was fine. However, check what you have.
  When we assembled things, we didn't even know the read length, which
  complicated things.
    
    * 151: 14042624, DTU2014 - can use auto
    * 100: F27, F159, F168 - can use auto
    * variable: S19 - need to fix
    
* Open one of the files with `less`. As you see it says HISEQ here.
    If you open the ESBL file, we will see M- which indicates Miseq data
    
* Re SPAdes manual: we need to make sure that we're looking at the manual
  for the right version.
    
  How to figure out what version of spades we have? 
    
  module load spades
  spades.py -v 
  
* QUAST results 


# scratch








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


    
    