# Creating User Defined Databases

## Introduction

Often the new/latest species/data/results (AMR/Virulence/MLST…etc,.) not included in the databases (DBs) which are used by tools. But it is possible to use the latest data with some of the tools. This documentation will walk you through "How to create your own database" for three of the most commonly used tools: 
1. [ARIBA](https://github.com/sanger-pathogens/ariba): ReadMapping Approach
2. [ResFinder and PointFinder](https://bitbucket.org/genomicepidemiology/resfinder.git/src): Local alignment Approach
3. [MLST](https://github.com/tseemann/mlst): Local alignment 

[Link to Introduction Slides](https://vetinst.sharepoint.com/:p:/r/sites/Bioinfadmins/Delte%20dokumenter/General/Bioinfo_Training/UD_DBs_Jeevan_300522.pptx?d=w43117335fecf47edb54fb34d3db21e22&csf=1&web=1&e=CJLpPE)

## ARIBA:
Tool Input: Assembly (.fasta) or reads (.fastq). 
For creating DB: fasta sequence with gene name, which wi

### Working location for hands-on in SAGA

Login to SAGA
```
ssh username@saga.sigma2.no
```

Activate ARIBA conda environment  
```
conda activate ARIBA
```

Go into test directory.  
```
cd /cluster/projects/nn9305k/tutorial/300522_User_Defined_DBs/

ls -lh 

total 1.0K
drwxrwsr-x 2 jeevka nn9305k 1 May 27 11:12 Additional_Virulence_Genes
drwxrwsr-x 2 jeevka nn9305k 2 May 27 11:11 Test_Fastq_Files
```
Additional_Virulence_Genes: contains the fasta sequences of additional virulence genes.  
Test_Fastq_Files: contains a test fastq paired end sequences


Creating data base from given fasta file 
```
ariba prepareref --all_coding yes -f Add_Virulence_Genes/Additional_Vir_genes.fasta ariba_virulence_genes_DB
```
Add_Virulence_Genes/Additional_Vir_genes.fasta : input files. 
ariba_virulence_genes_DB : output folder which contains the db


Mapping the reads to newly created db for new virulence genes
```
ariba run ariba_virulence_genes_DB/ Test_Fastq_Files/Vibrio_R1.fastq Test_Fastq_Files/Vibrio_R2.fastq Virulence_Output
```

**Virulence_Output/** contains the output files.

"report.tsv" - file contains the results

[Explanation for Output](https://github.com/sanger-pathogens/ariba/wiki/Task:-run)


Deactivate ARIBA conda environment 
```
conda deactivate 
```

## ResFinder/PointFinder
**Input**: Assembly (.fasta) or reads (.fastq)

**To ADD new genes to DB**: 
1. <resistance_phenotye_class>.fsa file: multifasta file of AMR genes. f. ex. aminoglycoside.fsa
```
head aminoglycoside.fsa

>aac(6')-Ib_2_M23634
ATGAGTATTCAACATTTCCAAAGAAAGTTAGGCATCACAAAGTACAGCATCGTGACCAAC
AGCAACGATTCCGTCACACTGCGCCTCATGACTGAGCATGACCTTGCGATGCTCTATGAG
TGGCTAAATCGATCTCATATCGTCGAGTGGTGGGGCGGAGAAGAAGCACGCCCGACACTT
GCTGACGTACAGGAACAGTACTTGCCAAGCGTTTTAGCGCAAGAGTCCGTCACTCCATAC
ATTGCAATGCTGAATGGAGAGCCGATTGGGTATGCCCAGTCGTACGTTGCTCTTGGAAGC
GGGGACGGATGGTGGGAAGAAGAAACCGATCCAGGAGTACGCGGAATAGACCAGTTACTG
GCGAATGCATCACAACTGGGCAAAGGCTTGGGAACCAAGCTGGTTCGAGCTCTGGTTGAG
TTGCTGTTCAATGATCCCGAGGTCACCAAGATCCAAACGGACCCGTCGCCGAGCAACTTG
CGAGCGATCCGATGCTACGAGAAAGCGGGGTTTGAGAGGCAAGGTACCGTAACCACCCCA
GATGGTCCAGCCGTGTACATGGTTCAAACACGCCAGGCATTCGAGCGAACACGCAGGTTT
GCCTAA
>aac(6')-Ib11_1_AY136758
ATGAAAAACACAATACATATCAACAGCAACGATTCCGTCACACTGCGCCTCATGACTGAG
CATGACCTTGCGATGCTCTATGAGTGGCTAAATCGATCTCATATCGTCGAGTGGTGGGGC
GGAGAAGAAGCACGCCCGACACTTGCTGACGTACAGGAACAGTACTTGCCAAGCGTTTTA
GCGCAAGAGTCCGTCACTCCATACATTGCAATGCTGAATGGAGAGCCGATTGGGTATGCC
CAGTCGTACGTTGCTCTTGGAAGCGGGGACGGATGGTGGGAAGAAGAAACCGATCCAGGA
GTACGCGGAATAGACCTGTCACTGGCGAATGCATCACAACTGGGCAAAGGCTTGGGAACC
AAGCTGGTTCGAGCTCTGGTTGAGTTGCTGTTCAATGATCCCGAGGTCACCAAGATCCAA
ACGGACCCGTCGCCGAGCAACTTGCGAGCGATCCGATGCTACGAGAAAGCGGGGTTTGAG
AGGCAAGGTACCGTAACCACCCCAGATGGTCCAGCCGTGTACATGGTTCAAACACGCCAG
GCATTCGAGCGAACACGCAGTGATGCCTAA
```

2. "phenotypes.txt" file : 
```
head phenotypes.txt

Gene_accession no.	Class	Phenotype	PMID	Mechanism of resistance	Notes	Required_gene
ant(2'')-Ia_1_X04555	Aminoglycoside	Gentamicin, Tobramycin	3024112	Enzymatic modification	Alternative name aadB
ant(2'')-Ia_10_HM367617	Aminoglycoside	Gentamicin, Tobramycin	21873033	Enzymatic modification
ant(2'')-Ia_11_HM367620	Aminoglycoside	Gentamicin, Tobramycin	21873033	Enzymatic modification
ant(2'')-Ia_12_HQ880250	Aminoglycoside	Gentamicin, Tobramycin	unpublished	Enzymatic modification
ant(2'')-Ia_13_DQ176450	Aminoglycoside	Gentamicin, Tobramycin	16304199	Enzymatic modification
ant(2'')-Ia_14_DQ266447	Aminoglycoside	Gentamicin, Tobramycin	unpublished	Enzymatic modification
ant(2'')-Ia_15_EF205594	Aminoglycoside	Gentamicin, Tobramycin	unpublished	Enzymatic modification
ant(2'')-Ia_16_HQ386848	Aminoglycoside	Gentamicin, Tobramycin	unpublished	Enzymatic modification
ant(2'')-Ia_17_JTTZ01000034	Aminoglycoside	Gentamicin, Tobramycin	unpublished	Enzymatic modification

```

3. "phenotype_panels.txt":

```
head phenotype_panels.txt

:Panel: Campylobacter
Ampicillin

:Panel: Campylobacter jejuni
:Include: Campylobacter
Ciprofloxacin
Erythromycin
Gentamicin
Nalidixic acid
Streptomycin
Tetracycline

:Panel: Campylobacter coli
:Include: Campylobacter
Ciprofloxacin
Erythromycin
Gentamicin
Streptomycin
Tetracycline
``` 

Activate conda environment for ResFinder
```
conda activate resfinder 
```

Go into the db_folder 
```
cd db_resfinder
ls -lh 

total 3.6M
-rw-rw-r-- 1 jeevka nn9305k 197K May 27 13:04 aminoglycoside.fsa
-rw-rw-r-- 1 jeevka nn9305k 2.8K May 27 13:04 antibiotic_classes.txt
-rw-rw-r-- 1 jeevka nn9305k 1.8M May 27 13:04 beta-lactam.fsa
-rwxrwxr-x 1 jeevka nn9305k 2.4K May 27 13:04 CHECK-entries.sh
-rw-rw-r-- 1 jeevka nn9305k  92K May 27 13:04 colistin.fsa
-rw-rw-r-- 1 jeevka nn9305k  912 May 27 13:04 config
-rw-rw-r-- 1 jeevka nn9305k  40K May 27 13:04 history.txt
-rwxrwxr-x 1 jeevka nn9305k 3.8K May 27 13:04 INSTALL.py
-rw-rw-r-- 1 jeevka nn9305k 173K May 27 13:04 macrolide.fsa
-rw-rw-r-- 1 jeevka nn9305k 6.9K May 27 13:04 nitroimidazole.fsa
-rw-rw-r-- 1 jeevka nn9305k  89K May 27 13:04 notes.txt
-rw-rw-r-- 1 jeevka nn9305k  45K May 27 13:04 oxazolidinone.fsa
-rw-rw-r-- 1 jeevka nn9305k  44K May 27 13:04 phenicol.fsa
-rw-rw-r-- 1 jeevka nn9305k 2.6K May 27 13:04 phenotype_panels.txt
-rw-rw-r-- 1 jeevka nn9305k 504K May 27 13:04 phenotypes.txt
-rw-rw-r-- 1 jeevka nn9305k 9.3K May 27 13:04 pseudomonicacid.fsa
-rw-rw-r-- 1 jeevka nn9305k  92K May 27 13:04 quinolone.fsa
-rw-rw-r-- 1 jeevka nn9305k 5.4K May 27 13:04 README.md
```

Make/Build DB
```
python3 INSTALL.py </path/to/kma_index> non_interactive
```


## MLST 
**Input File**:  reads (.fastq)
**To ADD new schemes**: Fasta file with allele sequences

Create a folder for new scheme
```
mkdir sareus
```

create/copy the necessary files insie the "sareus" folder
```
ls 
saureus.txt
arcC.tfa
aroE.tfa
glpF.tfa
gmk.tfa
pta.tfa
tpi.tfa
yqiL.tfa
```

```
head -n 5 saureus.txt
ST      arcC    aroE    glpF    gmk     pta     tpi     yqiL    clonal_complex
1       1       1       1       1       1       1       1
2       2       2       2       2       2       2       26
3       1       1       1       9       1       1       12
4       10      10      8       6       10      3       2
```


Alle sequences
```
head -n 20 arcC.tfa
>arcC_1
TTATTAATCCAACAAGCTAAATCGAACAGTGACACAACGCCGGCAATGCCATTGGATACT
TGTGGTGCAATGTCACAGGGTATGATAGGCTATTGGTTGGAAACTGAAATCAATCGCATT
TTAACTGAAATGAATAGTGATAGAACTGTAGGCACAATCGTTACACGTGTGGAAGTAGAT
AAAGATGATCCACGATTCAATAACCCAACCAAACCAATTGGTCCTTTTTATACGAAAGAA
GAAGTTGAAGAATTACAAAAAGAACAGCCAGACTCAGTCTTTAAAGAAGATGCAGGACGT
GGTTATAGAAAAGTAGTTGCGTCACCACTACCTCAATCTATACTAGAACACCAGTTAATT
CGAACTTTAGCAGACGGTAAAAATATTGTCATTGCATGCGGTGGTGGCGGTATTCCAGTT
ATAAAAAAAGAAAATACCTATGAAGGTGTTGAAGCG
>arcC_2
TTATTAATCCAACAAGCTAAATCGAACAGTGACACAACGCCGGCAATGCCATTGGATACT
TGTGGTGCAATGTCACAAGGTATGATAGGCTATTGGTTGGAAACTGAAATCAATCGCATT
TTAACTGAAATGAATAGTGATAGAACTGTAGGCACAATCGTAACACGTGTGGAAGTAGAT
AAAGATGATCCACGATTTGATAACCCAACTAAACCAATTGGTCCTTTTTATACGAAAGAA
GAAGTTGAAGAATTACAAAAAGAACAGCCAGGCTCAGTCTTTAAAGAAGATGCAGGACGT
GGTTATAGAAAAGTAGTTGCGTCACCACTACCTCAATCTATACTAGAACACCAGTTAATT
CGAACTTTAGCAGACGGTAAAAATATTGTCATTGCATGCGGTGGTGGCGGTATTCCAGTT
ATAAAAAAAGAAAATACCTATGAAGGTGTTGAAGCG
```


Activate conda environment for MLST


```
conda activate MLST
```

