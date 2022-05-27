# Creating User Defined Databases

## Introduction

Often the new/latest species/data/results (AMR/Virulence/MLST…etc,.) not included in the databases (DBs) which are used by tools. But it is possible to use the latest data with some of the tools. This documentation will walk you through "How to create your own database" for three of the most commonly used tools: 
1. [ARIBA](https://github.com/sanger-pathogens/ariba): ReadMapping Approach
2. [ResFinder and PointFinder](https://bitbucket.org/genomicepidemiology/resfinder.git/src): Local alignment Approach
3. [MLST](https://github.com/tseemann/mlst): Local alignment 

[Link to Introduction Slides](https://vetinst.sharepoint.com/:p:/r/sites/Bioinfadmins/Delte%20dokumenter/General/Bioinfo_Training/UD_DBs_Jeevan_300522.pptx?d=w43117335fecf47edb54fb34d3db21e22&csf=1&web=1&e=CJLpPE)

## ARIBA:
Tool Input: Assembly (.fasta) or reads (.fastq)
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


## ResFinder/PointFinder
**Input**: Assembly (.fasta) or reads (.fastq)

**To create DB**: 
1. <resistance_phenotye_class>.fsa file: multifasta file of AMR genes. f. ex. aminoglycoside.fsa
2. phenotypes.txt file : 
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

4.  
5. 

