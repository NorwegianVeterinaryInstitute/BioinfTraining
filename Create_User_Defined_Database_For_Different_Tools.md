# Creating User Defined Databases

## Introduction

Often the new/latest species/data/results (AMR/Virulence/MLST…etc,.) not included in the databases (DBs) which are used by tools. But it is possible to use the latest data with some of the tools. This documentation will walk you through "How to create your own database" for three of the most commonly used tools: 
1. [ARIBA](https://github.com/sanger-pathogens/ariba): ReadMapping Approach
2. [ResFinder and PointFinder](https://bitbucket.org/genomicepidemiology/resfinder.git/src): ReadMapping Approach
3. [MLST](https://github.com/tseemann/mlst): Local alignment 

## Read Mapping
 1. Read Coverage on target
 2. Some of the tools used for read mapping approach

## Local alignment (BLAST)
 1. How much % of the target is covered and similarity 

[Link to Introduction Slides](https://vetinst.sharepoint.com/:p:/r/sites/Bioinfadmins/Delte%20dokumenter/General/Bioinfo_Training/UD_DBs_Jeevan_300522.pptx?d=w43117335fecf47edb54fb34d3db21e22&csf=1&web=1&e=CJLpPE)

## ARIBA as an example. :
Input: Assembly (.fasta) or reads (.fastq)



## Working location for hands-on in SAGA

```

ssh username@saga.sigma2.no

conda activate ARIBA

cd /cluster/shared/vetinst/ARIBA_DB/

ariba prepareref --all_coding yes -f Add_Virulence_Genes/Additional_Vir_genes.fasta ariba_virulence_genes_DB

ariba run ariba_virulence_genes_DB/ Test_Fastq_Files/Vibrio_R1.fastq Test_Fastq_Files/Vibrio_R2.fastq Virulence_Output

```
