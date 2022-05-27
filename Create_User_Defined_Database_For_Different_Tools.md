# Creating User Defined Databases

## Introduction

Often the new/latest species/data/results (AMR/Virulence/MLST…etc,.) not included in the databases (DBs) which are used by tools. But it is possible to use the latest data with some of the tools. This documentation will walk you through "How to create your own database" for three of the most commonly used tools: 
1. [ARIBA](https://github.com/sanger-pathogens/ariba): ReadMapping Approach
2. [ResFinder and PointFinder](https://bitbucket.org/genomicepidemiology/resfinder.git/src): ReadMapping Approach
3. [MLST](https://github.com/tseemann/mlst): Local alignment 

[Link to Introduction Slides](https://vetinst.sharepoint.com/:p:/r/sites/Bioinfadmins/Delte%20dokumenter/General/Bioinfo_Training/UD_DBs_Jeevan_300522.pptx?d=w43117335fecf47edb54fb34d3db21e22&csf=1&web=1&e=CJLpPE)

## ARIBA as an example. :
Input: Assembly (.fasta) or reads (.fastq)


## Working location for hands-on in SAGA

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

```
ariba run ariba_virulence_genes_DB/ Test_Fastq_Files/Vibrio_R1.fastq Test_Fastq_Files/Vibrio_R2.fastq Virulence_Output
```

