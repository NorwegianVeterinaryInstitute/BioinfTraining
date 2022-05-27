# Introduction (PowerPoint slides)
## Read Mapping
 1. Read Coverage on target
 2. Some of the tools used for read mapping approach

## Local alignment (BLAST)
 1. How much % of the target is covered and similarity 

# Working location for hands-on in SAGA

```

ssh username@saga.sigma2.no

conda activate ARIBA

cd /cluster/shared/vetinst/ARIBA_DB/

ariba prepareref --all_coding yes -f Add_Virulence_Genes/Additional_Vir_genes.fasta ariba_virulence_genes_DB

ariba run ariba_virulence_genes_DB/ Test_Fastq_Files/Vibrio_R1.fastq Test_Fastq_Files/Vibrio_R2.fastq Virulence_Output

```

# ARIBA: ReadMapping Approach
Using Fastq (read) files illumina to search against user defined AMR, Virulence and MLST using ARIBA.

### AMR/Virulence
1. Live demo of creating database from amr gene sequences (fasta)



### MLST
1. Live demo of creating database from new MLST schemes 


# ResFinder and PointFinder: ReadMapping Approach


# MLST tool: Local alignment 
Adding new MLST schemes to the existing schemes

