# Technical University of Denmark (DTU) tools
DTU has developed a number of tools for detecting various genes and mutations related to resistance and virulence, determining serotypes, and also to detect plasmid incompatibility types. The tools are available both online and in the terminal. The online versions of the various tools can be found here: https://cge.cbs.dtu.dk/services/. A list of the command line tools can be found here: 
https://bitbucket.org/genomicepidemiology/workspace/repositories.

The DTU tools can take both reads and assemblies as input, depending on the user. When assemblies are supplied, as exemplified below, the program will use a BLAST method to identify the genes in the databases. However, when supplying reads, a mapping method will be used. In the older versions of these tools, the reads were assembled before the analysis. This process has been exchanged with the mapping procedure instead. In each tool, the genomes are compared against the respective databases (i.e. ResFinder compares the input genome to ResFinder_DB database).

## Data location in Saga
### Login to saga

```
ssh yourusername@saga.sigma2.no
     
ls -lh /cluster/projects/nn9305k/tutorial/20201019_DTU_Tools/data/
     
-rwxrwxr-x 1 jeevka     nn9305k 4.8M Oct 14 10:48 2016-02-522_S70.fasta
-rwxrwxr-x 1 jeevka     nn9305k 4.5M Oct 14 10:48 2016-02-620_S35.fasta
-rwxrwxr-x 1 jeevka     nn9305k 4.8M Oct 14 10:48 2016-17-164_S61.fasta
-rwxrwxr-x 1 jeevka     nn9305k 4.8M Oct 14 10:48 2016-17-292_S51.fasta
-rwxrwxr-x 1 jeevka     nn9305k 4.5M Oct 14 10:48 2016-17-363_S52.fasta
-rwxrwxr-x 1 jeevka     nn9305k 4.6M Oct 14 10:49 2016-17-550_S101.fasta
-rwxrwxrwx 1 hkaspersen nn9305k 5.4M Oct 15 11:22 2018-01-1095_assembly.fasta
```  

## Preparations for the tutorial
This is a temporary directory for this tutorial.

```
cd $USERWORK
mkdir dtu_tools
cd dtu_tools
```

Copy the data to your current directory 

```
rsync -rauWP /cluster/projects/nn9305k/tutorial/20201019_DTU_Tools/data/2018-01-1095_assembly.fasta .
ls -lh
```

screen command gives a possibility to logout and come back to same screen where we ran all the command before
```
screen
```

set up a interactive job, using the development mode (--qos=devel)
```
srun --account=nn9305k --qos=devel --mem-per-cpu=4800M --cpus-per-task=4 --time=0:30:00 --pty bash -i
```

**Note: Conda environment "cge_addons" contains all the dependencies for DTU tools. So, dont need to activate any other conda environment.**

## General usage
Most of the tools follow the same command pattern, with slight variations. In general, the command looks like this:

```
python path/to/script -i path/to/input-file -p path/to/database -m_p path/to/method -o output_folder -x
```
- -m_p (or in some cases, mp) describes which method is used in the analysis (mapping/BLAST) depending on your input. Some tools require you to specify which method you want to use with -m, then supply the method path with -mp.

- The -x commands is for generating extended output, such as a tab delimited file (easier for downstream analysis of results)

The use of ResFinder and PlasmidFinder is presented below.

## ResFinder
[ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) identifies acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria.

```
conda activate cge_addons

mkdir resfinder_output

python /cluster/projects/nn9305k/src/resfinder/resfinder.py -i 2018-01-1095_assembly.fasta -p /cluster/projects/nn9305k/src/resfinder_db/ -m_p  /cluster/software/BLAST+/2.10.1-gompi-2020a/bin/blastn -o resfinder_output/ -x

conda deactivate
```

## VirulenceFinder
[VirulenceFinder](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/) VirulenceFinder identifies virulence genes in total or partial sequenced isolates of bacteria - at the moment only E. coli, Enterococcus, S. aureus and Listeria are available.

```
conda activate cge_addons

mkdir virulencefinder_output

python /cluster/projects/nn9305k/src/virulencefinder/virulencefinder.py -i 2018-01-1095_assembly.fasta -p /cluster/projects/nn9305k/src/virulencefinder_db/ -o virulencefinder_output/ -mp /cluster/software/BLAST+/2.10.1-gompi-2020a/bin/blastn -x

conda deactivate
```



## Database management

The databases used in the tools mentioned above follow the same general structure, with some variations. For simplicity in this tutorial, we will look into detail in the virulencefinder database. 
The database contain fasta files (*.fsa) separated on species and/or gene type. Furthermore, it contains files that hold a short description of the gene (notes.txt), and a file that describes each fasta file in the database (config). It also holds a file called INSTALL.py, which is used when installing the database.
For this tutorial, we only have to think about the config and notes files. 

In this example, the database does not contain known virulence genes from Klebsiella pneumoniae. We will add a gene in a separate fasta file, give it a description in the notes file, and add the fasta file to the list in the config file.

First, copy the whole database into your $USERWORK directory above (it is not big):

```
cp /cluster/projects/nn9305k/src/virulencefinder_db $USERWORK/dtu_tools
ls $USERWORK/dtu_tools
cd $USERWORK/dtu_tools/virulencefinder_db
```

Now we can make alterations to the database (DO NOT alter the shared database in the src/ folder!)

The gene we want to add is listed below:
```
>iucB_1
ATGTCTAAGGCAAACATCGTTCACAGCGGATATGGACTGCGCTGTGAAAAACTCGACAAG
CCTCTGAATCTTAGCTGGGGGCTGGACAATAGCGCGGTGCTGCACTGGCCGGGGCCGTTG
CCGACAGGGTGGCTGCGCGACGCGCTGGAGCAGATATTTATCGCCGCACCGCAACTTTCA
GCGGTGGTTCTCCCTTGGGCCGAATGGCGTGAGGAGCCACAGGCGCTGACGCTTTTCGGG
CAGGTAAAAAGCGACATCATCCATCGCACCGCCTTCTGGCAGTTACCGCTATGGCTGAGT
TCTCCGGCAAACCGGGCCTCCGGCGAAATGGTTTTTGATGCAGAGCGTGAGATTTATTTC
CCGCTGCGCCCACCCCGTCCGCAGGGCGAGGTTTATCGCCGTTACGATCCACGCGTTCGC
AAGACGCTGAGTTTCCGCGTTGCCGATCCCGTTCTTGATGCAGAACGTTTCACCCGCTGG
ATGAACGATCCGCGCGTTGAGTATTTCTGGGAGCAAAGCGGCTCGCTGGAGGTACAGACC
GCCTATCTGGAACGCCAGCTCACCGGTAAACATGCGTTCCCGCTGATCGGCTGCTTCGAC
GATAGGCCGTTTAGCTATTTCGAAATCTACTGGGCGGCGGAAGACCGCATTGGCCGCCAC
TATTCATGGCAACCCTTTGACCGTGGCCTGCATCTGCTGGTTGGTGAACAGCAATGGCGC
GGAGCCCACTACGTACAAAGCTGGCTGCGAGGGCTGACACATTACCTGTTGCTGGATGAG
CCCCGCACGCAGCGCACCGTACTGGAGCCACGCGCCGATAACCAGCGCCTGTTCCGCCAC
CTTGAGCCTGCGGGATACCGGACAATTAAAGAGTTCGACTTCCCACACAAGCGCTCGCGC
ATGGTGATGGCGGATCGCCATCACTTCTTCACGGAGGTCGGTCTGTAA
```

Copy the fasta entry into a new file:

```
nano k.pneumoniae.fsa
```


Save the file. Then, open the config file with nano and add the following line:
(Note that the first column holds the file name without the .fsa ending, followed by a tab separator)
```
k.pneumoniae   Virulence genes for K. pneumoniae
```

Save the file, and open the notes.txt file with nano and add the following line:
```
iucB_1:aerobactin locus gene B type 1
```
Save the file. You are now ready to run VirulenceFinder with your new addition to the database!

```
conda activate cge_addons

mkdir virulencefinder__newdb_output

python /cluster/projects/nn9305k/src/virulencefinder/virulencefinder.py -i 2018-01-1095_assembly.fasta -p $USERWORK/dtu_tools/virulencefinder_db/ -o virulencefinder__newdb_output -mp /cluster/software/BLAST+/2.10.1-gompi-2020a/bin/blastn -x

conda deactivate
```

Take a look at the output files, and compare the results to the run without the added gene (above).
