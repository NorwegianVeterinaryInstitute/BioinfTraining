# Technical University of Denmark (DTU) tools
DTU has developed a number of tools for whole genome sequencing analysis. And, they are popular and useful. 

## Data location in Saga
### Login to saga

```
ssh yourusername@saga.sigma2.no
     
ls -lh /cluster/projects/nn9305k/tutorial/20201019_DTU_Tools/data/
     
-rwxrwxr-x 1 jeevka nn9305k 4.8M Oct 14 10:48 2016-02-522_S70.fasta
-rwxrwxr-x 1 jeevka nn9305k 4.5M Oct 14 10:48 2016-02-620_S35.fasta
-rwxrwxr-x 1 jeevka nn9305k 4.8M Oct 14 10:48 2016-17-164_S61.fasta
-rwxrwxr-x 1 jeevka nn9305k 4.8M Oct 14 10:48 2016-17-292_S51.fasta
-rwxrwxr-x 1 jeevka nn9305k 4.5M Oct 14 10:48 2016-17-363_S52.fasta
-rwxrwxr-x 1 jeevka nn9305k 4.6M Oct 14 10:49 2016-17-550_S101.fasta
```  

## Techincal preparation for tutorial
This is a temporary directory for this tutorial.

```
cd $USERWORK
mkdir dtu_tools
cd dtu_tools
```

Copy the data to your current directory 

```
rsync -rauWP /cluster/projects/nn9305k/tutorial/20201019_DTU_Tools/data/2016-02-522_S70.fasta .
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

## Points to remember
These tools can take both reads and assemblies as input, depending on the user. When assemblies are supplied, as exemplified below, the program will use a BLAST method to identify the genes in the databases. However, when supplying reads, a mapping method will be used. In the older versions of these tools, the reads were assembled before the analysis. This process has been exchanged with the mapping procedure instead. In each tool, the genomes are compared against the respective databases (f. ex. ResFinder compares the input genome to ResFinder_DB database)

## ResFinder
[ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) identifies acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria.

```
conda activate cge_addons

mkdir resfinder_output

python /cluster/projects/nn9305k/src/resfinder/resfinder.py -i 2016-02-522_S70.fasta -p /cluster/projects/nn9305k/src/resfinder_db/ -m_p  /cluster/software/BLAST+/2.10.1-gompi-2020a/bin/blastn -o resfinder_output/ -x

conda deactivate
```

## PointFinder
[PointFinder](https://bitbucket.org/genomicepidemiology/pointfinder/src/master/) service contains one python script PointFinder.py which is the script of the latest version of the PointFinder service. The method detects chromosomal mutations predictive of drug resistance based on WGS data.


```
conda activate cge_addons

mkdir pointfinder_output

python /cluster/projects/nn9305k/src/pointfinder/PointFinder.py -i 2016-02-522_S70.fasta -p /cluster/projects/nn9305k/src/pointfinder_db/ -o pointfinder_output/ -s escherichia_coli -m blastn -m_p /cluster/software/BLAST+/2.10.1-gompi-2020a/bin/blastn -x

conda deactivate
```

## VirulenceFinder
[VirulenceFinder](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/) service contains one python script virulencefinder.py which is the script of the latest version of the VirulenceFinder service. VirulenceFinder identifies virulence genes in total or partial sequenced isolates of bacteria - at the moment only E. coli, Enterococcus, S. aureus and Listeria are available.

```
conda activate cge_addons

mkdir virulencefinder_output

python /cluster/projects/nn9305k/src/virulencefinder/virulencefinder.py -i 2016-02-522_S70.fasta -p /cluster/projects/nn9305k/src/virulencefinder_db/ -o virulencefinder_output/ -mp /cluster/software/BLAST+/2.10.1-gompi-2020a/bin/blastn -x

conda deactivate
```

## SerotypeFinder
[SerotypeFinder](https://bitbucket.org/genomicepidemiology/serotypefinder/src/master/) service contains one python script serotypefinder.py which is the script of the latest version of the SerotypeFinder service. SerotypeFinder identifies the serotype in total or partial sequenced isolates of E. coli.

```
conda activate cge_addons

mkdir serotypefinder_output

python /cluster/projects/nn9305k/src/serotypefinder/serotypefinder.py -i 2016-02-522_S70.fasta -p /cluster/projects/nn9305k/src/serotypefinder_db/ -o serotypefinder_output/ -mp /cluster/software/BLAST+/2.10.1-gompi-2020a/bin/blastn -x

conda deactivate
``` 

## PlasmidFinder
[PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/) service contains one python script plasmidfinder.py which is the script of the latest version of the PlasmidFinder service. The service identifies plasmids in total or partial sequenced isolates of bacteria.


```
conda activate cge_addons

mkdir plasmidfinder_output

python /cluster/projects/nn9305k/src/plasmidfinder/plasmidfinder.py -i 2016-02-522_S70.fasta -p /cluster/projects/nn9305k/src/plasmidfinder_db/ -mp /cluster/software/BLAST+/2.10.1-gompi-2020a/bin/blastn -x -o plasmidfinder_output

conda deactivate
```

## List of all the tools available in DTU

https://bitbucket.org/genomicepidemiology/workspace/repositories


## Database management

Database-things goes here
