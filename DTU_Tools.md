# Technical University of Denmark (DTU) tools
DTU has developed number of tools for whole genome sequencing analysis. And, they are popular and useful. 

## Data location in Saga
### Login to saga
     `ssh yourusername@saga.sigma2.no`
     
     `ls -lh /cluster/projects/nn9305k/tutorial/20201019_DTU_Tools/data/`
     
     `
     -rwxrwxr-x 1 jeevka nn9305k 4.8M Oct 14 10:48 2016-02-522_S70.fasta
     -rwxrwxr-x 1 jeevka nn9305k 4.5M Oct 14 10:48 2016-02-620_S35.fasta
     -rwxrwxr-x 1 jeevka nn9305k 4.8M Oct 14 10:48 2016-17-164_S61.fasta
     -rwxrwxr-x 1 jeevka nn9305k 4.8M Oct 14 10:48 2016-17-292_S51.fasta
     -rwxrwxr-x 1 jeevka nn9305k 4.5M Oct 14 10:48 2016-17-363_S52.fasta
     -rwxrwxr-x 1 jeevka nn9305k 4.6M Oct 14 10:49 2016-17-550_S101.fasta
     `

## Points to remember
These tools can take FastQ reads as the input for their analysis. But they use [Velvet assembler](https://www.ebi.ac.uk/~zerbino/velvet/) to assemble a denova genome and do the downstream analysis using genome.
Since Vetvet assembler is not as good as [SPAdes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/), we suggest that use [SPAdes]((https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/)) or [Shovill](https://github.com/tseemann/shovill) pipeline first to assemble the genome and use the genome as input for these tools.

## ResFinder
[ResFinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) identifies acquired antimicrobial resistance genes in total or partial sequenced isolates of bacteria.

## PointFinder
[PointFinder](https://bitbucket.org/genomicepidemiology/pointfinder/src/master/) service contains one python script PointFinder.py which is the script of the latest version of the PointFinder service. The method detects chromosomal mutations predictive of drug resistance based on WGS data.

## VirulenceFinder
[VirulenceFinder](https://bitbucket.org/genomicepidemiology/virulencefinder/src/master/) service contains one python script virulencefinder.py which is the script of the latest version of the VirulenceFinder service. VirulenceFinder identifies viruelnce genes in total or partial sequenced isolates of bacteria - at the moment only E. coli, Enterococcus, S. aureus and Listeria are available.

## SerotypeFinder
[SerotypeFinder](https://bitbucket.org/genomicepidemiology/serotypefinder/src/master/) service contains one python script serotypefinder.py which is the script of the latest version of the SerotypeFinder service. SerotypeFinder identifies the serotype in total or partial sequenced isolates of E. coli.

## PlasmidFinder
[PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder/src/master/) service contains one python script plasmidfinder.py which is the script of the latest version of the PlasmidFinder service. The service identifies plasmids in total or partial sequenced isolates of bacteria.

## List of all the tools available in DTU

https://bitbucket.org/genomicepidemiology/workspace/repositories
