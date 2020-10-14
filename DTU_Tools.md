# Technical University of Denmark (DTU) tools
DTU has developed number of tools for whole genome sequencing analysis. And, they are popular and useful. 

## Data location in Saga
### Login to saga
     `ssh yourusername@saga.sigma2.no`
### `ls -lh /cluster/projects/nn9305k/`

## Points to remember
These tools can take FastQ reads as the input for their analysis. But they use [Velvet assembler](https://www.ebi.ac.uk/~zerbino/velvet/) to assemble a denova genome and do the downstream analysis using genome.
Since Vetvet assembler is not as good as [SPAdes] (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/), we suggest that use SPAdes or [Shovill](https://github.com/tseemann/shovill) pipeline first to assemble the genome and use the genome as input for these tools.


## ResFinder

## PointFinder

## VirulenceFinder

## SerotypeFinder

## PlasmidFinder


## List of all the tools available in DTU

https://bitbucket.org/genomicepidemiology/workspace/repositories
