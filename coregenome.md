# Whole and core genomes with chewbbaca

## Definitions

* allele - different variants of the same genes, i.e. might vary in terms of
  snps and indels
* MLST - multi locus sequence typing. A certain set of genes is chosen,
  normally 7. Different alleles of these genes are given a number. Note, this
  number is only a label, it does not reflect number of changes in relation to
  any other allele. The set of labels constitute a certain Sequence Type (ST).
* MLST schema - the set of genes that is chosen as the basis for MLST.
* whole genome MLST - a set of genomes is examined to find the set of genes that
  are shared by all the genomes. This set of genes forms the schema for a whole
  genome MLST analysis
* core genome MLST - a preselected set of genes that are deemed to be present
  in all genomes from that subtype or species. These form the core genome MLST
  schema. Note, this means that for some genomes, some of these genes might
  be missing.

## Links to resources

* [Chewbbaca paper]()
* [Chewbbaca github site](https://github.com/B-UMMI/chewBBACA)
* [Chewbbaca tutorial](https://github.com/B-UMMI/chewBBACA_tutorial)
* [Ready schemas](https://zenodo.org/communities/innuendo/?page=1&size=20)

## Practicals

### How to use chewBBACA

#### Request resources

When using this program, we have to remember to request compute resources,
in this case via a qlogin.

```
qlogin --account=nn9305k --time=02:00:00 --ntasks=6 --mem-per-cpu=4G
```

#### Activate the program

To use chewBBACA, you need to activate the conda environment that contains
the program. To do that, do

```
source activate chewbbaca
```
To use chewBBACA, the general syntax is

```
chewBBACA.py <subcommand> <options>
```

The available subcommands can be listed by doing

```
chewBBACA.py --help
```

### Create wgMLST and cgMLST schema

The chewBBACA developers have come up with a very good tutorial, which
is listed above.

A subset of the data used for that tutorial can be found here:

```
/projects/nn9305k/bioinf_course/compgenomics/coregenomes/subset_genomes
```

Make a directory in your project home area, and create a directory called
'genomes' in that directory. Copy the genomes in subset_directory into
your 'genomes' directory.

Use this directory when running through this tutorial.

#### CreateSchema step

````
chewBBACA.py CreateSchema -i genomes/ --cpu 6 -o schema_seed --ptf Streptococcus_agalactiae.trn
````
1. What are the options given above?
2. Which new files/directories were produced?
3. What is in the proteinID_Genome.tsv file?
4. What is inside the directory that was just made?
5. Can you figure out what the contents of the 'GCA-000636115-protein1.fasta' is?

#### AlleleCall step

```
chewBBACA.py AlleleCall -i genomes/ -g schema_seed/ -o results_cg --cpu 6 --ptf Streptococcus_agalactiae.trn
```
1. What do the options shown above mean?
2. Which bsr value was used?
3. How many exact matches were found?
4. Which genome had the most exact matches_
5. How many paralogs were found?

Go into the results_cg directory and into the results directory. Figure out
what the contents of each of the files there are.

#### ExtractCgMLST step

```
 chewBBACA.py ExtractCgMLST -i results_alleles.tsv -r RepeatedLoci.txt -o cgMLST_completegenomes -p 0.95
```

1. What does the options shown above do?
1. How many loci were deleted, and how many remained?
2. What are the files in the new directory that was created?
3. Use the awk command shown above on the cgMLST.tsv file. Can you find some
   columns that are different from the various genomes?









### Use a pre-prepared schema
