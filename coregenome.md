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

* [Chewbbaca paper](https://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166)
* [Chewbbaca github site](https://github.com/B-UMMI/chewBBACA)
* [Chewbbaca wiki](https://github.com/B-UMMI/chewBBACA/wiki)
* [Chewbbaca tutorial](https://github.com/B-UMMI/chewBBACA_tutorial)
* [Ready schemas](https://zenodo.org/communities/innuendo/?page=1&size=20)

## Practicals

### Goal for the exercise

The main goal when using this program is to start out with a set
of loci (aka the schema) and a set of genomes, and identify
a. if the genomes contain a gene that fits a locus in the schema, and
b. if so, which allele that that genome has for that locus.

We can do this through one of two ways: we can

1. define a schema ourselves
2. use a predefined schema

We will run through the first option first, and through that
show how we can do the second too.

For chewBBACA, a schema in this context is one of two things:

a. a directory with fasta files in it, one for each locus  
b. a text file containing the full path to the fasta files   

The locus fasta file contains fasta sequences for all alleles
found for that locus. In the beginning, after creating or downloading
a schema, the file only contains one sequence. However, when you use
a set of schema files with chewbbaca, any new alleles found for a
locus will be added to the locus file. That means that the schema
will be updated when people run the program.

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
is listed above. We will do an adapted run-through of that.

A subset of the data used for that tutorial can be found here:

```
/projects/nn9305k/bioinf_course/compgenomics/coregenomes/subset_genomes
```

Make a directory in your project home area, and create a directory called
'genomes' in that directory. Copy the genomes in subset_directory into
your 'genomes' directory.

Use this directory when running through this tutorial. For each step,
read throug the information found on the wiki page for that step to
figure out what it does.

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
4. Which genome had the most exact matches?
5. How many paralogs were found?
6. Have a look at the 'GCA-000636115-protein1.fasta' file again. Has it changed?

Go into the results_cg directory and into the results directory. Figure out
what the contents of each of the files there are.

#### ExtractCgMLST step

```
 chewBBACA.py ExtractCgMLST -i results_alleles.tsv -r RepeatedLoci.txt -o cgMLST_completegenomes -p 0.95
```

1. What does the options shown above do?
2. How many loci were deleted, and how many remained?
3. What are the files in the new directory that was created?
4. Have a look at the cgMLSTschema.txt file. This is the final schema file.

NOTE: to get a really good schema, we should have run through the evaluation step
as described in the wiki and the tutorial in addition to what we have done now.
This would have let us identify genomes that we should not include in the schema.
We should also have rerun the AlleleCall/ExtractCgMLST steps until we have
no paralogs. We are skipping this in this tutorial for the sake of time.

#### AlleleCall again

We now have a schema that is without paralogs, and which contains loci that are
present in 95% of the genomes we based the schema on. We will now use it for
allele calling.

##### Add path to allele file

We have a file that contains the fasta file names to the schema. To use it, we will
have to add the path name to the location where they are to the file.

1. Copy your schema file to the directory which also includes schema_seed.
1. Find out where your locus fasta files are. Figure out the full path.
2. Test that you can add the right path to the file by using the following
command:

```
head cgMLSTschema.txt |sed -e "s|^|fullpathhere/|g"
```
Test that you have the right path by cat-ing one of the files that show
up.

3. Save the results to a new file by doing the following:

```
cat cgMLSTschema.txt |sed -e "s|^|fullpathhere/|g" > fullpath_cgMLSTschema.txt
```

##### Do allele call with the new schema

We can now use the same allele call command as above, except with the
'fullpath_cgMLSTschema.txt' file for the '-g' option. Remember to use a
new output directory file name!


### Use a pre-prepared schema

As is clear from what is stated above, a schema can be a text file
containing a list of fasta filenames. When we download a schema from
the website shown above, we get a set of fasta files. These are for the
whole genome MLST schema from Enterobase. However, they also include a
text file for the cgMLST schema. Thus, to use a Enterobase schema,
what we need to do is to take that file, append the full directory
name to the beginning of each line in the file (as demonstrated above)
and use that for the allele calling as shown above.

Please note, using a schema from Enterobase requires quite a bit of
run time, and it might seem like nothing is happening at first.
However, be patient and it will work.
