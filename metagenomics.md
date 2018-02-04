# Introduction to metagenomics.

Metagenomics is a technique to study microbial diversity using environmental DNA samples.
The diversity can either be analyzed using gene-centric / amplicon methods (e.g. PCR fragements of the gene of interest) or via shotgun sequencing. During this course we will start with the bioinformatic analysis of an amplicon dataset.
This will gives us a chance to get used to data clean up methods and get introduced to diversity analysis in R-studio with the phyloseq package.
We will follow the [Mothur MiSeq SOP](https://www.mothur.org/wiki/MiSeq_SOP) for the analysis of illumina amplicon data.

## The course program

* Reducing sequencing and PCR errors
* Generating OTUs and analysis of diversity with Mothur.
* R-studio / Phyloseq tutorial for diversity analysis.
* Introduction to shotgun metagenomics

## Preperation step for the first lecture.
1. Create directory on the desktop in the Linux Virtual box called: `amplicons`
2. Download and save the following datasets in the folder `amplicons`:
	* [Example data from Schloss lab](https://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip)
	* [SILVA-based bacterial reference alignment](https://www.mothur.org/w/images/9/98/Silva.bacteria.zip)
	* [mothur-formated version of the RDP training set (v.9)](https://www.mothur.org/w/images/5/59/Trainset9_032012.pds.zip)

3. extract zipped files using the commandline: `unzipe FILENAME`. This creates two folders: 1) Miseq_SOP and 2) silva.bacteria. and two lose files.
4. move the directory silva.bacteria to the Miseq_SOP folder: mv silva.bacteria Miseq_SOP
5. move the files Trainset* to the Miseq_SOP folder: mv Trainset* Miseq_SOP
