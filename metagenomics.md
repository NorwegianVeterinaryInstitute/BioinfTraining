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

## Lecture 2018 – 5 – February

During this day we start with the Mothur SOP for [MiSeq data](https://www.mothur.org/wiki/MiSeq_SOP) and we will follow the tutorial using the Biolinux mothur version in our virtual box environment.

**important note.** 
You can work yourself through this tutorial, but note that the mothur version on your biolinux is not up to data with the latests mothur versions. This means that for a few commands we need to use an different command than is used in the mothur SOP. Below I will give a link to the commands that we used during the course. This also means that when you work through the tutorial you will not find the same numbers for reads kept or removed.

### Preperation step for the first lecture.
1. Create directory on the desktop in the Linux Virtual box called: `amplicons`
2. Download and save the following datasets in the folder `amplicons`:
	* [Example data from Schloss lab](https://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip)
	* [SILVA-based bacterial reference alignment](https://www.mothur.org/w/images/9/98/Silva.bacteria.zip)
	* [mothur-formated version of the RDP training set (v.9)](https://www.mothur.org/w/images/5/59/Trainset9_032012.pds.zip)

3. extract zipped files using the commandline: `unzipe FILENAME`. This creates two folders: 1) Miseq_SOP and 2) silva.bacteria. and two lose files.
4. move the directory silva.bacteria to the `Miseq_SOP` folder: `mv silva.bacteria Miseq_SOP`
5. move the files Trainset* to the `Miseq_SOP` folder: `mv Trainset* Miseq_SOP`

### The tutorial steps

Open file finder window and move to folder amplicons, check the files that you have downloaded.
Since we have a lot of fastq data and we are curious, we want to get an idea of the quality of the sequences. We can use the tool `fastqc` which is installed in the biolinux machine. To run this on all the MiSeq fastq datafiles of this tutorial do the following:

1. Open terminal window, and change to the directory with the fastq files
2. create a directory called:  `fastqc-results` with the `mkdir` command.
3. next we run the command: `fastqc -t 2 -o fastqc_results *.fastq` , to understand thie flags in this command type: `fastqc --help` 
4. open firefox and open a fastqc results file: `fastqc_report.html`. You can find for each dataset a folder inside the folder: `fastqc-results`. For instance open the folder for the forward sequences: `F3D0_S188_L001_R1_001_fastqc` and check the results there.

Note that this is amplicon data, why is that different from normal shotgun data, and how is that observed in the fastqc file? See also these documents: [Good sequences](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) or [Bad sequences](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)

Next we open mothur an start working on the MiSeq SOP produced by Pat Schloss.

We work until the step where we remove chimeric sequences.

The commands can be found here: [20180205 MiSeq_SOP commands](20180205_MiSeq_SOP_commands.md).

**Recommended literature**

* Preclustering: Huse et al., Ironing out the wrinkles in the rare biosphere through improved OTU clusteringEnv. Microbiol 2010, 1889 – 1898
* 

### homework
Read the manuscript: `Madness of Microbiome` download [here](http://aem.asm.org/content/early/2018/01/29/AEM.02627-17.abstract).

## Lecture 12 - February - 2018

During this day we discus the manuscript from [Pollock et al., 2018](http://aem.asm.org/content/early/2018/01/29/AEM.02627-17.abstract).

And we continue with the clean-up of the MiSeq SOP amplicon data. You can find the details of the mothur commands in this document [20180212 MiSeq_SOP commands](20180212_MiSeq_SOP_commands.md).

**Recommended literature**

* [Clustering in mothur](http://blog.mothur.org/2016/01/12/mothur-and-qiime/#Clustering)
* Westcott and Schloss, Peerj 2015.  [De novo clustering methods outperform reference-based methods for assigning 16S rRNA gene sequences to operational taxonomic units](https://peerj.com/articles/1487/) DOI: 10.7717/peerj.1487
* Westcott and Schloss, mSphere 2017. [OptiClust, an Improved Method for Assigning Amplicon-Based Sequence Data to Operational Taxonomic Units](http://msphere.asm.org/content/2/2/e00073-17) DOI: 10.1128/mSphereDirect.00073-17
* Edgar, PeerJ 2017. [Accuracy of microbial community diversity estimated by closed- and open-reference OTUs](https://peerj.com/articles/3889/) DOI: 10.7717/peerj.3889

### homework
Now that we have finished cleaning our dataset we can now put everything into a textfile, so we can analyse the data using the abel computing cluster.

Create a text file with all the mothur commands to clean-up your datafile.

You can test if this file works in your biolinux machine. The way to call the script is like this: `mothur yourfile_with_commands.txt`

For a little bit more info see the Mothur manual on [batch mode](https://www.mothur.org/wiki/Batch_mode)

Next I want you to create a slurm script that can will run mothur with thr script on abel.

Use the slurm commands you learned during the introduction course to set-up a computing allocation on able that will run your mothur script with the commands. For more on that see the [abel user guide](http://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/).

Make sure you ask for the right mothur module and call the script with your commands from your working directory on abel.

## Lecture 26 - February - 2018

We will start with the homework of the last time:

* Prepared questions
* script with mothur commands
* slurm script to run the mothur on abel

After that we will start with diversity analysis using mothur. Today we focus on alpha diversity and how to vizualize community evenness using R-studio.

**Recommended literature**

* Lozupone & Knight, 2005 [UniFrac: a New Phylogenetic Method for Comparing Microbial Communities](http://aem.asm.org/content/71/12/8228.long)
* Lozupone et al., 2007 [Quantitative and Qualitative β Diversity Measures Lead to Different Insights into Factors That Structure Microbial Communities.](http://aem.asm.org/content/73/5/1576.long)
* Ramette 2007 [Multivariate analyses in microbial ecology](https://academic.oup.com/femsec/article/62/2/142/434668) DOI:10.1111/j.1574-6941.2007.00375.x
* Buttigieg & Ramette 2014 [A guide to statistical analysis in microbial ecology: a community-focused, living review of multivariate data analyses.](http://onlinelibrary.wiley.com/doi/10.1111/1574-6941.12437/abstract;jsessionid=C7AF57B5E53898E4F0EEA4E177CCF2F7.f02t03) DOI:10.1111/1574-6941.12437
* Website: [GUide to STatistical Analysis in Microbial Ecology (GUSTA ME)](https://mb3is.megx.net/gustame)
* R-Tutorial (pdf) [Multivariate Analysis of Ecological
Communities in R: vegan tutorial] (http://cc.oulu.fi/~jarioksa/opetus/metodi/vegantutor.pdf)



 
## Future Lectures
* [Phyloseq tutorial](http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html)
