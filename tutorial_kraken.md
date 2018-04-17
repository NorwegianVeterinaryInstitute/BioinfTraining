# Introduction to classification with Kraken

**Goal**: Introduce the tool kraken, and how to use it on the computing cluster Abel.

**Kraken** is used for classification of shotgun reads. The reads can be from any shot-gun project, genomic or metagenomics.
A thing to realize is that kraken is  a very memory hungry tool. The normal kraken bacterial genome database (based on Refseq) has the size of >170 Gb and that needs to be loaded into memory to function correctly. So Kraken will run on abel most efficiently when using the hugemem nodes.
The kraken database contains all the unique kmer sequences found in the entire Refseq database, and those Kmers are matched to a taxon level. Some Kmers are unique for strains, while others are found in entire phyla. That is how the classification works.

**Modules** 

Kraken is loaded into memory using the module system on abel.
In order to load the modules add the following lines to the file: `bash_login`.

```
# adding software and modules found in /work/projects/nn9305k/bin 
export PATH=/work/projects/nn9305k/bin:$PATH
umask 0002
module use /work/projects/nn9305k/bin/modules
```
After adding these lines to the `bash_login` file, logout of abel and log in again. Now type: 
```
module avail kraken
```
When you see the following it works:

```
------ /work/projects/nn9305k/bin/modules ---------
kraken/1.0
```

**Database**

Database folder:  `minikraken_20171013_4GB`
located in folder: `/work/projects/nn9305k/db_flatfiles/kraken_databases/minikraken_20171013_4GB`

In this folder we find the mini kraken database. That database only contains a subset of the kmers (≈ 2-7 %). This 4GB database constructed from "dustmasked" bacterial, archaeal, and viral genomes in Refseq as of Oct. 18, 2017.

**Test data**.

To do a proper test, I have collected in one folder the `raw` fastq sequences of a *B.anthracis*  sequencing project. The total dataset has about 0.5 milion reads.
Test sequences folder:  `/work/projects/nn9305k/tmp/kraken_test`
Test sequence filenames: `Banthracis_seq_R1_001.fastq.gz` & `Banthracis_seq_R2_001.fastq.gz`
___

# Classification of test sequences
We are going to run kraken and classify test sequences using the interactive set-up available on abel.

#### setting up a and interactive job
login to abel and type on the command line: 

```
screen
```
Next we are going to set-up a qlogin with only 2 cpus and with a total RAM memory of 3800M * 2 = 7600M

```
qlogin --account=nn9305k --mem-per-cpu=3800M --cpus-per-task=2 --time=1:0:0
```
When your job is created then type to activate the job:

```
source /cluster/bin/jobsetup
```

#### loading the module for the tutorial
Next do we need to create our environment to run kraken. We first load the module kraken and then the module kronatools. So type on the commandline:

```
module load kraken/1.0
module load kronatools/2.7
module list   ## this will show the loaded modules
```

#### copy test data to home folder.
To test kraken we first need to copy the test data to a directory where we can work.

My directory is called: `kraken_test`

```
rsync -rv /work/projects/nn9305k/tmp/kraken_test/Banthracis_seq_R* ./
```
Now we are ready to run kraken from this test folder.

### Running kraken on the test data

let's type :`kraken`

Take a good look at the options you see now, and note which options you needs to use.

Let's build the kraken command:

```
kraken -db /work/projects/nn9305k/db_flatfiles/kraken_databases/minikraken_20171013_4GB \
--threads 2 \
--fastq-input \
--gzip-compressed \
--classified-out Ba_K1.classified.fastq \
--unclassified-out Ba_K1.unclassified.fastq \
--output Ba_K1.Kr_out.txt \
--preload \
--paired \
 Banthracis_seq_R1_001.fastq.gz \
 Banthracis_seq_R2_001.fastq.gz
```

Running this command with this database takes about 1 minute or so.

Now we have several files present in our directory. Let's look at them:

```
bash-4.1$ ls -l 
total 696740
-rw-rw-r-- 1 thhaverk nn9305k 160668069 Apr 16 15:42 Ba_K1.classified.fastq
-rw-rw-r-- 1 thhaverk nn9305k  69139255 Apr 16 15:42 Ba_K1.Kr_out.txt
-rw-rw-r-- 1 thhaverk nn9305k 201272106 Apr 16 15:42 Ba_K1.unclassified.fastq
-rwxrwx--- 1 thhaverk nn9305k 131638770 Apr 16 15:26 Banthracis_seq_R1_001.fastq.gz
-rwxrwx--- 1 thhaverk nn9305k 150742603 Apr 16 15:26 Banthracis_seq_R2_001.fastq.gz
```

We got two fastq files, with the classified and unclassied reads, and one file with the classifications: `Ba_K1.Kr_out.txt`.

Let's take a better look at that file with the command: `head Ba_K1.Kr_out.txt `.

```
C       M01132:116:000000000-ABE3F:1:1102:15814:1030    1386    603     A:1 1386:1 0:3 1386:1 0:41 A:526
U       M01132:116:000000000-ABE3F:1:1102:20200:1030    0       603     A:1 0:46 A:526
C       M01132:116:000000000-ABE3F:1:1102:9689:1030     1392    603     A:1 0:22 1392:1 0:23 A:526
U       M01132:116:000000000-ABE3F:1:1102:11068:1031    0       603     A:1 0:46 A:526
C       M01132:116:000000000-ABE3F:1:1102:20976:1032    86661   603     A:1 0:18 86661:1 0:11 86661:1 0:14 86661:1 A:526
U       M01132:116:000000000-ABE3F:1:1102:15402:1032    0       603     A:1 0:46 A:526
U       M01132:116:000000000-ABE3F:1:1102:14313:1033    0       603     0:47 A:526
U       M01132:116:000000000-ABE3F:1:1102:20557:1033    0       603     0:47 A:526
U       M01132:116:000000000-ABE3F:1:1102:16410:1033    0       603     0:47 A:526
U       M01132:116:000000000-ABE3F:1:1102:19939:1033    0       603     0:47 A:526
```
Each sequence classified by Kraken results in a single line of output. Output lines contain five tab-delimited fields; from left to right, they are:

1.   "C"/"U": one letter code indicating that the sequence was either classified or unclassified.
2. The sequence ID, obtained from the FASTA/FASTQ header.
3. The taxonomy ID Kraken used to label the sequence; this is 0 if the sequence is unclassified.
4. The length of the sequence in bp.
5. A space-delimited list indicating the LCA mapping of each k-mer in the sequence. 
	
For example, `A:1 1386:1 0:3 1386:1 0:41 A:526` would indicate that:
*  the first 1 k-mers contained an ambiguous nucleotide
*  the next 1 k-mers mapped to taxonomy ID #1386
*  the next 3 k-mers  
*  the next  k-mers mapped to taxonomy ID #1386
*  the next 41 k-mers were not in the database
*  the last 526 k-mers contained an ambiguous nucleotide

Note that this looks like a bad quality sequence --> So QC of your sequences is important for getting good classifications.

Let's get the proper names for all the classifications that we got.

```
kraken-translate --db /work/projects/nn9305k/db_flatfiles/kraken_databases/minikraken_20171013_4GB Ba_K1.Kr_out.txt > Ba_K1.Kr_out.labels.txt
```

That creates the file `Ba_K1.Kr_out.labels.txt`, which contains:

```
M01132:116:000000000-ABE3F:1:1102:15814:1030    root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus
M01132:116:000000000-ABE3F:1:1102:9689:1030     root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus cereus group;Bacillus anthracis
M01132:116:000000000-ABE3F:1:1102:20976:1032    root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus cereus group
M01132:116:000000000-ABE3F:1:1102:18528:1033    root;cellular organisms;Bacteria;Terrabacteria group;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;Bacillus cereus group
```
So now we know for all classified reads the classification. The next step is to generate a report of the kraken output, but we first might want to filter the classifications to get higher confident classifications. The command to do that is `kraken-filter` and the confidence threshold can be set between 0 and 1.

```
kraken-filter --db /work/projects/nn9305k/db_flatfiles/kraken_databases/minikraken_20171013_4GB --threshold 0 Ba_K1.Kr_out.txt > Ba_K1.Kr_out.conf_0.txt
```
The output looks like this:

```
C       M01132:116:000000000-ABE3F:1:1102:15814:1030    1386    603     P=0.043 A:1 1386:1 0:3 1386:1 0:41 A:526
U       M01132:116:000000000-ABE3F:1:1102:20200:1030    0       603     P=0.000 A:1 0:46 A:526
C       M01132:116:000000000-ABE3F:1:1102:9689:1030     1392    603     P=0.022 A:1 0:22 1392:1 0:23 A:526
U       M01132:116:000000000-ABE3F:1:1102:11068:1031    0       603     P=0.000 A:1 0:46 A:526
C       M01132:116:000000000-ABE3F:1:1102:20976:1032    86661   603     P=0.065 A:1 0:18 86661:1 0:11 86661:1 0:14 86661:1 A:526
```
Where we can see that P values are added, and based on our threshold C's are turned into U's. When you increase the threshold you will see that the kraken classification ID become 0 when a classification does not reach the threshold.

Now we have a table that can be used to build a report with counts per taxon.

```
kraken-report --db /work/projects/nn9305k/db_flatfiles/kraken_databases/minikraken_20171013_4GB Ba_K1.Kr_out.conf_0.txt > Ba_K1.Kr_out.conf_0.classifications.txt
```
This produces a count table:

```
55.61  309467  309467  U       0       unclassified
 44.39  247036  3216    -       1       root
 41.33  230018  0       -       131567    cellular organisms
 41.33  230016  2863    D       2           Bacteria
 40.76  226815  4       -       1783272       Terrabacteria group
 40.76  226807  10      P       1239            Firmicutes
 40.75  226797  190     C       91061             Bacilli
 40.72  226600  299     O       1385                Bacillales
 40.66  226294  6       F       186817                Bacillaceae
 40.66  226284  10863   G       1386                    Bacillus
 38.71  215396  170433  -       86661                     Bacillus cereus group
  7.76  43168   43129   S       1392                        Bacillus anthracis
  0.00  24      24      -       568206                        Bacillus anthracis str. CDC 684
  0.00  5       5       -       768494                        Bacillus anthracis str. H9401
  0.00  3       3       -       261591                        Bacillus anthracis str. Vollum
  0.00  3       3       -       1452727                       Bacillus anthracis str. Turkey32
  0.00  2       2       -       1412844                       Bacillus anthracis 52-G
  0.00  1       1       -       198094                        Bacillus anthracis str. Ames
  0.00  1       1       -       1412843                       Bacillus anthracis 8903-G
  0.22  1221    513     S       1396                        Bacillus cereus
  0.02  91      91      -       526968                        Bacillus cereus R309803
```
  
This shows that most of the reads in this dataset are classified to B.anthracis. Scroll through the results to see what else is found in the dataset. 

The final step is to create a report in the metaphlan style, which could be useful for comparisons.

```
kraken-mpa-report --db /work/projects/nn9305k/db_flatfiles/kraken_databases/minikraken_20171013_4GB Ba_K1.Kr_out.conf_0.txt > Ba_K1.Kr_out.conf_0.classifications.mpa.txt
```	

## Creating a krona plot
Krona plots are easy explorative plots that can be viewed in a web-browser. With the software kronatools we can create these from tabular data files. We can use our kraken output files to generate such a plot.

The first step is to extract the two columns of the file `Ba_K1.Kr_out.conf_0.txt` that are the sequence name and the classifications.
```
cat Ba_K1.Kr_out.conf_0.txt | cut -f2,3 > Ba_K1.Kr_out.conf_0.krona
```

And then we make a plot: 

```
ktImportTaxonomy Ba_K1.Kr_out.conf_0.krona -o krona.Ba_K1.html
```

Open a new terminal and use that to download the krona html page. First go to the location on your laptop where you wan to store the figure.

```
 rsync -rv abel.uio.no:/work/projects/nn9305k/home/thhaverk/kraken_test/krona.Ba_K1.htm* ./
```

and open the html page with your webbrowser.


Now we are done with generating our report for this one dataset.

You can run this also a a slurm script. See the folder `/work/projects/nn9305k/bin/kraken-1.0/slurm_scripts` where you find a script called: `kraken_classification_single_dataset.slurm`

If you only have one dataset, copy this script to your home area and modify it so that you can run the script on the job of interest.

___

# Running kraken on multiple samples
Kraken processes one sample at the time. If you want to analyze multiple samples, you have to split the process into multiple slurm jobs. In that way we can run jobs parallel and not after each other. On abel it is possible to start as a single user as many jobs as you like with an upper limit of ≈ 400 cpus being used.
In the folder:  `/work/projects/nn9305k/bin/kraken-1.0/slurm_scripts` we find multiple script to do this, which both use an array, but in a slightly different manner.

**Array Method 1**:

 Using a single script to start processing multiple jobs. The first script spawns new jobs, until all datasets have been processed. The example script `kraken_array_script_20180410.slurm` is an example of this.
 
 In order to use this method you simple add an SBATCH command line to the slurm file. It has to be the topline, like here:

```
    #SBATCH --array=0-9 
    #SBATCH --account=nn9305k
    #SBATCH --time=1:00:00
```
The top line needs to specify the number of jobs to run minus one. here is asks for 10 seperate jobs.

**Array Method 2**: 

Using a master script and a worker script to run the analysis on multiple jobs. The master script runs a small job with little requirements and that starts new jobs using the workerscript. The workscript is the heavy duty job and does the actual analysis.
Here the files are:
master script: array_submission_kraken_20180411.slurm.
worker script: workscript_kraken_20180411.slurm.

Both methods have their pro's and cons. The later option is prefered when you want to combine all the results from all your jobs after the last sample analysis has finished. If that is not important, than it might be easier to use the first option. Another advantage, is that you never have to set the number of jobs required.

To use this method the master script needs to contain the line with the `arrayrun` command. In the method below each sample is in a seperate folder collected in one main directory.

```
# collecting all the samples into an array called FOLDERS 
FOLDERS=($(ls -1))

# Than we count the number of folders (samples) to determine the number of tasks
NUM_TASKS=${#FOLDERS[*]}

# determine the maximum number of tasks
MAX_ID=$((NUM_TASKS - 1))

# running the arrayrun command with tasks id's ranging from 0-MAX_ID
arrayrun 0-$MAX_ID /work/users/thhaverk/slurm_scripts/workscript_kraken_20180411.slurm

```

The workscript needs to contain the following files to start the correct analysis of a specific sample.

```
FOLDERS=($(ls -1))   ##  FILES contains all the files for analysis
MYFOLDER=${FOLDERS[$TASK_ID]}  ## identify which dataset is analyzed
```


___


# Building a kraken database
Kraken databases are build with the commands found in the manual and this is the order that needs to be taken.

A slurm script can be found in the folder: `/work/projects/nn9305k/bin/kraken-1.0/slurm_scripts`.
The script is called: `kraken_build.slurm`

```
#### Load necessary modules
module load kraken/1.0

#### Download taxonomy
kraken-build --download-taxonomy --db kraken_db

#### Download and install standard genomic libraries
echo Downloading standard Kraken libraries

kraken-build --download-library plasmid --db kraken_db
kraken-build --download-library viral --db kraken_db
kraken-build --download-library archaea --db kraken_db
kraken-build --download-library bacteria --db kraken_db

echo Done with downloading and installing standard Kraken libraries
```

We could aslo add protist or other eukaryote genomes to the database.
For instance for protists you need to run these commands.

```
## Download protozoan and fungal genomes (no need to call perl, done in module)
echo "Downloading custom libraries (fungi and protozoa)"

## Download the taxonomy.
kraken-build --download-taxonomy --db kraken_db

## scripts found in the kraken directory.
$DOWNLOAD_SCRIPTS/download_fungi.pl
$DOWNLOAD_SCRIPTS/download_protozoa.pl

echo Done with downloading fungal and protozoan genomes


## Add protozoa and fungi to library
for dir in fungi protozoa; do
        for fna in `ls $dir/*.fna`; do
                kraken-build --add-to-library $fna --db kraken_db
        done
done
```

After all this it is time to build the database. That is a long and memory intensive job, so it is good to use a slurm script to run this process. Do not do it on a qlogin without using `screen`.

**Memory requirements:**

```
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=30000M
#SBATCH --cpus-per-task=8
#SBATCH --partition=hugemem
```
≈ 240 Gb of memory needed to build normal database. When adding larger genomes, more memory might be required.

```
## Build the actual database (from libraries)
kraken-build --build --db kraken_db --threads 8

## Remove all intermediate files
kraken-build --clean --db kraken_db --threads 8
```

Now a complete database is created.

___

# Example slurm scripts.

Below you find example slurm script to run kraken with one or more samples and to build a kraken database

* [Single sample slurm script](metagenomics_scripts/kraken_classification_single_dataset.slurm)
* [Array method 1 slurm script](metagenomics_scripts/kraken_array_script_20180410.slurm)
* [Array method 2 MASTER slurm script](metagenomics_scripts/array_submission_kraken_20180411.slurm)
* [Array method 2 WORKER slurm script](metagenomics_scripts/workscript_kraken_20180411.slurm)
* [Build a large kraken database slurm script](metagenomics_scripts/kraken_build.slurm)

