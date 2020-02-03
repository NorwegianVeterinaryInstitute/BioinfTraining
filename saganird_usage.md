# Saga and NIRD

From 2020 we are now using Saga and NIRD as our compute platform. We now have a more divided structure, with NIRD as long time storage, and Saga as our compute resource. This means that we will have to work a bit differently than what we used on Abel.


## Login

We use ssh to get access to both computers. The location to ssh for the two computers are:

* Saga: saga.sigma2.no
* NIRD: login.nird.sigma2.no

In case you have issues with username/passwords, that can be fixed via [the Sigma2 Metacenter](https://metacenter.no).

The operational log for Saga and NIRD can be found here: [Operations log](https://opslog.sigma2.no/).

Full documentation about Saga and NIRD can be found here: [Sigma2 documentation](https://documentation.sigma2.no/).


## Disk areas

Today (Jan 2020) we have 40TB available in NIRD, and 10TB in Saga. This means we have to be a bit careful regarding where we have data.

Each main disk area mentioned below has a README file in it. Please read these for more info. The main directories are mentioned below. 

### NIRD

NIRD is our master storage facility. We don't do computations on NIRD, we just store the data that we are not currently processing on that server.

On NIRD you have access to two different directories:

* login directory: the location you get into when you log in. You have 20 GBs there. Don't keep Vetinst data there.
* /projects/NS9305k: the main storage directory for the Institute. The main directory to know about is the one named datasets.
* /projects/NS9305K/datasets: This is the main raw data directory on NIRD. This should be the master storage for your sequenced reads.

<!--
wgs, trancriptomics and metagenomics folder structure within datasets above is not explained.
--!>

### Saga

Saga is our main compute node. Computations should be done using slurm. On Saga there are two different disk areas that we have access to:

* login directory: this is the directory you get into just after having logged in. This is commonly called your "home directory". You have 20 GBs of space there. Do not store vetinst data there.
* /cluster/projects/nn9305k: this is the main work directory for us on Saga. There are two main directories under this directory that you should know about, this is datasets and active
* /cluster/projects/nn9305k/datasets: this is where you can put your _copies_ of your datasets that you have stored in NIRD. 
* /cluster/projects/nn9305k/active: this is the directory where each of us should do our work.


## Working on Saga

Each user should have a directory in the active directory on Saga which is named after your username (/cluster/projects/nn9305k/active/<USERNAME>). This is where each user should do their work. This lets the administrators see how much data each of us is generating. Inside your active (under <USERNAME>) you are free to organize your data in any way you please. 

All raw sequencing data (fastq files) should live in the datasets directory. To use the data in analysis, there are two different options. You can either use the absolute path to the data itself, or you can softlink the data into a directory inside your active directory.

Once you copy a sequencing data into the datasets directory (/cluster/projects/nn9305k/datasets), please update the datasets registration file. Include your username under the list of users. This will let us (and you) to see who is working on that particular dataset. 

Please note: all data in the datasets directory should be a _copy_ of what is in the NIRD datasets directory. In case we run out of space in a critical situation, we might end up having to delete some of the raw datasets. Hence the importance of the Saga data being only a copy.


## Dataflow NIRD and SAGA

NIRD is for more long term storage, while Saga shuld be reserved for data that we are actually currently working on. This means that we will have to transfer data between Saga and NIRD. When you have worked on a project and have finished with it, please pack it up (use gzip and tar) and transfer it to your directory on NIRD.

<!--
Where should the analysed data go in NIRD?
--!>



