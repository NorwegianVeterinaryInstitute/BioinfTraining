# Access to computing resources: VIA qlogin and SLURM script

## Overview: requesting computing resources

Please read this: [Introduction to the queuing system](https://www.uio.no/english/services/it/research/events/2018b/abel_intro_march2018.pdf). 

A little summary:

 **SLURM** is the queue system and scheduler daemon (program/process which always runs in the background on Abel). It allows optimizing computing resources and schedule when jobs should run.

There are 2 ways to request for computing resources:

### SLURM script

1. By submitting a jobscript using `sbatch script_file.slurm` command.

The script is a list of commands (bash) that will be interpreted one after another and are necessary for your analyses to be performed.
The scripts will be scheduled to run on a computing node when resources are availables.

A script contains:
  - a description of the computing resources you ask for (time, memory..) and which project (nn9305k) should be charged for those computing resources
  - define transferred files (in/out) of the temporary folder created in the computing node when resources are allocated
  - what programs are run, which parameters are used, in which succession ...

- using slurm script is particularly appropriate when you have **heavy jobs** that need to run for a long time AND/OR
- want to run a serie of analyses one after the other (automatizing) without having to wait that one finishes before launching the next one.

### qlogin

2. By requesting resources to the queue system with the `qlogin` command. 

Doing so will provide you with an interactive way of using computing resources.

  - arguments in `qlogin` specify the description of computing resources you ask for (similar as done via slurm script) BUT:

When resources are attributed:
 - you actually are logged on a **computing node**(*) that now executes at once commands you type directly from the shell (Unless you misspelled the syntax!)
 - This is perfect for running small tests
 - It can also be used to run longer jobs.
   > in this case use `screen` before asking resources with qlogin. If you do not know what screen is, read this [section](techstuff.md#using-screen)

What you should do when you use `qlogin`
- use `hostname` to check which node you are on (be sure you are on a computing node and not a login node)
 > ex: check before and after asking for resources
- look at the `job_id` when you ask for resources (and find out how to do it when you forgot to look - see commands under)

- you need to move files you want to use (not databases) TO `$USERWORK` = `work/users/<username>` and FROM after your analyses are completed.

### Oops!
Whether you use a SLURM script or use resources via qlogin you will need to load required modules:
`module load module_name`

## Relevant commands: + viewing status, canceling jobs 

NB: you can view status: `PD	Pending`, `R	Running`, `CD	Completed` ... See [Abel Queue system] for other statuses.

See commands below.

| command     | relevant arguments            | does                    | example                       |
|:------------|:------------------------------|:------------------------|:------------------------------|
|qlogin       | <account, time, ntasks, mem > |ask for interactive computing resources | `qlogin --account=nn9305k --time=00:30:00 --ntasks=4 --mem-per-cpu=4G` |
|sbatch       | <script_name>                 |submit a script to slurm | `sbatch sample_slurm.slurm` | |
|squeue       |-u <your_user_name>            |checks your current queue status  - provide jobID | |
|scancel      |<jobid >                       |cancel job with the provided id  | |
|hostname     |                               |checks where your are on abel: `login node` `or computing node` | |


## Examples: exercises
### 1. Login to Abel 
```
ssh your_user_name@abel.uio.no
```

### 2. qlogin 
```
qlogin --account=nn9305k --ntasks-per-node=16
```
* How can you see if you are in a login node or a computing node? (2 ways)

Answers to questions can be found [here](Quiz_answers.md#access-to-computing-resources) (but you need to try answering yourself first)

### 3. Prepare SLURM Script
```
mkdir QSystem_Test
cp /work/projects/nn9305k/samplefiles/SLURM_Script_BioinfoCourse.slurm QSystem_Test/Trim.sh
```
Add Trimmomatic commands to SLURM script Trim.sh

``` 
# For Trimmomataic
module load java

$input_file1="/work/projects/nn9305k/samplefiles/Test1.fastq"
$input_file=2"/work/projects/nn9305k/samplefiles/Test2.fastq"

$output_file1="/usit/abel/u1/jeevka/QSystem_Test/Trimmed_Test1.fastq"
$output_file2="/usit/abel/u1/jeevka/QSystem_Test/Trimmed_Unpaired_Test1.fastq"
$output_file3="/usit/abel/u1/jeevka/QSystem_TestT/rimmed_Test2.fastq"
$output_file4="/usit/abel/u1/jeevka/QSystem_Test/Trimmed_Unpaired_Test2.fastq"


/work/projects/nn9305k/bin/trimmomatic-0.36.jar PE -threads 1 -trimlog vibrio_trimlog.log $input_file1 $input_file2 $output_file1 $output_file2 $output_file4 $output_file4 ILLUMINACLIP:adapter.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:20:15 MINLEN:60
```

### Submit SLURM Job
```
sbatch Trim.sh

squeue -u your_user_name

squeue 

scancel <jobid>

scancel -u jeevka 
```

## Going further:

[Abel Queue system]
in [Abel User guide](https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/)

[Slurm documentation]:(https://slurm.schedmd.com/)

[Abel Queue system]:https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/queue-system.html
