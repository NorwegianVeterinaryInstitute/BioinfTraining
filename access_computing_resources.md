# Access to computing resources: VIA squeue and slurm script

## Overview: requesting computing resources

Please read this: [Introduction to the queuing system](https://www.uio.no/english/services/it/research/events/2018b/abel_intro_march2018.pdf). 

A little summary:

 **SLURM** is the queue system and scheduler daemon (program/process which always runs in the background on Abel). It allows optimizing computing resources and schedule when jobs should run.

There are 2 ways to request for computing resources:

### Slurm script

1. By submitting a jobscript using `sbatch script_file.slurm` command.

The script is a list of commands (bash) that will be interpreted one after another and are necessary for your analyses to be performed.
The scripts will be scheduled to run on a computing node when resources are availables.

A script contains:
  - a description of the computing resources you ask for (time, memory..) and which project (nn9305k) should be charged for those computing resources
  - where on the disk are the files you want to analyses and programs you want to use are located
  - define transferred files (in/out) of the temporary folder created in the computing node when resources are allocated
  - what programs are run, which parameters are used, in which succession ...

- using slurm script is particularly appropriate when you have **heavy jobs** that need to run for a long time AND/OR
- want to run a serie of analyses one after the other (automatizing) without having to wait that one finishes before launching the next one.

## qLogin

2. By requesting resources to the queue system with the `qlogin` command. 

Doing so will provide you with an interactive way of using computing resources.

  - arguments in `qlogin` specify the description of computing resources you ask for (similar as done via slurm script) BUT:

When resources are attributed:
 - you actually are logged on a **computing node**(*) that now executes at once commands you type directly from the shell (Unless you misspelled the syntax!)
 - This is perfect for running small tests
 - It can also be used to run longer jobs.
   > in this case use `screen` before asking resources with qlogin. If you do not know what screen is, read the section in [techstuff] 

What you should do when you use `qlogin`
- use `hostname` to check which node you are on (be sure you are on a computing node and not a login node)
 > ex: check before and after asking for resources
- look at the `job_id` when you ask for resources (and find out how to do it when you forgot to look - see commands under)

- you need to move files you want to use (not databases) TO `$USERWORK` = `work/users/<username>` and FROM after your analyses are completed.

## Viewing status, canceling jobs

NB: you can view status: `PD	Pending`, `R	Running`, `CD	Completed` ... See [Abel Queue system] for other statuses.

See commands below.

### Relevant commands
| command     | relevant arguments            | does                    | example                       |
|:------------|:------------------------------|:------------------------|:------------------------------|
|qlogin       | <account, time, ntasks, mem > |ask for interactive computing resources | `qlogin --account=nn9305k --time=00:30:00 --ntasks=4 --mem-per-cpu=4G` |
|sbatch       | <script_name>                 |submit a script to slurm | `sbatch sample_slurm.slurm` | |
|squeue       |-u <your_user_name>            |checks your current queue status  - provide jobID | |
|scancel      |<jobid >                       |cancel job with the provided id  | |
|hostname     |                               |checks where your are on abel: `login node` `or computing node` | |
|screen       |                               |avoid dropped connections - see [techstuff] | |
|...^a+d      |                               |detaches screen -> softwares keep running if you close shell| |
|screen       | -ls                           |gives list of your running screens  & screen_id| |
|screen       |-DR <screen_id>                |reconnect to screen with id..| |
|...^a+k      |                               |kills (close the attached screen) | |

... means when you are in screen press `ctrl = ^` + keys `a+d` or `a+k`

## Examples: exercises

1. How can you see if you are in a login node or a computing node? (2 ways)

> Answers to [Access to computing resources] questions (but you need to try answering yourself first)


## Going further:

[Abel Queue system]
in [Abel User guide](https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/)

[Slurm documentation](https://slurm.schedmd.com/)

[Abel Queue system]:https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/queue-system.html

[techstuff]:https://github.com/NorwegianVeterinaryInstitute/BioinfTraining/blob/master/techstuff.md
