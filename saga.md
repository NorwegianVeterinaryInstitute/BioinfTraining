
# Practical differences Abel - Saga

- Combining storage: NIRD and computing power: SAGA
  - SAGA will be located in Trondheim. NIRD storage both in Tronheim and Tromsø (right?)
  - SAGA will be used for computing power only: at first there will be **NO BACKUP** from data hosted on Saga. It is planned to implement "slow file system where we could have some files at a later stage"

- **Disk space on Saga is much more limited on Saga than on Abel.** This means that we cannot have all our data permanently stored on Saga - only data required for ongoing analyses over some restricted period of time will have to be on Saga. 
> _in practice_: if you are finished with a certain type of data analysis (ie. your assembly), the the raw reads that you COPIED from NIRD to Saga to perform you analysis will have to be deleted/moved, and all your finished analyses will have to be moved to NIRD. If you do not do that, we might run out of space and nobody will be able to work

- **NIRD** data backup system is different than what was on Abel. On Abel, data were backed-up on tape. On NIRD this will be screenshot. 
- [ ] jeevan?
> if I understand well -> version control - Question: Does this have practical implications for users (said it will be more easy to messup)

- you will have the possibility to see your data directly from NIRD and Saga. (But this is not implemented yet). And you will **never be able to use data located on NIRD as input files**. 
  > This implies that you will have to transfer data back and forth from **NIRD** storage area to Saga to be able to do your analyses.
  
- **Disk areas on SAGA:**

$USER is your username variable 

  - `$HOME` or `/cluster/home/$USER` - 20GB - 
  - `$USERWORK` or `/cluster/work/users/$USER`
  - `$SCRATCH`-> created when job is running (same as Abel)
  - PROJECTS `/cluster/projects/nn9305k` -> for shared things
  - Per today there is NO /WORK area -> but this is planned at a later stage
  
  > $USERWORK and PROJECTS share the same file system - fast disks - no permanent storage. 
  - It is recommended to use $USERWORK for using computing resources (unless you have things that needs a lot of read/write there you should arrange that you use $SCRATCH, copying files to and back from $SCRATCH in your scripts).
  - There is actually no really difference between using $USERWORK and PROJECT for using computing resources. BUT files in $USERWORK will be deleted automatically after 20-42 days depending on usage. There is a information message when you log on SAGA (no email warning); and using $USERWORK will force you to clean your data regularly. 
  
- **Disk areas on NIRD:**

- NIRD HOME: `nird/home/<username>`. Quota 20GB and **100000 files**. 
- NIRD PROJECTS: `/nird/projects/NS9305K` VI quotat = 20TB. (also number of file limits)
- There is a small $SCRATCH `/scratch/username` area on NIRD -> 15TB 

Per today there is NO mounting between SAGA and NIRD (meaning you cant see your files that are on NIRD from SAGA). This will be fixed in the future. But anyway you wont be able to use NIRD files as input for your analyses. 


! add

 - [ ] does archives tar file of files counts for one or several? I guess one? 
  
- [ ] jeevan, did I get this right? 
- There is a part on Saga that can be reserved for bioinformatic, part reserved for machine learning on SAGA (ie. if you need better floating point for your analyses)

- Queue system in SAGA is **SLURM** - There are minor differences with Abel. Most important is for requesting **interactive ressources**. Please see bellow. 

# Extra features:
- It will be _easy_ to customize your own modules (ie. add libraries to python, R)
- You can have access to a graphical interface (ie. working with Rstudio, coding with python, see graphs you created ...)


# Administrators need to fix:
- find a way to have a common Conda environment 


# Where to find more detailed information

[SAGA documentation](https://documentation.sigma2.no/quick/saga.html)

[NIRD documentation](https://documentation.sigma2.no/storage/nird.html)

[Slurm documentation](https://slurm.schedmd.com/documentation.html)

[Password changing](https://www.metacenter.no/user/login/)

[Abel to Saga course link](https://www.uio.no/english/services/it/research/events/Migrate-to-HPC.html) - includes HPC course.

[HPC course](https://sabryr.github.io/hpc-intro/)

Support: `support@metacenter.no` indicate which resource. See informations in documentation

# Summary: how to: main commands

## VI project ID: 

`nn9305k` on SAGA, `NS9305K`on NIRD

## Login 

To saga 
```bash
ssh YOUR_USERNAME@saga.sigma2.no
```
To NIRD 
```bash
ssh YOUR_USERNAME@login.nird.sigma2.no
```

Check where you are: `hostname`


## Checkings disk usage for your quotas

SAGA - [ ] dusage does not work
- SAGA HOME: `du -h $HOME
- SAGA PRROJECTS: `du -h 

- NIRD HOME: `dusage`
- NIRD PROJECTS: `dusage -p NS9305K`

## Transfer data 

### Practical
Use `screen` to transfer data (to avoid that your data stops transfering when you shut down the window). See [techstuff](https://github.com/NorwegianVeterinaryInstitute/BioinfTraining/blob/master/techstuff.md)

As space both on SAGA and NIRD will be limited, and not everything will be totally functional at the beginning, a good organisation is vital. Therefore: please clean and sort your files before transfer of your data from ABEL to NIRD.

### One or fiew files: `scp <source> <destination>`
> Not for your initial transfer data from Abel to NIRD, but you can use that afterwards to tranfer between SAGA and NIRD

Examples: 
```bash
#path for projects: on NIRD /projects/NS9305K/

# For file
scp my_file.tar.gz <username>@login.nird.sigma2.no:/path
# For directory 
scp -r my_dir/ <username>@login.nird.sigma2.no:/path
```

### rsync: many files - not too many data
Typical transfer rate 4TB/day ~50MB/s 

Do it twice, to be sure all data has been transfered. Defined blocks that you are sure can be transfered, and check the log.

- [ ] Write here difference rsync / and not / in dest

Please do a little test with `--dry-run` to `rsync` : this will simulate the run without copying files: then you can check if they arrive in the folder you chose. There is a difference using `folder/` and `folder`.

When you are sure that files will arrive in the correct directory, you can remove `--dry-run` and launch the real moving of files. 

```bash
# IF you have many small files (ie. text files) compress:
tar -czvf archive_name.tar.gz /directory/data

#general usage: transfer
rsync -avxh <source> <destination_account/>
rsync -rauPWD <source> <destination> 

# Example: Logged on abel 
rsync -rauPWD /work/projects/nn9305k/yruck  <username>@login.nird.sigma2.no:/nird/projects/NS9305K

# Example: Logged on NIRD -> transfer SAGA 
rsync -rauPWD <username>@login.nird.sigma2.no:/nird/projects/NS9305K/path_myfiles \       <username>@saga.sigma2.no:/cluster/projects/nn9305k/mypath

# Example: Logged on SAGA -> transfer NIRD 
rsync -rauPWD <username>@saga.sigma2.no:/cluster/projects/nn9305k/mypath \ <username>@login.nird.sigma2.no:/nird/projects/NS9305K/path_myfiles        

# options -> we need to choose 
## OR: compress files on the flight 
..... -z 
## reates a directory on target where it puts files if something wrong appends -> can help to resume faster. 
......–partialdir=.partial-dir -> c

## Decompress archive
tar -xzvf archive_name.tar.gz
```

- [ ] Delete the file/folder in Abel. Be carrefull !

## parsyncfp : Not too many files but a lot of data. Divides in several tasks.
Available on Abel. 

- [ ] not tested

```
parsyncfp --rsyncopts="-ax –ignore-existing --relative” --NP=4 --startdir="/path” \
data username@login.nird.sigma2.no:/path
```

NP is the number of threads. They did recommend NOT to use too many as this will likely not go faster and might bug otherwise.


# Using SAGA 

## Computes nodes
4 different types [] - find back the nice table he had 

When you submit a job you will have to ask for a specific compute node. By default please use `normal` or `devel` -> for small tests. 
- jeevan correct? 



## Modules load
Identical abel

## Install modules / compile

## Install additional libraries to packages

## Use graphical interface
You can use a graphical interface on SAGA. It wont be the same as using your laptope, but you can use it for example with Rstudio, Phython scripts or to view graphics you created on SAGA.

To do so you need to login as such `ssh -Y username@saga.sigma2.no` and type the command `xeyes` on your shell. You will see a window with 2 pairs of eyes: this is the graphical display. 

## Shared databases
### blast: 
`/cluster/shared/databases/blast/latest` 

- [ ] there should be some guide here for shared things

## QUEUE SYSTEM : SLURM 

Very similar to ABEL -> some minor difference in commands for slurm scripts, **different for interactive login**

- Submitting a job to SLUM queue system `sbatch scriptname.slurm`. Examples of script [link]() -[ ] finish
- Interactive loggin: `squeue`does not exist anymore. 

Example: `srun --ntasks=1 --mem=8G --time=00:30:00 --qos=devel --cpus-per-task=2 --account=NN9305K --pty bash -i`

- [ ] check what it did correspond to
`--qos=` is the name of the computing node, chose between `devel` for small tests, `normal` ? was devel -> difference

--`partition=bigmem` for big jobs (see documentation)

Main commands to manage your jobs:

| command | does | additional information |
|:--------|:-----|:-----------------------|
| squeue -u $USER | show your jobs in queue | |
| scontrol show job <jobid> | show status of your job | [link]() |
| sacct -j <jobid> | show running status/info of OLD jobs (archive) |  
 


# Advanced SAGA commands

```bash
# check all the mounting points
df -h
# SAGA specific: indication load on SAGA - computing power nodes (pe for processor equivalence) 
freepe
```






