
Abel Disk areas that are relevant for your work at the Veterinary Institute

# Summary :

|directory|path                |what      |characteristics|
|:--------|:-------------------|:---------|:--------------|
|`$HOME`|`usit/abel/u1/username`|"home-home" area|500 GB - login area - NOT FOR Data Processing - backed-up|
|**work**|`/work/projects/nn9305k`|common resources (eg. softwares, databases)|10 TB - faster disc, closer to nodes connection|
|**project**|`/projects/nn9305k`|Project data area (backed-up, raw data, results)|30 TB - slower disc|
|**project-home**|`projects/nn9305k/home/username`|Non-collaborative projects|subdirectory of **project area**|
|`$SCRATCH`|`/work/jobs/JOBID.d`|temporary area on computing nodes: **FOR PROCESSING DATA = RUNNING JOBS**|fastest|
|`$USERWORK`|`/work/users/username`|sharing data among different jobs = **staging**|do not backed-up - removed after 45 days|

## Storage Areas

1. **$HOME**: "home-home" directory (on login nodes): 500GB - backed up regularly: slow and timewise costly.
  - If you need to store files that do not require backup: put them in a `nobackup` directory
  - Recommended for: software, job configurations, job preparation, debugging ...
  - **All Veterinary Institute data you work on has to be accessible to the Veterinary Institute user-group. So you can’t store them in your "home-home"**.
 
2. **work area**: `/work/projects/nn9305k` will be used to contain:
  - Softwares and libraries
  - Public reference data, template files, adapters files ...

3. **project area**: `/projects/nn9305k` will be used to contain:
  - Project data, including raw sequencing files (zipped/tared:compressed and decompressed)
  - Data analyses that have been performed
  - ... this is where the work will NOW be done.
 
4. **project-home area**: `projects/nn9305k/home/username`: sub-directory of the Veterinary Institute **project area**. Here should all data from projects that do not require collaboration with other users go. NB: Data from collaborative projects must be placed in a directory directly under `project area`. 

<img src="https://docs.google.com/drawings/d/e/2PACX-1vRMKYHFeAQz-c03eQMKXXwKuUFxUao0Fe_2UHDBAZkLVcdGJzUon8nS6xW7DjeOBXfx3zJXkhgZCvxs/pub?w=960&amp;h=720">

**NB**: There is also a **version of NCBI databases** hosted on Abel: at `/work/databases/bio/ncbi`

--------------------------------------------------------------------------------------------------------------------------------

## Computing Areas (computing nodes)

5. **$SCRATCH**: when a job is started, a **temporary area** on $SCRATCH is automatically created. It is also automatically deleted when the job finishes. This allows the jobs to run faster without interfering with the work of others. You can access $SCRATCH for status monitoring of running jobs: with `/work/jobs/JOBID.d` where JOBID.d is the ID that has been provided when you submitted the job in the queue system SLURM (thus accessible only during job's lifetime).

6. **$USERWORK** If files are needed for more than one job, they are **staged** in here. Bifrost is configured to set files in $USERWORK by itself, BUT using some other softwares/pipelines might requires that you copy needed files here. You should run analysis from this area (via SLRUM script or other mechanisms). Path: `/work/users/username` OR type `$USERWORK` instead of path.
  - > NB: Abel managing system automatically delete files older than 45 days (those files are not backed-up). 
 
For more details you can look at: [Managing Data on Abel] and more generally at [Abel User Guide]. You can also look at [Computer Resources at CEES].

<img src="https://docs.google.com/drawings/d/e/2PACX-1vSEmwfZ_3Meo_GHKmRi0aaUK316j84oYEHy5qqDW-lXKR8wkgNjNUBchvDk9aLQllpN607Mq271g1uJ/pub?w=960&amp;h=720">

**SLURM** is the queue system and scheduler daemon (program/process which always runs in the background on Abel). It allows optimizing computing resources and schedule when jobs are to be run on $SCRATCH.- **All long jobs should be sumbited to to the queue system SLURM (and are run on $SCRATCH)**. 
 - There are 2 ways to use slurm SLURM queue system: 1. **interactively** with the command: `qlogin` and 2. by **submitting a script**: using `sbatch`. 
 - > when you use `qlogin` you actually are logged on a computing node that now executes commands from the shell. Usage: small tests (testing a working process OR using interactive programs).
 - > when you submit a script to SLURM, it will be scheduled to run on a computing node when resources are availables. Usage: heavy jobs. 
 - > **NB**: Small tests outside SLURM are tolerated (processes taking more than 30 min are automatically killed)

[Abel User Guide]:https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/
[Managing Data on Abel]:https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/data.html
[Workflow]:https://docs.google.com/drawings/d/e/2PACX-1vSY_KCj3fubTH1zk6ZkOL6eLhoOOuAbp4bfu1YkOAvkadHhPfbuZrsepwHCUpEqwr45Zqt2hlEoCwVk/pub?w=960&h=720
[Computer Resources at CEES]:https://github.com/uio-cees/hpc/wiki/Computer-resources
