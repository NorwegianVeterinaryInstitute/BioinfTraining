
Abel Disk areas that are relevant for your work at the Veterinary Institute

# Summary :

|directory|path                |what      |characteristics|
|:--------|:-------------------|:---------|:--------------|
|`$HOME`|`usit/abel/u1/username`|"home-home" area|500GB - login area- NOT FOR Data Processing -backed-up|
|**work**|`/work/projects/nn9305k`|common resources (eg. softwares, databases)|10 TB - faster disc, closer to nodes connection|
|**project**|`/projects/nn9305k`|Project data area (backed-up, raw data, results)|30 TB - slower disc|
|`$SCRATCH`|`/work/jobs/JOBID.d`|computing nodes: **FOR PROCESSING DATA = RUNNING JOBS**|fastest|
|`$USERWORK`|`/work/users/username`|sharing data among different jobs = **staging**|do not backed-up - removed after 45 days|

1. **$HOME**: "home-home" directory (on login nodes): 500GB - backed up regularly: slow and costly.
  - If you need to store files that do not require backup: put them in a `nobackup` directory
  - recommended for: software, job configurations, job preparation, debugging ...
  - DO NOT USE to run jobs => long jobs should be run on $SCRATCH (submitted before hand to the queue system **SLURM**)

QQQ: -> **login nodes** NB: small tests are tolerated (processes taking more than 30 min are automatically killed) automatically <=> so in home-home

**work area**: `/work/projects/nn9305k` will be used to contain:
- softwares and libraries
- public reference data, template files, adapters files ...

2. **project area**: `/projects/nn9305k` will be used to contain:
 - project data, including raw sequencing files (zipped/tared:compressed and decompressed)
 - your `project/home/username`
 - data analyses that have been performed
 - ... this is where the work will NOW be done.

QQQ - [ ]-> **so our(vetinst users) /work/projects/nn9305k/home/username (incl. all subfolder should be moved here)** -> `project/home/username`

3. **$SCRATCH**: when a job is started, a temporary directory on $SCRATCH is automatically created. It is also automatically deleted when the job finishes. This allows the jobs to run faster without interfering with the work of others. You can access $SCRATCH for status monitoring of running jobs: with `/work/jobs/JOBID.d` where JOBID.d is the ID that has been provided when you submitted the job in the queue system SLURM.

4. **$USERWORK** If files are needed for more than one job, they are **staged** in here. Files are automatically deleted after 45 days, and there is no backup. Access: `/work/users/username`

For more details you can look at: [Managing Data on Abel] and more generally at [Abel User Guide]. You can also look at [Computer Resources at CEES].

#### [Workflow] between areas

<img src="https://docs.google.com/drawings/d/e/2PACX-1vSY_KCj3fubTH1zk6ZkOL6eLhoOOuAbp4bfu1YkOAvkadHhPfbuZrsepwHCUpEqwr45Zqt2hlEoCwVk/pub?w=960&amp;h=720">

NB: **SLURM** is the queue system and scheduler daemon (program/process which always runs in the background on Abel). It allows optimizing computing resources and schedule when jobs are to be run on $SCRATCH.

NB: There is also a **version of NCBI databases** hosted on Abel: at `/work/databases/bio/ncbi`

QQQ: Thought we could also have a summary checklist overview (of what was written in README f.eks.)-> what to do (at least easier for beginners..I guess)  bellow example

# MEMO: Guidelines:

```
QQQ
*this part should be checked-and modified after moving files in Abel*-> NEED to be sure areas ok
- [ ] the README should
- [ ] and New_user.txt
```

**New users**
- [ ] README in `/work/projects/nn9305k/`
- [ ] `/work/projects/nn9305k/samplefiles/new_user.txt` - follow instructions

**Everybody**
- [ ] read NEWS regularly: `/work/projects/nn9305k/`

#### All **projects directories** should be organized as such:

- [ ] : README copy from README_datafile - in `samplefiles` directory and filled about => rename to README_tarfilename
- [ ] directory structure
  - [ ] rawdata
  - [ ] analysis
  - [ ] scripts
  - [ ] logs
- [ ] please fill in the data registry file : `/projects/nn9305k/sys/DataRegistry.txt`

[Abel User Guide]:https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/
[Managing Data on Abel]:https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/data.html
[Workflow]:https://docs.google.com/drawings/d/e/2PACX-1vSY_KCj3fubTH1zk6ZkOL6eLhoOOuAbp4bfu1YkOAvkadHhPfbuZrsepwHCUpEqwr45Zqt2hlEoCwVk/pub?w=960&h=720
[Computer Resources at CEES]:https://github.com/uio-cees/hpc/wiki/Computer-resources
