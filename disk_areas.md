
Disc areas on Abel which are relevant for your work at the Veteriary Institute

# Summary :

|directory|path                |what      |characteristics|
|:--------|:-------------------|:---------|:--------------|
|`$HOME`|`usit/abel/u1/username`|"home-home" area|500GB - login area- NOT FOR DataProcessing -backed-up|
|**work**|`/work/projects/nn9305k`|common resources (eg. softwares, databases)|10 TB - faster disc, closer to nodes connection|
|**project**|`/projects/nn9305k`|Project data area (backed-up, raw data, results)|30 TB - slower disc|
|`$SCRATCH`|`/work/jobs/JOBID.d`|computing nodes: **FOR PROCESSING DATA = RUNNIng JOBS**|fastest|
|`$USERWORK`|`/work/users/username`|sharing data among different jobs = **staging**|do not backed-up - removed after 45 days|

1. **$HOME**: "home-home" directory (on login nodes): 500GB - backed up regularly: slow and costly.
  - If you need to store files which do not require backup: put them in a `nobackup` directory
  - recommended for: software, job configurations, job preparation, debugging ...
  - DO NOT USE to run jobs -> use $SCRATCH

QQQ: -> **login nodes** no big jobs but tolerated test small for jobs.-> processes more taking more than 30 min are killed automatically => so in home-home

**work area**: `/work/projects/nn9305k` will only be used to:
- contain programs and librairies for programs to work
- public reference data, template files, adapters files ...

2. **project area**: `/projects/nn9305k` that contains
 - project data, including raw sequencing files (zipped/tared and un-Z/T)
 - data analyses that are performed
 - this is where the work will NOW be done.

QQQ - [ ]-> **so our(vetinst users) /work/projects/nn9305k/home/username (incl. all subfolder should be moved here)**

3. **$SCRATCH**: when you start a job a temporary directory on SCRATCH is automatically created. It is automatically deleted when the job finishes. This allows the jobs to run faster without interfering with the work of others. You can access for status monitoring of your jobs: with `/work/jobs/JOBID.d`.

4. **$USERWORK** If files are needed for more than one job, they are staged in here. Files are deøeted after 45 days, and there is no backup. Access: `/work/users/username`

For more details you can look at: [Managing Data on Abel] and more generaly at [Abel User Guide]. You can also look at [Computer Ressources at CEES].
#### [Workflow] between areas

<img src="https://docs.google.com/drawings/d/e/2PACX-1vSY_KCj3fubTH1zk6ZkOL6eLhoOOuAbp4bfu1YkOAvkadHhPfbuZrsepwHCUpEqwr45Zqt2hlEoCwVk/pub?w=960&amp;h=720">

NB: **SLURM** is the queue system and shelduler deamon (program/process which always run in the background on Abel). It allows optimizing computing ressource usage and sheldule when runed on $SCRATCH.

NB: There is also a **version of NCBI databases** hosted on Abel: at `/work/databases/bio/ncbi`

QQQ: Thought we could also have a summary checklist overview -> what to do (at least easier for beginners..I guess)  beloow example
# MEMO: Guidelines:

```
QQQ
*this part shoud be checked-and modified after moving*-> NEEd to be sure areas ok before
- [ ] the README should
- [ ] and New_user.txt
```

**New users**
- [ ] README in `/work/projects/nn9305k/`
- [ ] `/work/projects/nn9305k/samplefiles/new_user.txt` - follow instructions

**Everybody**
- [ ] read NEWS regularly: `/work/projects/nn9305k/`

#### All **projects directories** should be organized as such:

- [ ] : README copy from README_datafile - in samplefiles directory and filled about -> rename to README_tarfilename
- [ ] directory restructure
  - [ ] rawdata
  - [ ] analysis
  - [ ] scripts
  - [ ] logs
- [ ] please fill in the data registry file : /projects/nn9305k/sys/DataRegistry.txt

[Abel User Guide]:https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/
[Managing Data on Abel]:https://www.uio.no/english/services/it/research/hpc/abel/help/user-guide/data.html
[Workflow]:https://docs.google.com/drawings/d/e/2PACX-1vSY_KCj3fubTH1zk6ZkOL6eLhoOOuAbp4bfu1YkOAvkadHhPfbuZrsepwHCUpEqwr45Zqt2hlEoCwVk/pub?w=960&h=720
[Computer Ressources at CEES]:https://github.com/uio-cees/hpc/wiki/Computer-resources
