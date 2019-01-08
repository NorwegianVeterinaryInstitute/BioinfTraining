|directory|path                |what      |characteristics|
|:--------|:-------------------|:---------|:--------------|
|~        |usit/abel/u1/username|home area |"home-home" area -> when                |
|work     |/work/projects/nn9305k|common resources (eg. softwares, databases)|10 TB - faster disc, closer to nodes connection|
|project   |/projects/nn9305k   |Project data area (backup, raw data..)|30 TB - slower disc|
|   |   |   |computing nodes: scratch ...|
|   |   |   |bifrost and slurm?   |
|   |   |   |   |
|other good to know:|   |   |   |
|/work/databases/bio/ncbi|   |version of NCBI databases on Abel|   |

starting jobs...

cd /work/jobs/1000. number that you have been assigned
or cd /work/users/username <=> $USERWORK  -> faster but do not save data -> removed after 45days

the work area: `/work/projects/nn9305k` will only be used to:
- contain programs and librairies for programs to work
- public reference data, template files, adapters files ...
- the README should and New_user.txt should be checked-modified -> be sure areas ok

- add a "workflow schema" between areas

the project area: `/projects/nn9305k` that contains
- project data, including raw sequencing files (zipped/tared and un-Z/T)
- data analyses that are performed
- this is where the work will NOW be done.-> **so our /work/projects/nn9305k/home/username (incl. all subfolder should be move here)**

- iformation should be rechecked and corrected if necessary:
 - [ ] Bioinformatics course:  working_with_hpc.md
 - [ ] organizational/wiki
 - [ ] annotation.md
 - [ ] data_pre_processing.md
 - [ ] run_mothur.md
 - [ ] specificgene.md
 - [ ] tutorial_kraken.md
 
 - SLURM area: $
