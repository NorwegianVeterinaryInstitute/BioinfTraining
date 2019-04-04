**Executing VirulanceFinder**
-----------------------------
Don't need to mention about BLAST. Contact Jeevan if there is any issue.


**The following code should be used inside a SLURM script.**
**We are preparing a full help page on how to run using SLURM script**


```
DB="/work/projects/nn9305k/src/VirulanceFinder/Virulance_DB/virulencefinder_db"

conda activate VirulanceFinder

python /work/projects/nn9305k/src/VirulanceFinder/src/virulencefinder.py -p $DB 

conda deactivate

```
