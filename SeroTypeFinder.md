**Executing SeroTypeFinder**
============================
https://bitbucket.org/genomicepidemiology/serotypefinder/src/master/ 


**The following code should be used inside a SLURM script.**
**We are preparing a full help page on how to run using SLURM script**


```
conda activate SeroTyperFinder
DB=/work/projects/nn9305k/src/SeroTypeFinder/serotypefinder_db/
cd /work/projects/nn9305k/src/SeroTypeFinder/
python serotypefinder.py --help 
conda deactivate 
```
