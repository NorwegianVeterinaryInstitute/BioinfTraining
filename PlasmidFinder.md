
**Executing PlasmidFinder**
-------------------------
Use the below code to execute PlasmidFinder Abel. Dont need to mention the BLAST location. 
Contact Jeevan if you have any issues.

**The following code should be used inside a SLURM script.**
**We are preparing a full help page on how to run using SLURM script**

```
conda activate PlasmidFinder

DB="/work/projects/nn9305k/src/PlasmidFinder/PlasmidFinder_DB/plasmidfinder_db/"

python /work/projects/nn9305k/src/PlasmidFinder/plasmidfinder/plasmidfinder.py -p $DB 

conda deactivate PlasmidFinder
```
