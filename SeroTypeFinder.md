**Executing SeroTypeFinder**
============================
SerotypeFinder identifies the serotype in total or partial sequenced isolates of E. coli.

Contact Jeevan in slack if you have any issues or further assistance (F. ex. run the tool for multiple isolates).

#### User Manual 
https://bitbucket.org/genomicepidemiology/serotypefinder/src/master/


#### Here is the EXAMPLE SLURM script for Abel to excute the tool.
Important rules to follow
* Refer the user manual for all the parameters in the tool
* Keep your data in /project/nn9305k/
* Store your resutls also in /project/nn9305k/
* Execute the script from your home directory

```
#!/bin/bash
#SBATCH --job-name=DontKillMe
#SBATCH --account=nn9305k
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=32G

## Set up job environment:
source /cluster/bin/jobsetup

conda activate SeroTyperFinder
DB=/work/projects/nn9305k/src/SeroTypeFinder/serotypefinder_db/
python /work/projects/nn9305k/src/SeroTypeFinder/serotypefinder.py -i input_file -o output -p $DB -mp /work/projects/nn9305k/src/kma/kma
conda deactivate 
```
