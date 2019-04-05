**Executing ResFinder** 
----------------------
Use the below code to execute ResFinder Abel. Dont need to mention BLAST location.
Contact Jeevan if you have any issues. 

PlasmidFinder identifies plasmids in total or partial sequenced isolates of bacteria.

Contact Jeevan in slack if you have any issues or further assistance (F. ex. run the tool for multiple isolates).

#### User Manual 
https://bitbucket.org/genomicepidemiology/pointfinder

#### Here is the EXAMPLE SLURM script for Abel to excute the tool.
Important rules to follow
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

module load Miniconda3/4.4.10

conda activate ResFinder 

# Location of PointFinder DB
PF_DB="/work/projects/nn9305k/src/ResFinder/ResFinderDB/src/"

python /work/projects/nn9305k/src/ResFinder/src/resfinder.py -i <Input File> -p /work/projects/nn9305k/src/ResFinder/ResFinderDB/src/ -k /work/projects/nn9305k/src/kma/ -o Output

conda deactivate
```
