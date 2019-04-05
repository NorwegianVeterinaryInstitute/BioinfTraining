**Executing PointFinder**
-------------------------
The tool detects chromosomal mutations predictive of drug resistance based on WGS data.

Contact Jeevan in slack if you have any issues or further assistance (F. ex. run the tool for multiple isolates).

#### User Manual 
https://bitbucket.org/genomicepidemiology/pointfinder

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

conda activate PointFinder

# Database location
PointFinder_DB="/work/projects/nn9305k/src/PointFinder_DB/src/"
python PointFinder.py -p $PointFinder_DB -i input_file -o output_file

conda deactivate 
```
