**Executing VirulanceFinder**
-----------------------------
VirulenceFinder identifies viruelnce genes in total or partial sequenced isolates of bacteria - at the moment only E. coli, Enterococcus, S. aureus and Listeria are available.

Contact Jeevan in slack if you have any issues or further assistance (F. ex. run the tool for multiple isolates).

#### User Manual 
https://bitbucket.org/genomicepidemiology/virulencefinder

#### For further reading
https://www.ncbi.nlm.nih.gov/pubmed/24574290

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

# Activate Conda environment 
conda activate VirulanceFinder
DB="/work/projects/nn9305k/src/VirulanceFinder/Virulance_DB/virulencefinder_db"
python /work/projects/nn9305k/src/VirulanceFinder/src/virulencefinder.py -p $DB 
conda deactivate
```
