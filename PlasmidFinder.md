**PlasmidFinder**
-------------------------
PlasmidFinder identifies plasmids in total or partial sequenced isolates of bacteria.

Contact Jeevan in slack if you have any issues or further assistance (F. ex. run the tool for multiple isolates).

#### User Manual 
https://bitbucket.org/genomicepidemiology/plasmidfinder/src

#### For further reading
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068535/ 

#### Here is the EXAMPLE SLURM script for Abel to excute the tool.
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
conda activate PlasmidFinder

# Database location
DB="/work/projects/nn9305k/src/PlasmidFinder/PlasmidFinder_DB/plasmidfinder_db/"

# Note: Dont need to mention the BLAST location
python /work/projects/nn9305k/src/PlasmidFinder/plasmidfinder/plasmidfinder.py -p $DB -i input_file -o output_file

# deactivate Conda environment 
conda deactivate
```
