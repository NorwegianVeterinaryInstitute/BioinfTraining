**ShighaTyper**
---------------
ShigaTyper is a quick and easy tool designed to determine Shigella serotype using Illumina paired end reads with low computation requirement.

Contact Jeevan in slack if you have any issues or further assistance (F. ex. run the tool for multiple isolates).

#### User Manual 
https://bitbucket.org/genomicepidemiology/plasmidfinder/src

#### For further reading
https://aem.asm.org/content/85/7/e00165-19

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

conda activate ShigaTyper
python /work/projects/nn9305k/src/ShigaTyper/shigatyper/shigatyper.py Read1 Read2 -n sample_name 
conda deactivate 
```
