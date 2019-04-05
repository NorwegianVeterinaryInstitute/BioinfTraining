**Excuting SeqSero**
--------------------
SeqSero is a pipeline for Salmonella serotype determination from raw sequencing reads or genome assemblies. 

Contact Jeevan in slack if you have any issues or further assistance (F. ex. run the tool for multiple isolates).

#### User Manual 
https://github.com/denglab/SeqSero 

#### For further reading
http://jcm.asm.org/content/early/2015/03/05/JCM.00323-15 

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

conda activate SeqSero_Shared
python /work/projects/nn9305k/src/SeqSero/SeqSero/SeqSero.py -m "1" -i input_file -b "mem"
conda deactivate 
```
