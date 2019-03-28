**Executing ResFinder** 
----------------------
Use the below code to execute ResFinder Abel. Dont need to mention BLAST location.
Contact Jeevan if you have any issues. 

```
conda activate ResFinder 

# Location of PointFinder DB
PF_DB="/work/projects/nn9305k/src/ResFinder/ResFinderDB/src/"

python /work/projects/nn9305k/src/ResFinder/src/resfinder.py -i <Input File> -p /work/projects/nn9305k/src/ResFinder/ResFinderDB/src/ -k /work/projects/nn9305k/src/kma/ -o Output

conda deactivate
```
