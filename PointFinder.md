*Executing PointFinder*
=========================
Location of PointFinder. Contact Jeevan if there is any issues. 



**The following code should be used inside a SLURM script
We are preparing a full help page on how to run using SLURN script.** 

```
/work/projects/nn9305k/src/PointFinder
```

```
conda activate PointFinder

PointFinder_DB="/work/projects/nn9305k/src/PointFinder_DB/src/"
python PointFinder.py -p $PointFinder_DB

conda deactivate 
```
