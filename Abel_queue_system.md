# 1. Login to Abel 
```
ssh your_user_name@abel.uio.no
```

# 2. Qlogin 
```
qlogin --account=nn9305k --ntasks-per-node=16
```
# 3. Prepare SLURM Script
```
mkdir QSystem_Test
cp /work/projects/nn9305k/samplefiles/Sample_Slurm_Script.slurm QSystem_Test/Trim.sh
```
Add Trmmomatic commands to SLURM script 

``` 
# For Trimmomataic
module load java

$input_file1=""
$input_file=2""

$output_file1=""
$output_file2=""
$output_file3=""
$output_file4=""


/work/projects/nn9305k/bin/trimmomatic-0.36.jar PE -threads 1 -trimlog vibrio_trimlog.log $input_file1 $input_file2 $output_file1 $output_file2 $output_file4 $output_file4 ILLUMINACLIP:adapter.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:20:15 MINLEN:60
```





