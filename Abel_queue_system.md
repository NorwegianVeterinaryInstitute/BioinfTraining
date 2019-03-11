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

$input_file1="/work/projects/nn9305k/samplefiles/Test1.fastq"
$input_file=2"/work/projects/nn9305k/samplefiles/Test2.fastq"

$output_file1="Trimmed_Test1.fastq"
$output_file2="Trimmed_Unpaired_Test1.fastq"
$output_file3="Trimmed_Test2.fastq"
$output_file4="Trimmed_Unpaired_Test2.fastq"


/work/projects/nn9305k/bin/trimmomatic-0.36.jar PE -threads 1 -trimlog vibrio_trimlog.log $input_file1 $input_file2 $output_file1 $output_file2 $output_file4 $output_file4 ILLUMINACLIP:adapter.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:20:15 MINLEN:60
```

# Submit SLURM Job
```
sbatch Trim.sh

squeue -u your_user_name

squeue 

scancel <jobid>

scancel -u jeevka 
```

