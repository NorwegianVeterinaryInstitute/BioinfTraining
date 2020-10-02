# Removing contamination

## Sample data location
Saga: /cluster/projects/nn9305k/tutorial/

## Removing PhiX

bbduk.sh threads=5 ref=phix_location,adapter_location in1=INF1 in2=INF2 out=OF1 out2=OF2 k=31 ktrim=r mink=11 hdist=1 tbo=f tpe=f qtrim=r trimq=15 maq=15 minlen=36 forcetrimright=149 stats=stats.txt > log_file

## Removing Human DNA 

conda activate BBTools

IF1=""
IF2=""
OF1=""
OF2=""
bbduk.sh threads=5 ref=phix_location,adapter_location in1=$IF1 in2=$IF2 out=$OF1 out2=$OF2 k=31 ktrim=r mink=11 hdist=1 tbo=f tpe=f qtrim=r trimq=15 maq=15 minlen=36 forcetrimright=149 stats=stats.txt > log_file

conda deactivate

## Other Contamination

screen

srun --account=nn9305k --qos=devel --mem-per-cpu=4800M --cpus-per-task=4 --time=0:30:00 --pty bash -i

conda activate kraken2

kraken2 --db /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312 --threads 4 --output 18-Cjejuni-927.kraken2.out --report 18-Cjejuni-927.kraken2_report.txt --minimum-base-quality 20 --paired --gzip-compressed 18-Cjejuni-927_Subsampled_L008_R1_0089.fastq.gz 18-Cjejuni-927_Subsampled_L008_R2_0089.fastq.gz

kraken2 --db /cluster/shared/biobases/classification_dbs/minikraken_2_8GB_20200312 --threads 4 --output 17-Cjejuni-CCUG11284T.kraken2.out --report 17-Cjejuni-CCUG11284T.kraken2_report.txt --minimum-base-quality 20 --paired --gzip-compressed 17-Cjejuni-CCUG11284T_Subsampled_L008_R1_0087.fastq.gz 17-Cjejuni-CCUG11284T_Subsampled_L008_R2_0087.fastq.gz

conda deactivate
